/* //////////////////////////////////////////////////////////////////////////////////////
 * includes
 */
#include <assert.h>
#include <cglm/cglm.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cglm/types.h"
#include "cglm/vec4.h"
#include "morton.h"
#include "octree.h"
#include "radix_tree.h"
#include "tbox/prefix/type.h"
#include "tbox/tbox.h"

// /*
// //////////////////////////////////////////////////////////////////////////////////////
//  * Step 1: Convert 3D point to morton code (32-bit)
//  */

// static morton_t single_point_to_code_v2(tb_float_t x,
//                                         tb_float_t y,
//                                         tb_float_t z,
//                                         const tb_float_t min_coord,
//                                         const tb_float_t range) {
//   // const uint32_t bitscale = 0xFFFFFFFFu >> (32 - (MORTON_BITS / 3));  //
//   1023

//   const uint32_t bitscale = 1024;

//   x = (x - min_coord) / range;
//   y = (y - min_coord) / range;
//   z = (z - min_coord) / range;

//   return m3D_e_magicbits((coord_t)(x * bitscale),
//                          (coord_t)(y * bitscale),
//                          (coord_t)(z * bitscale));
// }

// static void morton32_to_xyz(vec4* ret,
//                             const morton_t code,
//                             const float min_coord,
//                             const float range) {
//   // const uint32_t bitscale = 0xFFFFFFFFu >> (32 - (MORTON_BITS / 3));  //
//   1023 const uint32_t bitscale = 1024;

//   coord_t dec_raw_x[3];
//   // libmorton::morton3D_64_decode(code, dec_raw_x, dec_raw_y, dec_raw_z);
//   m3D_d_magicbits(code, dec_raw_x);

//   float dec_x = ((float)dec_raw_x[0] / bitscale) * range + min_coord;
//   float dec_y = ((float)dec_raw_x[1] / bitscale) * range + min_coord;
//   float dec_z = ((float)dec_raw_x[2] / bitscale) * range + min_coord;

//   // vec4 result = {dec_x, dec_y, dec_z, 1.0f};
//   // glm_vec4_copy(result, *ret);
//   (*ret)[0] = dec_x;
//   (*ret)[1] = dec_y;
//   (*ret)[2] = dec_z;
//   (*ret)[3] = 1.0f;
// }

// // functor for uint32_t, used in qsort
// tb_int_t compare_uint32_t(const void* a, const void* b) {
//   const tb_uint32_t value1 = *(const tb_uint32_t*)a;
//   const tb_uint32_t value2 = *(const tb_uint32_t*)b;

//   if (value1 < value2) return -1;
//   if (value1 > value2) return 1;
//   return 0;
// }

// void convert_xyz_to_morton_code(const vec4* data,
//                                 tb_uint32_t* morton_keys,
//                                 const tb_size_t n,
//                                 const tb_float_t min_coord,
//                                 const tb_float_t range) {
// #pragma omp parallel for
//   for (tb_int_t i = 0; i < n; i++) {
//     morton_keys[i] = single_point_to_code_v2(
//         data[i][0], data[i][1], data[i][2], min_coord, range);
//   }
// }

/* //////////////////////////////////////////////////////////////////////////////////////
 * Step 3: Remove consecutive duplicates in morton
 */
tb_size_t remove_consecutive_duplicates(tb_uint32_t* array, tb_size_t size) {
  if (size == 0) return 0;

  tb_size_t index = 0;
  for (tb_size_t i = 1; i < size; ++i) {
    if (array[index] != array[i]) {
      ++index;
      array[index] = array[i];
    }
  }
  return index + 1;
}

// https://en.cppreference.com/w/cpp/algorithm/unique
tb_size_t unique(tb_uint32_t* array, tb_ptrdiff_t first, tb_ptrdiff_t last) {
  if (first == last) return last;

  tb_ptrdiff_t result = first;
  while (++first != last) {
    if (!(array[result] == array[first]) && ++result != first) {
      array[result] = array[first];
    }
  }

  return ++result;
}

void create_radix_tree(radix_tree_t* tree,
                       const tb_uint32_t* morton_keys,
                       const tb_size_t n_unique_keys,
                       const tb_float_t min_coord,
                       const tb_float_t max_coord) {
  tree->n_pts = n_unique_keys;
  tree->n_nodes = n_unique_keys - 1;
  tree->min_coord = min_coord;
  tree->max_coord = max_coord;
  tree->d_tree.morton_codes = morton_keys;
  tree->d_tree.hasLeafLeft = tb_nalloc_type(n_unique_keys, tb_bool_t);
  tree->d_tree.hasLeafRight = tb_nalloc_type(n_unique_keys, tb_bool_t);
  tree->d_tree.prefixN = tb_nalloc_type(n_unique_keys, tb_uint8_t);
  tree->d_tree.leftChild = tb_nalloc_type(n_unique_keys, int);
  tree->d_tree.parent = tb_nalloc_type(n_unique_keys, int);

  tb_assert_and_check_return(tree->d_tree.hasLeafLeft);
  tb_assert_and_check_return(tree->d_tree.hasLeafRight);
  tb_assert_and_check_return(tree->d_tree.prefixN);
  tb_assert_and_check_return(tree->d_tree.leftChild);
  tb_assert_and_check_return(tree->d_tree.parent);
}

void destroy_radix_tree(radix_tree_t* tree) {
  tb_free(tree->d_tree.hasLeafLeft);
  tb_free(tree->d_tree.hasLeafRight);
  tb_free(tree->d_tree.prefixN);
  tb_free(tree->d_tree.leftChild);
  tb_free(tree->d_tree.parent);
}

/* //////////////////////////////////////////////////////////////////////////////////////
 * main
 */
tb_int_t main(const tb_int_t argc, tb_char_t** argv) {
  if (!tb_init(tb_null, tb_null)) return -1;

  // read n from command line
  tb_int_t n = 640 * 480;
  if (argc > 1) {
    n = tb_atoi(argv[1]);
  }

  // printf("n = %d\n", n);
  tb_trace_i("n = %d", n);

  tb_int_t num_threads = 16;
  if (argc > 2) {
    num_threads = tb_atoi(argv[2]);
  }
  // printf("num_threads = %d\n", num_threads);
  tb_trace_i("num_threads = %d", num_threads);

  omp_set_num_threads(num_threads);

  tb_srand(114514);

  // allocate n vec4 elements
  vec4* data = tb_nalloc_type(n, vec4);
  tb_assert_and_check_return_val(data, EXIT_FAILURE);

  // allocate n morton code
  tb_uint32_t* morton_keys = tb_nalloc_type(n, tb_uint32_t);
  tb_assert_and_check_return_val(morton_keys, EXIT_FAILURE);

  // print total memory used
  tb_size_t total_memory = 0;
  total_memory += n * sizeof(vec4);
  total_memory += n * sizeof(tb_uint32_t);
  printf("Total memory used: %f MB\n", (float)total_memory / 1024 / 1024);

  const tb_float_t min_coord = 0.0f;
  const tb_float_t max_coord = 1024.0f;
  const tb_float_t range = max_coord - min_coord;

  for (tb_int_t i = 0; i < n; i++) {
    data[i][0] = (tb_float_t)tb_rand() / RAND_MAX * range + min_coord;
    data[i][1] = (tb_float_t)tb_rand() / RAND_MAX * range + min_coord;
    data[i][2] = (tb_float_t)tb_rand() / RAND_MAX * range + min_coord;
    data[i][3] = 1.0f;
  }

  // peek 32 points
  for (tb_size_t i = 0; i < 32; i++) {
    printf("data[%lu] = (%f, %f, %f)\n", i, data[i][0], data[i][1], data[i][2]);
  }

  // step 1: convert 3D points to morton code
  //
  convert_xyz_to_morton_code(data, morton_keys, n, min_coord, range);

  // step 2: sort morton code
  //
  qsort(morton_keys, n, sizeof(tb_uint32_t), compare_uint32_t);

  // peek 32 morton
  for (tb_size_t i = 0; i < 32; i++) {
    printf("morton_keys[%lu] = %u\n", i, morton_keys[i]);
  }

  // step 3: remove consecutive duplicates
  const tb_size_t n_unique_keys = remove_consecutive_duplicates(morton_keys, n);

  tb_trace_i("n_unique_keys = %lu", n_unique_keys);
  tb_trace_i("n - n_unique_keys = %lu", n - n_unique_keys);

  // step 4: build radix tree (allocate memory)
  radix_tree_t* tree = tb_nalloc_type(1, radix_tree_t);
  tb_assert_and_check_return_val(tree, EXIT_FAILURE);

  create_radix_tree(tree, morton_keys, n_unique_keys, min_coord, max_coord);
  build_radix_tree(tree);

  // peek 32 nodes
  for (int i = 0; i < 32; ++i) {
    printf(
        "idx = %d, code = %u, prefixN = %d, left = %d, parent = %d, "
        "leftLeaf=%d, rightLeft=%d\n",
        i,
        morton_keys[i],
        tree->d_tree.prefixN[i],
        tree->d_tree.leftChild[i],
        tree->d_tree.parent[i],
        tree->d_tree.hasLeafLeft[i],
        tree->d_tree.hasLeafRight[i]);
  }

  // allocate temporary memory for octree

  tb_int_t* edge_count = tb_nalloc_type(tree->n_nodes, tb_int_t);
  tb_assert_and_check_return_val(edge_count, EXIT_FAILURE);

  // step 5: count edges
  count_edges(
      tree->d_tree.prefixN, tree->d_tree.parent, edge_count, tree->n_nodes);

  // peek 32 edges
  for (int i = 0; i < 32; ++i) {
    printf("edge_count[%d] = %d\n", i, edge_count[i]);
  }

  // step 5.5: compute the prefix sum
  // allocate 1 extra element for the last element
  tb_int_t* edge_count_prefix_sum = tb_nalloc_type(tree->n_nodes + 1, tb_int_t);
  tb_assert_and_check_return_val(edge_count_prefix_sum, EXIT_FAILURE);

  // print how much temporary memory is used
  total_memory = 0;
  total_memory += tree->n_nodes * sizeof(tb_int_t);
  total_memory += tree->n_nodes * sizeof(tb_int_t);
  printf("Total memory used for temporary memory: %f MB\n",
         (float)total_memory / 1024 / 1024);

  tb_int_t sum = 0;
  // allocate 1 extra element for the last element
  for (tb_int_t i = 0; i < tree->n_nodes + 1; ++i) {
    edge_count_prefix_sum[i] = sum;
    sum += edge_count[i];
  }

  // peek 32 prefix sum
  for (int i = 0; i < 32; ++i) {
    printf("edge_count_prefix_sum[%d] = %d\n", i, edge_count_prefix_sum[i]);
  }

  const tb_int_t num_oct_nodes = edge_count_prefix_sum[tree->n_nodes];
  tb_trace_i("num_oct_nodes = %d", num_oct_nodes);

  // allocate memory for octree
  oct_node_t* octree_nodes = tb_nalloc_type(num_oct_nodes, oct_node_t);
  tb_assert_and_check_return_val(octree_nodes, EXIT_FAILURE);

  // print how much memory is used
  total_memory = 0;
  total_memory += num_oct_nodes * sizeof(oct_node_t);
  printf("Total memory used for octree: %f MB\n",
         (float)total_memory / 1024 / 1024);

  // step 6: make octree

  float tree_range = max_coord - min_coord;
  int root_level = tree->d_tree.prefixN[0] / 3;
  morton_t root_prefix =
      tree->d_tree.morton_codes[0] >> (MORTON_BITS - (3 * root_level));

  morton32_to_xyz(&octree_nodes[0].corner,
                  root_prefix << (MORTON_BITS - (3 * root_level)),
                  min_coord,
                  tree_range);
  octree_nodes[0].cell_size = tree_range;

  // free temporary memory
  tb_free(edge_count);
  tb_free(edge_count_prefix_sum);

  // free memory
  destroy_radix_tree(tree);
  tb_free(tree);
  tb_free(data);
  tb_free(morton_keys);

  // morton_t test_code =
  //     single_point_to_code_v2(12.0f, 34.0f, 56.0f, min_coord, tree_range);
  // printf("test_code = %u\n", test_code);

  // vec4 result;
  // morton32_to_xyz(&result, test_code, min_coord, tree_range);
  // printf("result = (%f, %f, %f)\n", result[0], result[1], result[2]);

  tb_exit();
  return EXIT_SUCCESS;
}
