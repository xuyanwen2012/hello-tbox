/* //////////////////////////////////////////////////////////////////////////////////////
 * includes
 */
#include <cglm/cglm.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "cglm/types.h"
#include "morton.h"
#include "octree.h"
#include "radix_tree.h"
#include "tbox/prefix/type.h"
#include "tbox/tbox.h"

/* //////////////////////////////////////////////////////////////////////////////////////
 * Step 3: Remove consecutive duplicates in morton
 */
static tb_size_t remove_consecutive_duplicates(tb_uint32_t* array,
                                               tb_size_t size) {
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
static tb_size_t std_unique(tb_uint32_t* array,
                            tb_ptrdiff_t first,
                            const tb_ptrdiff_t last) {
  if (first == last) return last;

  tb_ptrdiff_t result = first;
  while (++first != last) {
    if (!(array[result] == array[first]) && ++result != first) {
      array[result] = array[first];
    }
  }

  return ++result;
}

static void create_radix_tree(radix_tree_t* tree,
                              const tb_uint32_t* morton_keys,
                              const tb_size_t n_unique_keys,
                              const tb_float_t min_coord,
                              const tb_float_t max_coord) {
  tree->n_pts = (tb_int_t)n_unique_keys;
  tree->n_nodes = (tb_int_t)n_unique_keys - 1;
  tree->min_coord = min_coord;
  tree->max_coord = max_coord;
  tree->d_tree.morton_codes = morton_keys;
  tree->d_tree.hasLeafLeft = tb_nalloc_type(n_unique_keys, tb_bool_t);
  tree->d_tree.hasLeafRight = tb_nalloc_type(n_unique_keys, tb_bool_t);
  tree->d_tree.prefixN = tb_nalloc_type(n_unique_keys, tb_uint8_t);
  tree->d_tree.leftChild = tb_nalloc_type(n_unique_keys, tb_int_t);
  tree->d_tree.parent = tb_nalloc_type(n_unique_keys, tb_int_t);

  tb_assert_and_check_return(tree->d_tree.hasLeafLeft);
  tb_assert_and_check_return(tree->d_tree.hasLeafRight);
  tb_assert_and_check_return(tree->d_tree.prefixN);
  tb_assert_and_check_return(tree->d_tree.leftChild);
  tb_assert_and_check_return(tree->d_tree.parent);
}

static void destroy_radix_tree(const radix_tree_t* tree) {
  tb_free(tree->d_tree.hasLeafLeft);
  tb_free(tree->d_tree.hasLeafRight);
  tb_free(tree->d_tree.prefixN);
  tb_free(tree->d_tree.leftChild);
  tb_free(tree->d_tree.parent);
}

// https://en.cppreference.com/w/cpp/algorithm/partial_sum#Version_1
static tb_int_t* std_partial_sum(const tb_int_t* data,
                                 tb_ptrdiff_t first,
                                 const tb_ptrdiff_t last,
                                 tb_int_t* d_first) {
  if (first == last) return d_first;
  tb_int_t sum = data[first];
  *d_first = sum;
  while (++first != last) {
    sum += data[first];
    *++d_first = sum;
  }
  return ++d_first;
}

void print_bits(const uint32_t x) {
  for (int i = 31; i >= 0; i--) {
    uint32_t bit = (x >> i) & 1;
    // printf("%u", bit);
    putchar(bit ? '1' : '0');
    if (i % 3 == 0) {
      putchar(' ');  // Add space every 3 bits
    }
  }
}

tb_int_t extract_sub_index(const morton_t code, const tb_int_t level) {
  const tb_int_t bit_position = level * 3;  // 30 -
  return (code >> bit_position) & 0x07;     // 0b111
}

void print_octree_path(uint32_t x) {
  for (int i = 0; i < 10; i++) {
    tb_int_t sub_index = extract_sub_index(x, i);
    printf("%d ", sub_index);
  }
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
  printf("Total memory used: %f MB\n", (tb_float_t)total_memory / 1024 / 1024);

  const tb_float_t min_coord = 0.0f;
  const tb_float_t max_coord = 1024.0f;
  const tb_float_t range = max_coord - min_coord;

  for (tb_int_t i = 0; i < n; i++) {
    data[i][0] =
        (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min_coord;
    data[i][1] =
        (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min_coord;
    data[i][2] =
        (tb_float_t)tb_rand() / (tb_float_t)RAND_MAX * range + min_coord;
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
    printf("morton_keys[%lu] =\t%u\t", i, morton_keys[i]);
    print_bits(morton_keys[i]);
    putchar('\t');
    print_octree_path(morton_keys[i]);
    putchar('\n');
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
  for (tb_int_t i = 0; i < 32; ++i) {
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
  for (tb_int_t i = 0; i < 32; ++i) {
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
         (tb_float_t)total_memory / 1024 / 1024);

  // do prefix sum (std version)
  std_partial_sum(edge_count, 0, tree->n_nodes, edge_count_prefix_sum + 1);
  edge_count_prefix_sum[0] = 0;

  // peek 32 prefix sum
  for (tb_int_t i = 0; i < 32; ++i) {
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
         (tb_float_t)total_memory / 1024 / 1024);

  // step 6: make octree

  const tb_float_t tree_range = max_coord - min_coord;
  const tb_int_t root_level = tree->d_tree.prefixN[0] / 3;
  const morton_t root_prefix =
      tree->d_tree.morton_codes[0] >> (MORTON_BITS - (3 * root_level));

  morton32_to_xyz(&octree_nodes[0].corner,
                  root_prefix << (MORTON_BITS - (3 * root_level)),
                  min_coord,
                  tree_range);
  octree_nodes[0].cell_size = tree_range;

  make_oct_nodes(octree_nodes,
                 edge_count_prefix_sum,
                 edge_count,
                 tree->d_tree.morton_codes,
                 tree->d_tree.prefixN,
                 tree->d_tree.parent,
                 min_coord,
                 tree_range,
                 tree->n_nodes);

  // // peek 32 octree nodes
  // for (tb_int_t i = 0; i < 32; ++i) {
  //   printf("octree_nodes[%d].corner = (%f, %f, %f)\n",
  //          i,
  //          octree_nodes[i].corner[0],
  //          octree_nodes[i].corner[1],
  //          octree_nodes[i].corner[2]);
  // }

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
