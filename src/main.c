/* //////////////////////////////////////////////////////////////////////////////////////
 * includes
 */
#include <cglm/cglm.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "brt.h"
#include "morton.h"
#include "tbox/tbox.h"

/* //////////////////////////////////////////////////////////////////////////////////////
 * Step 1: Convert 3D point to morton code (32-bit)
 */

morton_t single_point_to_code_v2(tb_float_t x,
                                 tb_float_t y,
                                 tb_float_t z,
                                 const tb_float_t min_coord,
                                 const tb_float_t range) {
  x = (x - min_coord) / range;
  y = (y - min_coord) / range;
  z = (z - min_coord) / range;

  return m3D_e_magicbits(
      (coord_t)(x * 1024), (coord_t)(y * 1024), (coord_t)(z * 1024));
}

tb_int_t compare_uint32_t(const void* a, const void* b) {
  const tb_uint32_t value1 = *(const tb_uint32_t*)a;
  const tb_uint32_t value2 = *(const tb_uint32_t*)b;

  if (value1 < value2) return -1;
  if (value1 > value2) return 1;
  return 0;
}

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

/* //////////////////////////////////////////////////////////////////////////////////////
 * main
 */
tb_int_t main(const tb_int_t argc, tb_char_t** argv) {
  if (!tb_init(tb_null, tb_null)) return -1;

  // read n from command line
  tb_int_t n = 1024 * 1024;
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

  // #pragma omp parallel for
  for (tb_int_t i = 0; i < n; i++) {
    data[i][0] = (tb_float_t)tb_rand() / RAND_MAX * range + min_coord;
    data[i][1] = (tb_float_t)tb_rand() / RAND_MAX * range + min_coord;
    data[i][2] = (tb_float_t)tb_rand() / RAND_MAX * range + min_coord;
    data[i][3] = 1.0f;
  }

  // peek 10 points
  for (tb_size_t i = 0; i < 10; i++) {
    printf("data[%lu] = (%f, %f, %f)\n", i, data[i][0], data[i][1], data[i][2]);
  }

#pragma omp parallel for
  for (tb_int_t i = 0; i < n; i++) {
    morton_keys[i] = single_point_to_code_v2(
        data[i][0], data[i][1], data[i][2], min_coord, range);
  }

  qsort(morton_keys, n, sizeof(tb_uint32_t), compare_uint32_t);

  const tb_size_t n_unique_keys = remove_consecutive_duplicates(morton_keys, n);
  // const tb_size_t n_unique_keys = unique(morton_keys, 0, n);
  tb_trace_i("n_unique_keys = %lu", n_unique_keys);
  tb_trace_i("n - n_unique_keys = %lu", n - n_unique_keys);

  radix_tree_t tree;
  init_radix_tree(&tree, morton_keys, n_unique_keys, min_coord, max_coord);

  // peek 10 morton
  for (tb_size_t i = 0; i < 10; i++) {
    printf("morton_keys[%lu] = %u\n", i, morton_keys[i]);
  }

  free_radix_tree(&tree);
  tb_free(data);
  tb_free(morton_keys);

  tb_exit();
  return EXIT_SUCCESS;
}
