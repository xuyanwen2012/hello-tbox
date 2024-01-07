/* //////////////////////////////////////////////////////////////////////////////////////
 * includes
 */
#include <cglm/cglm.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "morton.h"
#include "tbox/tbox.h"

morton single_point_to_code_v2(tb_float_t x,
                               tb_float_t y,
                               tb_float_t z,
                               const tb_float_t min_coord,
                               const tb_float_t range) {
  x = (x - min_coord) / range;
  y = (y - min_coord) / range;
  z = (z - min_coord) / range;

  return m3D_e_magicbits(
      (coord)(x * 1024), (coord)(y * 1024), (coord)(z * 1024));
}

static void functionToMeasure() {}

tb_int_t compare_uint32_t(const void* a, const void* b) {
  tb_uint32_t value1 = *((const tb_uint32_t*)a);
  tb_uint32_t value2 = *((const tb_uint32_t*)b);

  if (value1 < value2)
    return -1;
  else if (value1 > value2)
    return 1;
  else
    return 0;
}

/* //////////////////////////////////////////////////////////////////////////////////////
 * main
 */
tb_int_t main(tb_int_t argc, tb_char_t** argv) {
  if (!tb_init(tb_null, tb_null)) return -1;

  // read n from command line
  tb_int_t n = 1024 * 1024;
  if (argc > 1) {
    n = tb_atoi(argv[1]);
  }
  printf("n = %d\n", n);

  tb_int_t num_threads = 16;
  if (argc > 2) {
    num_threads = atoi(argv[2]);
  }
  printf("nthreads = %d\n", num_threads);

  omp_set_num_threads(num_threads);

  srand(114514);

  // allocate n vec4 elements
  vec4* data = tb_nalloc_type(n, vec4);
  tb_assert_and_check_return_val(data, EXIT_FAILURE);

  const float min_coord = 0.0f;
  const float max_coord = 1024.0f;
  const float range = max_coord - min_coord;

  // #pragma omp parallel for
  for (tb_int_t i = 0; i < n; i++) {
    data[i][0] = (float)rand() / RAND_MAX * range + min_coord;
    data[i][1] = (float)rand() / RAND_MAX * range + min_coord;
    data[i][2] = (float)rand() / RAND_MAX * range + min_coord;
    data[i][3] = 1.0f;
  }

  // peek 10 points
  for (tb_size_t i = 0; i < 10; i++) {
    printf("data[%lu] = (%f, %f, %f)\n", i, data[i][0], data[i][1], data[i][2]);
  }

  // allocate n morton code
  tb_uint32_t* morton_keys = tb_nalloc_type(n, tb_uint32_t);
  tb_assert_and_check_return_val(morton_keys, EXIT_FAILURE);

#pragma omp parallel for
  for (tb_int_t i = 0; i < n; i++) {
    morton_keys[i] = single_point_to_code_v2(
        data[i][0], data[i][1], data[i][2], min_coord, range);
  }

  qsort(morton_keys, n, sizeof(tb_uint32_t), compare_uint32_t);

  // peek 10 morton
  for (tb_size_t i = 0; i < 10; i++) {
    printf("morton_keys[%lu] = %u\n", i, morton_keys[i]);
  }

  tb_free(data);
  tb_free(morton_keys);

  tb_exit();
  return EXIT_SUCCESS;
}
