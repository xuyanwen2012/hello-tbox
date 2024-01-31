#include "morton.h"

static morton_t single_point_to_code_v2(tb_float_t x,
                                        tb_float_t y,
                                        tb_float_t z,
                                        const tb_float_t min_coord,
                                        const tb_float_t range) {
  const tb_float_t bit_scale = 1024.0f;

  x = (x - min_coord) / range;
  y = (y - min_coord) / range;
  z = (z - min_coord) / range;

  return m3D_e_magicbits((coord_t)(x * bit_scale),
                         (coord_t)(y * bit_scale),
                         (coord_t)(z * bit_scale));
}

void morton32_to_xyz(vec4* ret,
                     const morton_t code,
                     const tb_float_t min_coord,
                     const tb_float_t range) {
  const tb_float_t bit_scale = 1024.0f;

  coord_t dec_raw_x[3];
  m3D_d_magicbits(code, dec_raw_x);

  const tb_float_t dec_x =
      ((tb_float_t)dec_raw_x[0] / bit_scale) * range + min_coord;
  const tb_float_t dec_y =
      ((tb_float_t)dec_raw_x[1] / bit_scale) * range + min_coord;
  const tb_float_t dec_z =
      ((tb_float_t)dec_raw_x[2] / bit_scale) * range + min_coord;

  // vec4 result = {dec_x, dec_y, dec_z, 1.0f};
  // glm_vec4_copy(result, *ret);
  (*ret)[0] = dec_x;
  (*ret)[1] = dec_y;
  (*ret)[2] = dec_z;
  (*ret)[3] = 1.0f;
}

// functor for uint32_t, used in qsort
tb_int_t compare_uint32_t(const void* a, const void* b) {
  const tb_uint32_t value1 = *(const tb_uint32_t*)a;
  const tb_uint32_t value2 = *(const tb_uint32_t*)b;

  if (value1 < value2) return -1;
  if (value1 > value2) return 1;
  return 0;
}

void convert_xyz_to_morton_code(const vec4* data,
                                tb_uint32_t* morton_keys,
                                const tb_size_t n,
                                const tb_float_t min_coord,
                                const tb_float_t range) {
#pragma omp parallel for schedule(static)
  for (tb_int_t i = 0; i < n; i++) {
    morton_keys[i] = single_point_to_code_v2(
        data[i][0], data[i][1], data[i][2], min_coord, range);
  }
}