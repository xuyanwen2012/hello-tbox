#include "morton.h"

static morton_t single_point_to_code_v2(tb_float_t x,
                                        tb_float_t y,
                                        tb_float_t z,
                                        const tb_float_t min_coord,
                                        const tb_float_t range) {
  // const uint32_t bitscale = 0xFFFFFFFFu >> (32 - (MORTON_BITS / 3));  // 1023

  const tb_uint32_t bitscale = 1024;

  x = (x - min_coord) / range;
  y = (y - min_coord) / range;
  z = (z - min_coord) / range;

  return m3D_e_magicbits((coord_t)(x * bitscale),
                         (coord_t)(y * bitscale),
                         (coord_t)(z * bitscale));
}

void morton32_to_xyz(vec4* ret,
                     const morton_t code,
                     const float min_coord,
                     const float range) {
  // const uint32_t bitscale = 0xFFFFFFFFu >> (32 - (MORTON_BITS / 3));  // 1023
  const tb_uint32_t bitscale = 1024;

  coord_t dec_raw_x[3];
  // libmorton::morton3D_64_decode(code, dec_raw_x, dec_raw_y, dec_raw_z);
  m3D_d_magicbits(code, dec_raw_x);

  float dec_x = ((float)dec_raw_x[0] / bitscale) * range + min_coord;
  float dec_y = ((float)dec_raw_x[1] / bitscale) * range + min_coord;
  float dec_z = ((float)dec_raw_x[2] / bitscale) * range + min_coord;

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
#pragma omp parallel for
  for (tb_int_t i = 0; i < n; i++) {
    morton_keys[i] = single_point_to_code_v2(
        data[i][0], data[i][1], data[i][2], min_coord, range);
  }
}