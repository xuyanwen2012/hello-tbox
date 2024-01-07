#ifndef MORTON_H
#define MORTON_H

#include <tbox/tbox.h>

#include "cglm/types.h"

#define MORTON_BITS 31

typedef tb_uint32_t morton_t;
typedef tb_uint16_t coord_t;

static inline morton_t morton3D_SplitBy3bits(const coord_t a) {
  morton_t x = ((morton_t)a) & 0x000003ff;
  x = (x | x << 16) & 0x30000ff;
  x = (x | x << 8) & 0x0300f00f;
  x = (x | x << 4) & 0x30c30c3;
  x = (x | x << 2) & 0x9249249;
  return x;
};

static inline morton_t m3D_e_magicbits(const coord_t x,
                                       const coord_t y,
                                       const coord_t z) {
  return morton3D_SplitBy3bits(x) | (morton3D_SplitBy3bits(y) << 1) |
         (morton3D_SplitBy3bits(z) << 2);
}

static inline coord_t morton3D_GetThirdBits(const morton_t m) {
  morton_t x = m & 0x9249249;
  x = (x ^ (x >> 2)) & 0x30c30c3;
  x = (x ^ (x >> 4)) & 0x0300f00f;
  x = (x ^ (x >> 8)) & 0x30000ff;
  x = (x ^ (x >> 16)) & 0x000003ff;
  return (coord_t)x;
}

static inline void m3D_d_magicbits(const morton_t m, coord_t* xyz) {
  xyz[0] = morton3D_GetThirdBits(m);
  xyz[1] = morton3D_GetThirdBits(m >> 1);
  xyz[2] = morton3D_GetThirdBits(m >> 2);
}

/* //////////////////////////////////////////////////////////////////////////////////////
 * Step 1: Convert 3D point to morton code (32-bit)
 */

static morton_t single_point_to_code_v2(tb_float_t x,
                                        tb_float_t y,
                                        tb_float_t z,
                                        const tb_float_t min_coord,
                                        const tb_float_t range);

void morton32_to_xyz(vec4* ret,
                     const morton_t code,
                     const float min_coord,
                     const float range);

tb_int_t compare_uint32_t(const void* a, const void* b);

void convert_xyz_to_morton_code(const vec4* data,
                                tb_uint32_t* morton_keys,
                                const tb_size_t n,
                                const tb_float_t min_coord,
                                const tb_float_t range);

#endif  // MORTON_H