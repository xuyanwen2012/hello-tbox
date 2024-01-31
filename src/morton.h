#ifndef MORTON_H
#define MORTON_H

#include <tbox/tbox.h>

#include "cglm/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// Thus the maximum depth of the octree is 10 for 32-bit morton code.
// for 32-bit morton, each 3-bit is used to encode one coordinate.
// we can only use 10 chunk of 3-bits, so 2 bits are wasted.
// for 64-bit morton,
// we can use 21 chunk of 3-bits, so 63. 1 bit is wasted.
enum { MORTON_BITS = 30 };

typedef tb_uint32_t morton_t;
typedef tb_uint16_t coord_t;

static morton_t morton3D_SplitBy3bits(const coord_t a) {
  morton_t x = ((morton_t)a) & 0x000003ff;
  x = (x | x << 16) & 0x30000ff;
  x = (x | x << 8) & 0x0300f00f;
  x = (x | x << 4) & 0x30c30c3;
  x = (x | x << 2) & 0x9249249;
  return x;
}

static morton_t m3D_e_magicbits(const coord_t x,
                                const coord_t y,
                                const coord_t z) {
  return morton3D_SplitBy3bits(x) | (morton3D_SplitBy3bits(y) << 1) |
         (morton3D_SplitBy3bits(z) << 2);
}

static coord_t morton3D_GetThirdBits(const morton_t m) {
  morton_t x = m & 0x9249249;
  x = (x ^ (x >> 2)) & 0x30c30c3;
  x = (x ^ (x >> 4)) & 0x0300f00f;
  x = (x ^ (x >> 8)) & 0x30000ff;
  x = (x ^ (x >> 16)) & 0x000003ff;
  return x;
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
                                        tb_float_t min_coord,
                                        tb_float_t range);

void morton32_to_xyz(vec4* ret,
                     morton_t code,
                     tb_float_t min_coord,
                     tb_float_t range);

tb_int_t compare_uint32_t(const void* a, const void* b);

void convert_xyz_to_morton_code(const vec4* data,
                                tb_uint32_t* morton_keys,
                                tb_size_t n,
                                tb_float_t min_coord,
                                tb_float_t range);

#ifdef __cplusplus
}
#endif

#endif  // MORTON_H