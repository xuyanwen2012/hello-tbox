#ifndef MORTON_H
#define MORTON_H

#include <stdint.h>

#define MORTON_BITS 31

typedef uint32_t morton_t;
typedef uint16_t coord_t;

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

#endif  // MORTON_H