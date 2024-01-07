#pragma once

#include <stdint.h>

typedef uint32_t morton;
typedef uint16_t coord;

static morton morton3D_SplitBy3bits(const coord a) {
  morton x = ((morton)a) & 0x000003ff;
  x = (x | x << 16) & 0x30000ff;
  x = (x | x << 8) & 0x0300f00f;
  x = (x | x << 4) & 0x30c30c3;
  x = (x | x << 2) & 0x9249249;
  return x;
};

static morton m3D_e_magicbits(const coord x, const coord y, const coord z) {
  return morton3D_SplitBy3bits(x) | (morton3D_SplitBy3bits(y) << 1) |
         (morton3D_SplitBy3bits(z) << 2);
}
