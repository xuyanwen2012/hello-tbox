#ifndef BASELINE_H
#define BASELINE_H

#include <tbox/tbox.h>

#include "cglm/types.h"
#include "morton.h"

struct Baseline {
  vec4* u_input;
  morton_t* u_morton;
  morton_t* u_morton_alt;
};

#endif  // BASELINE_H