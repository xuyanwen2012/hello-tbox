#ifndef RADIX_TREE_H
#define RADIX_TREE_H

#include <tbox/tbox.h>

#include "morton.h"

/* //////////////////////////////////////////////////////////////////////////
 * Types
 */

typedef struct brt_nodes {
  const morton_t* morton_codes;

  // Flags determining whether the left and right children are leaves
  tb_bool_t* hasLeafLeft;
  tb_bool_t* hasLeafRight;

  // The number of bits in the mortonCode this node represents
  // Corresponds to delta_node in [Karras]
  tb_uint8_t* prefixN;

  // Index of left child of this node
  // Right child is leftChild + 1
  tb_int_t* leftChild;

  // Index of parent
  tb_int_t* parent;
} brt_nodes_t;

// initialization and memory allocation
void init_brt_nodes(brt_nodes_t* nodes,
                    const morton_t* morton_codes,
                    tb_int_t n_nodes);

void free_brt_nodes(const brt_nodes_t* nodes);

typedef struct {
  brt_nodes_t d_tree;  // radix tree on GPU
  tb_int_t n_pts;      // number of points (n = number of Unique morton codes)
  tb_int_t n_nodes;    // number of radix tree nodes (n - 1)

  // minimum and maximum coordinate in points.
  // Represents the scaling factor for the morton codes
  tb_float_t min_coord;
  tb_float_t max_coord;
} radix_tree_t;

void init_radix_tree(radix_tree_t* tree,
                     const morton_t* sorted_morton_keys,
                     tb_int_t n_unique_keys,
                     tb_float_t min_coord,
                     tb_float_t max_coord);

void build_radix_tree(radix_tree_t* tree);

void free_radix_tree(const radix_tree_t* tree);

/* //////////////////////////////////////////////////////////////////////////
 */

tb_int_t log2_ceil_u32(tb_uint32_t x);

tb_uint8_t delta_u32(tb_uint32_t a, tb_uint32_t b);

#endif  // RADIX_TREE_H