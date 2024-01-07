#ifndef BRT_H
#define BRT_H

#include <tbox/tbox.h>

#include "morton.h"

typedef struct brt_nodes {
  const morton_t* morton_codes;

  // Flags determining whether the left and right children are leaves
  tb_bool_t* hasLeafLeft;
  tb_bool_t* hasLeafRight;

  // The number of bits in the mortonCode this node represents
  // Corresponds to delta_node in [Karras]
  uint8_t* prefixN;

  // Index of left child of this node
  // Right child is leftChild + 1
  int* leftChild;

  // Index of parent
  int* parent;
} brt_nodes_t;

void init_brt_nodes(brt_nodes_t* nodes,
                    const morton_t* morton_codes,
                    int n_nodes);

void free_brt_nodes(brt_nodes_t* nodes);

typedef struct radix_tree {
  brt_nodes_t d_tree;  // radix tree on GPU
  int n_pts;           // number of points (n = number of Unique morton codes)
  int n_nodes;         // number of radix tree nodes (n - 1)

  // minimum and maximum coordinate in points.
  // Represents the scaling factor for the morton codes
  float min_coord;
  float max_coord;
} radix_tree_t;

void init_radix_tree(radix_tree_t* tree,
                     const morton_t* sorted_mortons_keys,
                     int n_unique_keys,
                     float min_coord,
                     float max_coord);

void free_radix_tree(radix_tree_t* tree);

#endif  // BRT_H