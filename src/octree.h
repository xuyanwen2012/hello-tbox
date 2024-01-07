#ifndef OCTREE_H
#define OCTREE_H

#include <tbox/tbox.h>

#include "cglm/types.h"
#include "morton.h"
#include "tbox/prefix/type.h"

/* //////////////////////////////////////////////////////////////////////////
 * Types
 */

typedef struct {
  // TODO: This is overkill number of pointers
  tb_int_t children[8];

  // bounding box of cell is defined by minimum (x, y, z) coordinate of cell and
  // depth in tree.
  vec4 corner;
  tb_float_t cell_size;

  // For bit position i (from the right):
  //     If 1, children[i] is the index of a child octree node
  //     If 0, the ith child is either absent, or children[i] is the index of a
  //     leaf.
  tb_int_t child_node_mask;

  // For bit position i (from the right):
  //     If 1, children[i] is the index of a leaf (in the corresponding points
  //     array) If 0, the ith child is either absent, or an octree node.
  tb_int_t child_leaf_mask;
} oct_node_t;

// Set a child
//     child: index of octree node that will become the child
//     my_child_idx: which of my children it will be [0-7]
void set_child(oct_node_t* node, tb_int_t child, tb_int_t my_child_idx);

// Set a leaf child
//     leaf: index of point that will become the leaf child
//     my_child_idx; which of my children it will be [0-7]
void set_leaf(oct_node_t* node, tb_int_t leaf, tb_int_t my_child_idx);

/* //////////////////////////////////////////////////////////////////////////
 * Step 5: edge counting
 */

void count_edges(const tb_uint8_t* prefixN,
                 const tb_int_t* parents,
                 tb_int_t* edge_count,
                 tb_int_t n_brt_nodes);

void make_oct_nodes(oct_node_t* oct_nodes,
                    const tb_int_t* node_offsets,    // prefix sum
                    const tb_int_t* rt_node_counts,  // edge count
                    const morton_t* codes,
                    const tb_uint8_t* rt_prefixN,
                    const tb_int_t* rt_parents,
                    const tb_float_t min_coord,
                    const tb_float_t range,
                    const tb_int_t N  // number of brt nodes
);

#endif  // OCTREE_H
