#include "octree.h"

void set_child(oct_node_t* node,
               const tb_int_t child,
               const tb_int_t my_child_idx) {
  node->children[my_child_idx] = child;
  node->child_node_mask |= (1 << my_child_idx);
}

void set_child_node(oct_node_t* node,
                    const tb_int_t leaf,
                    const tb_int_t my_child_idx) {
  node->children[my_child_idx] = leaf;
  node->child_node_mask &= ~(1 << my_child_idx);
}

void count_edges(const tb_uint8_t* prefixN,
                 const tb_int_t* parents,
                 tb_int_t* rt_edge_counts,
                 tb_int_t n_brt_nodes) {
#pragma omp parallel for schedule(static)
  for (tb_size_t i = 1; i < n_brt_nodes; ++i) {
    tb_int_t my_depth = prefixN[i] / 3;
    tb_int_t parent_depth = prefixN[parents[i]] / 3;
    rt_edge_counts[i] = my_depth - parent_depth;
  }
}
