#include "octree.h"

#include "morton.h"

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
                 tb_int_t* edge_count,
                 tb_int_t n_brt_nodes) {
  // root
  edge_count[0] = 0;

  // root has no parent, so don't do for index 0
#pragma omp parallel for schedule(static)
  for (tb_size_t i = 1; i < n_brt_nodes; ++i) {
    tb_int_t my_depth = prefixN[i] / 3;
    tb_int_t parent_depth = prefixN[parents[i]] / 3;
    edge_count[i] = my_depth - parent_depth;
  }
}

void make_oct_nodes(oct_node_t* oct_nodes,
                    const tb_int_t* node_offsets,    // prefix sum
                    const tb_int_t* rt_node_counts,  // edge count
                    const morton_t* codes,
                    const tb_uint8_t* rt_prefixN,
                    const tb_int_t* rt_parents,
                    const tb_float_t min_coord,
                    const tb_float_t range,
                    const tb_int_t N  // number of brt nodes
) {
  // the root doesn't represent level 0 of the "entire" octree
  for (int i = 1; i < N; ++i) {
    int root_level = rt_prefixN[0] / 3;
    int oct_idx = node_offsets[i];
    // int n_new_nodes = node_offsets[i] - node_offsets[i - 1];
    int n_new_nodes = rt_node_counts[i];
    for (int j = 0; j < n_new_nodes - 1; ++j) {
      int level = rt_prefixN[i] / 3 - j;
      morton_t node_prefix = codes[i] >> (MORTON_BITS - (3 * level));
      int child_idx = node_prefix & 0b111;
      int parent = oct_idx + 1;

      //   oct_nodes[parent].setChild(oct_idx, child_idx);
      set_child(&oct_nodes[parent], oct_idx, child_idx);

      // calculate corner point
      //   (less significant bits have already been shifted off)
      //   oct_nodes[oct_idx].corner = codeToPoint(
      //   node_prefix << (CODE_LEN - (3 * level)), min_coord, range);

      // each cell is half the size of the level above it
      oct_nodes[oct_idx].cell_size = range / (float)(1 << (level - root_level));
      oct_idx = parent;
    }
  }
}