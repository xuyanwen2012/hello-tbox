#include "octree.h"

#include "cglm/vec4.h"
#include "morton.h"

// basically, turn on the bit at the specific 'my_child_idx', and save the
// child's node id to the children array
void set_child(oct_node_t* node,
               const tb_int_t child,
               const tb_int_t my_child_idx) {
  node->children[my_child_idx] = child;
  node->child_node_mask |= (1 << my_child_idx);
}

// basically, turn off the bit at the specific 'my_child_idx', and save the
// child's leaf id to the children array
void set_child_node(oct_node_t* node,
                    const tb_int_t leaf,
                    const tb_int_t my_child_idx) {
  node->children[my_child_idx] = leaf;
  node->child_node_mask &= ~(1 << my_child_idx);
}

// count the number of edges in the radix tree
void count_edges(const tb_uint8_t* prefixN,
                 const tb_int_t* parents,
                 tb_int_t* edge_count,
                 tb_int_t n_brt_nodes) {
  // root
  edge_count[0] = 0;

  // root has no parent, so don't do for index 0
#pragma omp parallel for schedule(static)
  for (tb_int_t i = 1; i < n_brt_nodes; ++i) {
    const tb_int_t my_depth = prefixN[i] / 3;
    const tb_int_t parent_depth = prefixN[parents[i]] / 3;
    edge_count[i] = my_depth - parent_depth;
  }
}

// tb_int_t extract_sub_index(const morton_t code, const tb_int_t level) {
//   const tb_int_t bit_position = level * 3;
//   return (code >> bit_position) & 0x07;  // 0b111
// }

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
  const tb_int_t root_level = rt_prefixN[0] / 3;

  // the root doesn't represent level 0 of the "entire" octree
  for (tb_int_t i = 1; i < N; ++i) {
    tb_trace_i("i: %d", i);

    // tb_int_t error_array[11] = {77118,
    //                             77119,
    //                             230604,
    //                             230605,
    //                             307134,
    //                             307135,
    //                             307137,
    //                             307138,
    //                             307139,
    //                             307143,
    //                             307144};
    // // skip all i in error_array
    // tb_int_t error = 0;
    // for (tb_int_t j = 0; j < 11; ++j) {
    //   if (i == error_array[j]) {
    //     error = 1;
    //     break;
    //   }
    // }
    // if (error) {
    //   continue;
    // }

    tb_int_t oct_idx = node_offsets[i];
    const tb_int_t n_new_nodes = rt_node_counts[i];

    // for each nodes, make n_new_nodes
    for (tb_int_t j = 0; j < n_new_nodes - 1; ++j) {
      const tb_int_t level = rt_prefixN[i] / 3 - j;

      const morton_t node_prefix = codes[i] >> (MORTON_BITS - (3 * level));
      const tb_int_t child_idx = node_prefix & 0b111;
      const tb_int_t parent = oct_idx + 1;

      set_child(&oct_nodes[parent], oct_idx, child_idx);

      vec4 ret;
      morton32_to_xyz(
          &ret, node_prefix << (MORTON_BITS - (3 * level)), min_coord, range);
      glm_vec4_copy(ret, oct_nodes[oct_idx].corner);

      // each cell is half the size of the level above it
      oct_nodes[oct_idx].cell_size =
          range / (tb_float_t)(1 << (level - root_level));
      oct_idx = parent;
    }

    if (n_new_nodes > 0) {
      tb_int_t rt_parent = rt_parents[i];

      if (i == 77118) {
        // print all information of this node
        printf("i: %d\n", i);
        printf("rt_prefixN[i]: %d\n", rt_prefixN[i]);
        printf("rt_node_counts[i]: %d\n", rt_node_counts[i]);
        printf("rt_parents[i]: %d\n", rt_parents[i]);
        printf("node_offsets[i]: %d\n", node_offsets[i]);
        printf("codes[i]: %d\n", codes[i]);
        printf("oct_idx: %d\n", oct_idx);
        // print rt_node_counts[rt_parent]
        printf("rt_node_counts[rt_parent]: %d\n", rt_node_counts[rt_parent]);
      }

      while (rt_node_counts[rt_parent] == 0) {
        rt_parent = rt_parents[rt_parent];
      }

      if (i == 77118) {
        printf("-------------------: %d\n", rt_parent);
      }

      const tb_int_t oct_parent = node_offsets[rt_parent];
      const tb_int_t top_level = rt_prefixN[i] / 3 - n_new_nodes + 1;
      const morton_t top_node_prefix =
          codes[i] >> (MORTON_BITS - (3 * top_level));

      const tb_int_t child_idx = top_node_prefix & 0b111;

      set_child(&oct_nodes[oct_parent], oct_idx, child_idx);

      morton32_to_xyz(&oct_nodes[oct_idx].corner,
                      top_node_prefix << (MORTON_BITS - (3 * top_level)),
                      min_coord,
                      range);
      oct_nodes[oct_idx].cell_size =
          range / (tb_float_t)(1 << (top_level - root_level));
    }
  }
}