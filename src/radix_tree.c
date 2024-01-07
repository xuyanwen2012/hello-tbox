#include "radix_tree.h"

#include <stdio.h>

#include "assert.h"
#include "morton.h"
#include "tbox/prefix/type.h"

#define TEST_LOG2_CEIL(x)               \
  do {                                  \
    tb_int_t a = log2_ceil_u32(x);      \
    tb_float_t b = tb_ceil(tb_log2(x)); \
    assert(a == b);                     \
  } while (0)

void init_brt_nodes(brt_nodes_t* nodes,
                    const morton_t* morton_codes,
                    int n_nodes) {
  nodes->morton_codes = morton_codes;

  // allocate memory for brt nodes
  nodes->hasLeafLeft = tb_nalloc_type(n_nodes, tb_bool_t);
  nodes->hasLeafRight = tb_nalloc_type(n_nodes, tb_bool_t);
  nodes->prefixN = tb_nalloc_type(n_nodes, tb_uint8_t);
  nodes->leftChild = tb_nalloc_type(n_nodes, int);
  nodes->parent = tb_nalloc_type(n_nodes, int);

  // print total memory used
  tb_size_t total_memory = 0;
  total_memory += n_nodes * sizeof(tb_bool_t);
  total_memory += n_nodes * sizeof(tb_bool_t);
  total_memory += n_nodes * sizeof(tb_uint8_t);
  total_memory += n_nodes * sizeof(int);
  total_memory += n_nodes * sizeof(int);
  printf("Total memory used for brt_nodes: %f MB\n",
         (float)total_memory / 1024 / 1024);
}

void free_brt_nodes(brt_nodes_t* nodes) {
  // dont free morton_codes
  tb_free(nodes->hasLeafLeft);
  tb_free(nodes->hasLeafRight);
  tb_free(nodes->prefixN);
  tb_free(nodes->leftChild);
  tb_free(nodes->parent);
}

void init_radix_tree(radix_tree_t* tree,
                     const morton_t* sorted_mortons_keys,
                     int n_unique_keys,
                     float min_coord,
                     float max_coord) {
  tree->n_pts = n_unique_keys;
  tree->n_nodes = n_unique_keys - 1;
  tree->min_coord = min_coord;
  tree->max_coord = max_coord;

  init_brt_nodes(&tree->d_tree, sorted_mortons_keys, tree->n_nodes);
}

void free_radix_tree(radix_tree_t* tree) { free_brt_nodes(&tree->d_tree); }

/* /////////////////////////////////////////////////////////////////////////////
 * radix tree construction
 */

uint32_t ceil_div_u32(uint32_t a, uint32_t b) {
  tb_assert(b != 0);
  return (a + b - 1) / b;
}

int_fast8_t delta_u32(tb_uint32_t a, tb_uint32_t b) {
  const tb_uint32_t bit1_mask = ((tb_uint32_t)1) << (sizeof(a) * 8 - 1);
  tb_assert((a & bit1_mask) == 0);
  tb_assert((b & bit1_mask) == 0);
  return __builtin_clz(a ^ b) - 1;
}

tb_int_t log2_ceil_u32(tb_uint32_t x) {
  tb_assert_static(sizeof(x) ==
                   sizeof(unsigned int));  // __clz(x) is for long long int"

  // Counting from LSB to MSB, number of bits before last '1'
  // This is floor(log(x))
  const int n_lower_bits = (8 * sizeof(x)) - __builtin_clz(x) - 1;

  // Add 1 if 2^n_lower_bits is less than x
  //     (i.e. we rounded down because x was not a power of 2)
  return n_lower_bits + ((1 << n_lower_bits) < x);
}

void build_radix_tree(radix_tree_t* tree) {
  // alias
  const int n = tree->n_pts;
  const morton_t* codes = tree->d_tree.morton_codes;
  tb_bool_t* has_leaf_left = tree->d_tree.hasLeafLeft;
  tb_bool_t* has_leaf_right = tree->d_tree.hasLeafRight;
  tb_uint8_t* prefix_n = tree->d_tree.prefixN;
  int* left_child = tree->d_tree.leftChild;
  int* parent = tree->d_tree.parent;

#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; i++) {
    const morton_t code_i = codes[i];
    // Determine direction of the range (+1 or -1)
    int d;
    if (i == 0) {
      d = 1;
    } else {
      const int_fast8_t delta_diff_right = delta_u32(code_i, codes[i + 1]);
      const int_fast8_t delta_diff_left = delta_u32(code_i, codes[i - 1]);
      const int direction_difference = delta_diff_right - delta_diff_left;
      d = (direction_difference > 0) - (direction_difference < 0);
    }

    // Compute upper bound for the length of the range

    morton_t l = 0;
    if (i == 0) {
      // First node is root, covering whole tree
      l = n - 1;
    } else {
      const int_fast8_t delta_min = delta_u32(code_i, codes[i - d]);
      morton_t l_max = 2;
      // Cast to ptrdiff_t so in case the result is negative (since d is +/- 1),
      // we can catch it and not index out of bounds
      while (i + (tb_ptrdiff_t)l_max * d >= 0 && i + l_max * d <= n &&
             delta_u32(code_i, codes[i + l_max * d]) > delta_min) {
        l_max *= 2;
      }
      const int l_cutoff = (d == -1) ? i : n - i;
      int t;
      int divisor;
      // Find the other end using binary search
      for (t = l_max / 2, divisor = 2; t >= 1;
           divisor *= 2, t = l_max / divisor) {
        if (l + t <= l_cutoff &&
            delta_u32(code_i, codes[i + (l + t) * d]) > delta_min) {
          l += t;
        }
      }
    }

    int j = i + l * d;

    // Find the split position using binary search
    const int_fast8_t delta_node = delta_u32(codes[i], codes[j]);
    prefix_n[i] = delta_node;
    int s = 0;
    const int max_divisor = 1 << log2_ceil_u32(l);
    int divisor = 2;
    const int s_cutoff = (d == -1) ? i - 1 : n - i - 1;
    for (int t = ceil_div_u32(l, 2); divisor <= max_divisor;
         divisor <<= 1, t = ceil_div_u32(l, divisor)) {
      if (s + t <= s_cutoff &&
          delta_u32(code_i, codes[i + (s + t) * d]) > delta_node) {
        s += t;
      }
    }

    // Split position
    const int gamma = i + s * d + tb_min(d, 0);
    left_child[i] = gamma;
    has_leaf_left[i] = (tb_min(i, j) == gamma);
    has_leaf_right[i] = (tb_max(i, j) == gamma + 1);
    // Set parents of left and right children, if they aren't leaves
    // can't set this node as parent of its leaves, because the
    // leaf also represents an internal node with a differnent parent
    if (!has_leaf_left[i]) {
      parent[gamma] = i;
    }
    if (!has_leaf_right[i]) {
      parent[gamma + 1] = i;
    }
  }

  tb_trace_i("build_radix_tree done");
}