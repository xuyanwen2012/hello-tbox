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

#if defined(__GNUC__) || defined(__clang__)
#define CLZ(x) __builtin_clz(x)
#elif defined(_MSC_VER)
#define CLZ(x) _lzcnt_u32(x)
#else
#error "CLZ not supported on this platform"
#endif

void init_brt_nodes(brt_nodes_t* nodes,
                    const morton_t* morton_codes,
                    const tb_int_t n_nodes) {
  nodes->morton_codes = morton_codes;

  // allocate memory for brt nodes
  nodes->hasLeafLeft = tb_nalloc_type(n_nodes, tb_bool_t);
  nodes->hasLeafRight = tb_nalloc_type(n_nodes, tb_bool_t);
  nodes->prefixN = tb_nalloc_type(n_nodes, tb_uint8_t);
  nodes->leftChild = tb_nalloc_type(n_nodes, tb_int_t);
  nodes->parent = tb_nalloc_type(n_nodes, tb_int_t);

  // print total memory used
  tb_size_t total_memory = 0;
  total_memory += n_nodes * sizeof(tb_bool_t);
  total_memory += n_nodes * sizeof(tb_bool_t);
  total_memory += n_nodes * sizeof(tb_uint8_t);
  total_memory += n_nodes * sizeof(tb_int_t);
  total_memory += n_nodes * sizeof(tb_int_t);
  printf("Total memory used for brt_nodes: %f MB\n",
         (tb_float_t)total_memory / 1024 / 1024);
}

void free_brt_nodes(const brt_nodes_t* nodes) {
  // dont free morton_codes
  tb_free(nodes->hasLeafLeft);
  tb_free(nodes->hasLeafRight);
  tb_free(nodes->prefixN);
  tb_free(nodes->leftChild);
  tb_free(nodes->parent);
}

void init_radix_tree(radix_tree_t* tree,
                     const morton_t* sorted_morton_keys,
                     const tb_int_t n_unique_keys,
                     const tb_float_t min_coord,
                     const tb_float_t max_coord) {
  tree->n_pts = n_unique_keys;
  tree->n_nodes = n_unique_keys - 1;
  tree->min_coord = min_coord;
  tree->max_coord = max_coord;

  init_brt_nodes(&tree->d_tree, sorted_morton_keys, tree->n_nodes);
}

void free_radix_tree(const radix_tree_t* tree) {
  free_brt_nodes(&tree->d_tree);
}

/* /////////////////////////////////////////////////////////////////////////////
 * radix tree construction
 */

tb_uint32_t ceil_div_u32(const tb_uint32_t a, const tb_uint32_t b) {
  tb_assert(b != 0);
  return (a + b - 1) / b;
}

tb_uint8_t delta_u32(const tb_uint32_t a, const tb_uint32_t b) {
  const tb_uint32_t bit1_mask = ((tb_uint32_t)1) << (sizeof(a) * 8 - 1);
  tb_assert((a & bit1_mask) == 0);
  tb_assert((b & bit1_mask) == 0);
  return CLZ(a ^ b) - 1;
}

tb_int_t log2_ceil_u32(const tb_uint32_t x) {
  // Counting from LSB to MSB, number of bits before last '1'
  // This is floor(log(x))
  const tb_int_t n_lower_bits = (8 * sizeof(x)) - CLZ(x) - 1;

  // Add 1 if 2^n_lower_bits is less than x
  //     (i.e. we rounded down because x was not a power of 2)
  return n_lower_bits + ((1 << n_lower_bits) < x);
}

void build_radix_tree(radix_tree_t* tree) {
  // alias
  const tb_int_t n = tree->n_pts;
  const morton_t* codes = tree->d_tree.morton_codes;
  tb_bool_t* has_leaf_left = tree->d_tree.hasLeafLeft;
  tb_bool_t* has_leaf_right = tree->d_tree.hasLeafRight;
  tb_uint8_t* prefix_n = tree->d_tree.prefixN;
  tb_int_t* left_child = tree->d_tree.leftChild;
  tb_int_t* parent = tree->d_tree.parent;

#pragma omp parallel for schedule(static)
  for (tb_int_t i = 0; i < n; i++) {
    const morton_t code_i = codes[i];
    // Determine direction of the range (+1 or -1)
    tb_int_t d;
    if (i == 0) {
      d = 1;
    } else {
      const tb_uint8_t delta_diff_right = delta_u32(code_i, codes[i + 1]);
      const tb_uint8_t delta_diff_left = delta_u32(code_i, codes[i - 1]);
      const tb_int_t direction_difference = delta_diff_right - delta_diff_left;
      d = (direction_difference > 0) - (direction_difference < 0);
    }

    // Compute upper bound for the length of the range

    morton_t l = 0;
    if (i == 0) {
      // First node is root, covering whole tree
      l = n - 1;
    } else {
      const tb_uint8_t delta_min = delta_u32(code_i, codes[i - d]);
      morton_t l_max = 2;
      // Cast to ptrdiff_t so in case the result is negative (since d is +/- 1),
      // we can catch it and not index out of bounds
      while (i + (tb_ptrdiff_t)l_max * d >= 0 && i + l_max * d <= n &&
             delta_u32(code_i, codes[i + l_max * d]) > delta_min) {
        l_max *= 2;
      }
      const tb_int_t l_cutoff = (d == -1) ? i : n - i;
      tb_int_t t;
      tb_int_t divisor;
      // Find the other end using binary search
      for (t = l_max / 2, divisor = 2; t >= 1;
           divisor *= 2, t = l_max / divisor) {
        if (l + t <= l_cutoff &&
            delta_u32(code_i, codes[i + (l + t) * d]) > delta_min) {
          l += t;
        }
      }
    }

    const tb_int_t j = i + l * d;

    // Find the split position using binary search
    const tb_uint8_t delta_node = delta_u32(codes[i], codes[j]);
    prefix_n[i] = delta_node;
    tb_int_t s = 0;
    const tb_int_t max_divisor = 1 << log2_ceil_u32(l);
    tb_int_t divisor = 2;
    const tb_int_t s_cutoff = (d == -1) ? i - 1 : n - i - 1;
    for (tb_int_t t = ceil_div_u32(l, 2); divisor <= max_divisor;
         divisor <<= 1, t = ceil_div_u32(l, divisor)) {
      if (s + t <= s_cutoff &&
          delta_u32(code_i, codes[i + (s + t) * d]) > delta_node) {
        s += t;
      }
    }

    // Split position
    const tb_int_t gamma = i + s * d + tb_min(d, 0);
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