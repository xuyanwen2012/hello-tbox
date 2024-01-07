#include "brt.h"

#include "assert.h"
#include "morton.h"
#include "tbox/prefix/type.h"

#define TEST_LOG2_CEIL(x)               \
  do {                                  \
    tb_int_t a = log2_ceil(x);          \
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

uint32_t ceil_div_uint32(uint32_t a, uint32_t b) {
  assert(b != 0);
  return (a + b - 1) / b;
}

int_fast8_t delta_u32(tb_uint32_t a, tb_uint32_t b) {
  const tb_uint32_t bit1_mask = ((tb_uint32_t)1) << (sizeof(a) * 8 - 1);
  assert((a & bit1_mask) == 0);
  assert((b & bit1_mask) == 0);
  return __builtin_clz(a ^ b) - 1;
}

tb_int_t log2_ceil(tb_uint32_t x) {
  // Counting from LSB to MSB, number of bits before last '1'
  // This is floor(log(x))
  const int n_lower_bits = (8 * sizeof(x)) - __builtin_clz(x) - 1;

  // Add 1 if 2^n_lower_bits is less than x
  //     (i.e. we rounded down because x was not a power of 2)
  return n_lower_bits + ((1 << n_lower_bits) < x);
}

void build_brt_nodes(brt_nodes_t* nodes,
                     const morton_t* morton_codes,
                     int n_nodes) {}