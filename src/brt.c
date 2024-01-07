#include "brt.h"

#include "tbox/prefix/type.h"

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
