#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/decoder_main.h"

int main(){
  srand(494839);
  int lsize = 100; // size in each dimension
  float p_erasure = 0.00; // probability of erasure

  Graph g = get_2d_triangular_toric_code(lsize);
  int r = validate_graph(&g);

  int num_err = 20;
  for(int i=0; i<num_err; i++){
    float p_err = 0.0 + 0.2 * i / num_err;
    int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
    get_even_clusters_bfs_skip(&g, num_syndromes);
    Forest f = get_forest(&g);
    peel_forest(&f, &g, false);
    r |= check_correction(&g);
    free_forest(&f);
  }

  free_graph(&g);
  return r;
}
