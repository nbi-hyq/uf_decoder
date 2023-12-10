#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(494839);
  int lsize = 10; // size in each dimension
  int num_err = 20;
  int num_erasure = 8;

  Graph g = get_2d_toric_code(lsize);
  int r = validate_graph(&g);

  for(int e=0; e<num_erasure; e++){
    float p_erasure = 0.0 + 0.8 * e / num_erasure; // probability of erasure
    for(int i=0; i<num_err; i++){
      float p_err = 0.0 + 0.2 * i / num_err;
      int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
      int num_bfs = get_even_clusters_bfs(&g, num_syndromes);
      Forest f = get_forest(&g, num_bfs);
      peel_forest(&f, &g, false);
      visualize_decode(&g, lsize);
      int c = check_correction(&g);
      printf("%f %i\n", p_err, c);
      r |= c;
      free_forest(&f);
    }
  }

  free_graph(&g);
  return r;
}

