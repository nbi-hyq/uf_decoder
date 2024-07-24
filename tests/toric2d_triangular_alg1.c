#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/decoder_main.h"

int main(){
  srand(5637934);
  int lsize = 10; // size in each dimension
  int n_rep = 20; // repeat several times
  int num_err = 20;
  int num_erasure = 16;

  Graph g = get_2d_triangular_toric_code(lsize);
  int r = validate_graph(&g);

  for(int rep=0; rep<n_rep; rep++){
  for(int e=0; e<num_erasure; e++){
    float p_erasure = 0.0 + 0.8 * e / num_erasure; // probability of erasure
    for(int i=0; i<num_err; i++){
      float p_err = 0.0 + 0.2 * i / num_err;
      int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
      get_even_clusters_bfs(&g, num_syndromes); // Algorithm 1
      Forest f = get_forest(&g);
      peel_forest(&f, &g, false);
      visualize_decode(&g, lsize);
      int c = check_correction(&g);
      printf("%f %i\n", p_err, c);
      r |= c;
      free_forest(&f);
    }
  }
  }

  free_graph(&g);
  return r;
}

