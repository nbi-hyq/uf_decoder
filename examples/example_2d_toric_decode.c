#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/decoder_main.h"

/* simulate square lattice toric code (using Algorithm 2) */
int main(){
  srand(83959790);
  float p_erasure = 0.0; // probability of erasure
  int lsize = 10;

  int num_err = 20;
  for(int i=0; i<num_err; i++){
    float p_err = 0.0 + 0.15 * i / num_err;
    printf("---------- p = %f\n", p_err);
    Graph g = get_2d_toric_code(lsize);
    int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
    get_even_clusters_bfs_skip_store_root(&g, num_syndromes); // Algorithm 2
    Forest f = get_forest(&g);
    visualize_error(&g, lsize);
    visualize_forest(&f, lsize);
    peel_forest(&f, &g, true);
    visualize_decode(&g, lsize);
    free_graph(&g);
    free_forest(&f);
  }

  return 0;
}
