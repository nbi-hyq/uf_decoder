#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* simulate square lattice toric code */
int main(){
  srand(23578838);
  float p_erasure = 0.0; // probability of erasure
  int lsize = 50;

  int num_err = 5;
  for(int i=0; i<num_err; i++){
    float p_err = 0.01 + 0.05 * i / num_err;
    printf("---------- p = %f\n", p_err);
    Graph g = get_2d_toric_code(lsize);
    for(int rep=0; rep<100000; rep++){
      int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
      get_even_clusters_bfs(&g, num_syndromes);
      Forest f = get_forest(&g);
      peel_forest(&f, &g, false);
      free_forest(&f);
    }
    free_graph(&g);
  }

  return 0;
}
