#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/decoder_main.h"

/* simulate square lattice toric code, check how speed scales with size (for Algorithm 2) */
int main(){
  srand(23578838);
  float p_erasure = 0.0; // probability of erasure
  int nrep = 10000;
  int nerr = 9;

  for(int i=0; i<nerr; i++){
    float p_err = 0.01 + 0.01 * i;
    for(int lsize=2; lsize<20; lsize+=2){
      Graph g = get_2d_toric_code(lsize);
      struct timeval t0, t1;
      gettimeofday(&t0, NULL);

      for(int rep=0; rep<nrep; rep++){
        int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
        get_even_clusters_bfs_skip_store_root(&g, num_syndromes); // Algorithm 2
        Forest f = get_forest(&g);
        peel_forest(&f, &g);
        free_forest(&f);
      }

      gettimeofday(&t1, NULL);
      printf("%i %i %li\n", g.n_qbt, g.n_syndr, 1000000*(t1.tv_sec - t0.tv_sec) + t1.tv_usec - t0.tv_usec);
      free_graph(&g);
    }
  }

  return 0;
}
