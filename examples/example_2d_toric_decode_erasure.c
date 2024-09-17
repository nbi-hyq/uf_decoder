#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/decoder_main.h"

/* simulate square lattice toric code with only erasures (using Algorithm 4) */
int main(){
  srand(37577);
  int num_rep = 10000; // number of repetitions for averaging
  int num_erasure = 100; // number of different erasure values
  float p_err = 0.0;

  for(int lsize=10; lsize<41; lsize+=10){
    Graph g = get_2d_toric_code(lsize);
    for(int e=0; e<num_erasure; e++){
      float p_erasure = 0.0 + 0.8 * e / num_erasure; // probability of erasure
      int cnt_error = 0;
      for(int rep=0; rep<num_rep; rep++){
        int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
        /* apply decoding algorithm */
        get_even_clusters_bfs_skip(&g, num_syndromes); // Algorithm 4
        Forest f = get_forest(&g);
        peel_forest(&f, &g, false);
        free_forest(&f);
        /* check for logical errors */
        for(int l=0; l<g.num_logicals; l++){
          bool logical_error = 0;
          for(int c=0; c<g.logical_weight[l]; c++){
            logical_error ^= g.decode[g.logicals[l][c]];
            logical_error ^= g.error[g.logicals[l][c]];
          }
          if(logical_error){
            cnt_error++;
            break;
          }
        }
      }
      printf("%f %f\n", p_erasure, (double)cnt_error / (double)num_rep);
    }
    free_graph(&g);
  }
  return 0;
}
