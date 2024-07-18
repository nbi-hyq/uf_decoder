#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/decoder_main.h"

/* simulate square lattice toric code with Pauli errors and erasure (using Algorithm 2) */
int main(){
  srand(37577);
  int num_rep = 10000; // number of repetitions for averaging
  int num_erasure = 20; // number of different erasure values
  int num_err = 120; // number of different error values

  for(int lsize=10; lsize<41; lsize+=10){
    Graph g = get_2d_toric_code(lsize);
    for(int e=0; e<num_erasure; e++){
      float p_erasure = 0.0 + 0.5 * e / num_erasure; // probability of erasure
      for(int i=0; i<num_err; i++){
        float p_err = 0.0 + 0.15 * i / num_err; // probability of error
        int cnt_error = 0;
        for(int rep=0; rep<num_rep; rep++){
          int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
          get_even_clusters_bfs_skip_store_root(&g, num_syndromes); // Algorithm 2
          Forest f = get_forest(&g);
          peel_forest(&f, &g, false);
          bool logical_error_x = 0;
          bool logical_error_y = 0;
          for(int c=0; c<g.num_crr_x; c++){
            logical_error_x ^= g.decode[g.crr_surf_x[c]];
            logical_error_x ^= g.error[g.crr_surf_x[c]];
          }
          for(int c=0; c<g.num_crr_y; c++){
            logical_error_y ^= g.decode[g.crr_surf_y[c]];
            logical_error_y ^= g.error[g.crr_surf_y[c]];
          }
          if(logical_error_x || logical_error_y){
            cnt_error++;
            //visualize_decode(&g, lsize);
          }
          free_forest(&f);
        }
        printf("%f %f %f\n", p_erasure, p_err, (double)cnt_error / (double)num_rep);
      }
    }
    free_graph(&g);
  }
  return 0;
}
