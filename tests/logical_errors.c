#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/decoder_main.h"

/* test square lattice toric code threshold without erasures (for Algorithm 2) */
int main(){
  srand(67242);
  int num_rep = 150; // number of repetitions for averaging
  int num_err = 45; // number of different error values
  float p_erasure = 0.0; // probability of erasure
  double* logicalArry = malloc(num_err * sizeof(double)); 
  int r = 0;

  for(int lsize=10; lsize<21; lsize+=10){
    Graph g = get_2d_toric_code(lsize);
    for(int i=0; i<num_err; i++){
      float p_err = 0.0 + 0.15 * i / num_err; // probability of error
      int cnt_error = 0;
      for(int rep=0; rep<num_rep; rep++){
        int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
        /* apply decoding algorithm */
        get_even_clusters_bfs_skip_store_root(&g, num_syndromes);
        Forest f = get_forest(&g);
        peel_forest(&f, &g);
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
      double pLogical = (double)cnt_error / (double)num_rep;
      if(lsize == 10) logicalArry[i] = pLogical; // store to compare to value at larger size
      if(p_err < 0.05 && lsize == 20 && pLogical > logicalArry[i]) r |= 1;
      if(p_err < 0.035 && pLogical > p_err) r |= 1;
      if(p_err < 0.05 && lsize == 10 && pLogical > 0.1) r |= 1;
      if(p_err < 0.07 && lsize == 20 && pLogical > 0.1) r |= 1;
      if(p_err > 0.093 && lsize == 10 && pLogical < 0.1) r |= 1;
      if(p_err > 0.095 && lsize == 20 && pLogical < 0.1) r |= 1;
    }
    free_graph(&g);
  }
  free(logicalArry);
  return r;
}
