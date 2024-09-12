#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/decoder_main.h"
#include "inc/decoder_ldpc.h"

/* simulate square lattice toric code with Pauli errors and erasure, use general ldpc decoder (Algorithm 3) for this */
int main(){
  srand(358938);
  int num_rep = 2000; // number of repetitions for averaging
  int num_erasure = 10; // number of different erasure values
  int num_err = 50; // number of different error values

  for(int lsize=10; lsize<21; lsize+=10){
    Graph g = get_2d_toric_code(lsize);
    for(int e=0; e<num_erasure; e++){
      float p_erasure = 0.0 + 0.4 * e / num_erasure; // probability of erasure
      for(int i=0; i<num_err; i++){
        float p_err = 0.0 + 0.10 * i / num_err; // probability of error
        int cnt_error = 0;
        for(int rep=0; rep<num_rep; rep++){
          int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
          ldpc_syndrome_validation_and_decode(&g, num_syndromes); // Algorithm 3
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
        printf("%f %f %f\n", p_erasure, p_err, (double)cnt_error / (double)num_rep);
      }
    }
    free_graph(&g);
  }
  return 0;
}
