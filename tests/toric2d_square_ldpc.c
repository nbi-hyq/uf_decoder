#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/decoder_main.h"
#include "../inc/decoder_ldpc.h"

int main(){
  srand(1749347);
  int lsize = 10; // size in each dimension
  int n_rep = 4; // repeat several times
  int num_err = 15;
  int num_erasure = 6;

  Graph g = get_2d_toric_code(lsize);
  int r = validate_graph(&g);

  for(int rep=0; rep<n_rep; rep++){
  for(int e=0; e<num_erasure; e++){
    float p_erasure = 0.0 + 0.6 * e / num_erasure; // probability of erasure
    for(int i=0; i<num_err; i++){
      float p_err = 0.0 + 0.15 * i / num_err;
      int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
      ldpc_syndrome_validation_and_decode(&g, num_syndromes); // Algorithm 3
      visualize_decode(&g, lsize);
      int c = check_correction_general(&g);
      printf("%f %i\n", p_err, c);
      r |= c;
    }
  }
  }

  free_graph(&g);
  return r;
}

