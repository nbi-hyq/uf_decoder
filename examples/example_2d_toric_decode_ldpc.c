#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/decoder_main.h"
#include "inc/decoder_ldpc.h"

/* simulate square lattice toric code (decoding via general ldcp decoder) */
int main(){
  srand(83959790);
  float p_erasure = 0.1; // probability of erasure
  int lsize = 10;
  int num_err = 20;
  for(int i=0; i<num_err; i++){
    float p_err = 0.0 + 0.15 * i / num_err;
    printf("---------- p = %f\n", p_err);
    Graph g = get_2d_toric_code(lsize);
    int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
    ldpc_syndrome_validation_and_decode(&g, num_syndromes);
    visualize_error(&g, lsize);
    visualize_decode(&g, lsize);
    free_graph(&g);
  }
  return 0;
}
