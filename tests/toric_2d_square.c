#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(494839);
  int lsize = 1000; // size in each dimension
  float p_erasure = 0.02; // probability of erasure
  float p_error = 0.05; // probability of error;

  Graph g = get_2d_toric_code(lsize);
  int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_error);
  int num_bfs = get_even_clusters_bfs(&g, num_syndromes);
  printf("%i %i\n", num_syndromes, num_bfs);

  int r = validate_graph(&g);
  free_graph(&g);
  return r;
}

