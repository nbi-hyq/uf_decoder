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
  srand(time(NULL));
  float p_erasure = 0.00; // probability of erasure
  int num_rep = 20; // number of repetitions for averaging

  for(int lsize = 200; lsize<900; lsize+=200){
    Graph g = get_2d_toric_code(lsize);
    int num_err = 100;
    for(int i=0; i<num_err; i++){
      float p_err = 0.0 + 0.2 * i / num_err;
      float avg = 0;
      for(int rep=0; rep<num_rep; rep++){
        int num_syndromes = apply_erasure_and_error(&g, p_erasure, p_err);
        get_even_clusters_bfs(&g, num_syndromes);
        avg += (float)g.big / (float)g.nnode;
      }
      printf("%f %f\n", p_err, avg / num_rep);
    }
    free_graph(&g);
  }
  return 0;
}
