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
  int64_t lsize = 100; // size in each dimension
  Graph g = get_2d_toric_code(lsize);

  free_graph(&g);
  return 0;
}
