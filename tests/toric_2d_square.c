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
  int64_t lsize = 1000; // size in each dimension
  int r = 0;

  Graph g = get_2d_toric_code(lsize);

  free_graph(&g);

  return r;
}

