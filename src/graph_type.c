#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "../inc/graph_type.h"

// one-sided Tanner graph with syndrome and qubit nodes
Graph new_graph(int nnode, uint8_t num_nb_max){
  Graph g;
  g.ptr = malloc(nnode * sizeof(int));  // several meanings: (if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: syndrome parity of component, qubits and syndromes
  g.nn = malloc(nnode * (size_t)num_nb_max * sizeof(int)); // neighbors of a node (has a lot of zeros for tannder graph -- improve)
  g.len_nb = malloc(nnode); // until which index there are neighbors (255 neighbors max)
  g.is_qbt = malloc(nnode * sizeof(bool)); // 0: syndrome, 1: qubit
  g.num_nb_max = num_nb_max; // maximum number of neighbors per node
  g.nnode = nnode; // number of nodes or central qubits

  g.syndrome = malloc(nnode * sizeof(bool)); // syndrome (for node type 1)
  g.erasure = malloc(nnode * sizeof(bool)); // erasure (for node type 0)
  g.error = malloc(nnode * sizeof(bool)); // erasure (for node type 0)
  g.parity = malloc(nnode * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.num_parity = 0; // number of unpaired syndromes
  g.big = 0; // size of largest connected cluster
  return g;
}

void free_graph(Graph* g){
  free(g->ptr);
  free(g->nn);
  free(g->len_nb);
  free(g->is_qbt);
  free(g->syndrome);
  free(g->erasure);
  free(g->error);
  free(g->parity);
}
