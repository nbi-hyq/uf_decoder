#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "../inc/graph_type.h"

/* one-sided Tanner graph with syndrome and qubit nodes */
Graph new_graph(int nnode, uint8_t num_nb_max){
  Graph g;
  g.ptr = malloc(nnode * sizeof(int)); // several meanings: (if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: syndrome parity of component, qubits and syndromes
  g.nn = malloc(nnode * (size_t)num_nb_max * sizeof(int)); // neighbors of a node (TBD: has a lot of zeros for tanner graph due to different vertex degrees)
  g.len_nb = malloc(nnode); // until which index there are neighbors (255 neighbors max)
  g.is_qbt = malloc(nnode * sizeof(bool)); // 0: syndrome, 1: qubit
  g.num_nb_max = num_nb_max; // maximum number of neighbors per node
  g.nnode = nnode; // number of nodes (qubits + syndromes)
  g.visited = malloc(nnode * sizeof(bool)); // node visited (added to bfs_list)
  g.syndrome = malloc(nnode * sizeof(bool)); // syndrome (for node type 0)
  g.erasure = malloc(nnode * sizeof(bool)); // erasure (for node type 1)
  g.error = malloc(nnode * sizeof(bool)); // error (for node type 1)
  g.parity = malloc(nnode * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.decode = malloc(nnode * sizeof(bool)); // decoder output
  g.crr_surf_x = NULL; // coorelation surafce 1 (for checking logical error)
  g.crr_surf_y = NULL; // coorelation surafce 2 (for checking logical error)
  g.num_parity = 0; // number of unpaired syndromes
  return g;
}

/* spanning forest for all errors */
Forest new_forest(int nnode){
  Forest f;
  f.nnode = nnode;
  f.parent = malloc(nnode * sizeof(int));
  f.visited = malloc(nnode * sizeof(bool));
  f.leaf = malloc(nnode * sizeof(bool));
  return f;
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
  free(g->decode);
  free(g->visited);
  free(g->crr_surf_x);
  free(g->crr_surf_y);
}

void free_forest(Forest* f){
  free(f->parent);
  free(f->visited);
  free(f->leaf);
}
