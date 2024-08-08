#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "../inc/graph_type.h"

/* one-sided Tanner graph with syndrome and qubit nodes
   (in data-structures such as ptr that contain data qubity and syndromes, the data qubits must come first)*/
Graph new_graph(int n_qbt, uint8_t num_nb_max_qbt, int n_syndr, uint8_t num_nb_max_syndr){
  Graph g; // data qubits come before syndromes in indexing
  int nnode = n_qbt + n_syndr;
  g.ptr = malloc(nnode * sizeof(int)); // if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: size of cluster, qubits and syndromes
  g.visited = malloc(nnode * sizeof(bool)); // node visited (added to bfs_list)
  g.parity = malloc(nnode * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.num_qbt = malloc(nnode * sizeof(int)); // number of data qubits in cluster (only for ldpc decoder)
  g.len_nb = malloc(nnode); // until which index there are neighbors (255 neighbors max)
  g.nn_qbt = malloc(n_qbt * (size_t)num_nb_max_qbt * sizeof(int)); // neighbors of a data qubit
  g.nn_syndr = malloc(n_syndr * (size_t)num_nb_max_syndr * sizeof(int)); // neighbors of a syndrome
  g.syndrome = malloc(n_syndr * sizeof(bool)); // syndrome (for node type 0)
  g.erasure = malloc(n_qbt * sizeof(bool)); // erasure (for node type 1)
  g.error = malloc(n_qbt * sizeof(bool)); // error (for node type 1)
  g.decode = malloc(n_qbt * sizeof(bool)); // decoder output
  g.crr_surf_x = NULL; // coorelation surafce 1 (for checking logical error)
  g.crr_surf_y = NULL; // coorelation surafce 2 (for checking logical error)
  g.num_parity = 0; // number of unpaired syndromes
  g.num_nb_max_qbt = num_nb_max_qbt; // maximum number of neighbors per node (of data qubits)
  g.num_nb_max_syndr = num_nb_max_syndr; // maximum number of neighbors per node (of syndromes)
  g.n_qbt = n_qbt; // number of data qubits
  g.n_syndr = n_syndr; // number of syndromes
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
  free(g->nn_qbt);
  free(g->nn_syndr);
  free(g->len_nb);
  free(g->syndrome);
  free(g->erasure);
  free(g->error);
  free(g->parity);
  free(g->decode);
  free(g->visited);
  free(g->num_qbt);
  free(g->crr_surf_x);
  free(g->crr_surf_y);
}

void free_forest(Forest* f){
  free(f->parent);
  free(f->visited);
  free(f->leaf);
}
