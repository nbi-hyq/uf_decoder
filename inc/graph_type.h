#ifndef GRAPH_TYPE
#define GRAPH_TYPE

typedef struct {
  int* ptr;
  int* nn_qbt;
  int* nn_syndr;
  uint8_t* len_nb;
  bool* syndrome;
  bool* erasure;
  bool* error;
  bool* parity;
  bool* visited;
  bool* decode;
  int* num_qbt;
  int** logicals;
  int* logical_weight;
  int n_qbt, n_syndr, num_edges, num_parity, num_logicals;
  uint8_t num_nb_max_qbt, num_nb_max_syndr;
} Graph;

typedef struct {
  int* parent;
  bool* visited;
  bool* leaf;
  int nnode;
} Forest;

Graph new_graph(int n_qbt, uint8_t num_nb_max_qbt, int n_syndr, uint8_t num_nb_max_syndr);
void free_graph(Graph* g);
Forest new_forest(int nnode);
void free_forest(Forest* f);

#endif
