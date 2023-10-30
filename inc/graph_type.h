#ifndef GRAPH_TYPE
#define GRAPH_TYPE

typedef struct {
  int* ptr;
  int* nn;
  uint8_t* len_nb;
  bool* is_qbt;
  bool* syndrome;
  bool* erasure;
  bool* error;
  bool* parity;
  int nnode, num_edges, num_parity, big;
  uint8_t num_nb_max;
} Graph;

Graph new_graph(int nnode, uint8_t num_nb_max);
void free_graph(Graph* g);

#endif
