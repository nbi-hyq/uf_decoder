#ifndef GRAPH_TYPE
#define GRAPH_TYPE

typedef struct {
  int* ptr;
  int* nn;
  int* bfs_list;
  uint8_t* len_nb;
  bool* is_qbt;
  bool* syndrome;
  bool* erasure;
  bool* error;
  bool* parity;
  bool* visited;
  int nnode, num_edges, num_parity, big;
  uint8_t num_nb_max;
} Graph;

typedef struct {
  int* root;
  bool* visited;
  bool* leaf;
  int nnode;
} Forest;

Graph new_graph(int nnode, uint8_t num_nb_max);
void free_graph(Graph* g);
Forest new_forest(int nnode);
void free_forest(Forest* f);

#endif
