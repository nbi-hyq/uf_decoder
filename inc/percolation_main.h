#ifndef PERCOLATION
#define PERCOLATION

int apply_erasure_and_error(Graph* g, float p_erasure, float p_error);
int get_even_clusters_bfs(Graph* g, int num_syndromes);
Forest get_forest(Graph* g);
int peel_forest(Forest* f, Graph* g, bool print);
int check_correction(Graph* g);
void collect_graph_and_decode(int nnode, uint8_t num_nb_max, int* nn, uint8_t* len_nb, bool* is_qbt, bool* syndrome, bool* erasure, bool* decode);

#endif
