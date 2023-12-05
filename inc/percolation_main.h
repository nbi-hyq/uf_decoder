#ifndef PERCOLATION
#define PERCOLATION

int apply_erasure_and_error(Graph* g, float p_erasure, float p_error);
int get_even_clusters_bfs(Graph* g, int num_syndromes);
Forest get_forest(Graph* g, int num_bfs);
int peel_forest(Forest* f, Graph* g, bool print);

#endif
