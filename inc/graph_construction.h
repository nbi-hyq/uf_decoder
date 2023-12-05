#ifndef GRAPH_CONSTRUCT
#define GRAPH_CONSTRUCT

Graph get_2d_toric_code(int lsize);
Graph get_2d_triangular_toric_code(int lsize);

int validate_graph(Graph* g);
int visualize_error(Graph* g, int size);
int visualize_decode(Graph* g, int size);
int visualize_forest(Forest* f, int size);

#endif
