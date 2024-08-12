#ifndef DECODER_MAIN
#define DECODER_MAIN

int apply_erasure_and_error(Graph* g, float p_erasure, float p_error);
int get_even_clusters_bfs(Graph* g, int num_syndromes);
int get_even_clusters_bfs_skip(Graph* g, int num_syndromes);
int get_even_clusters_bfs_skip_store_root(Graph* g, int num_syndromes);
Forest get_forest(Graph* g);
int peel_forest(Forest* f, Graph* g, bool print);
int check_correction(Graph* g);
void collect_graph_and_decode(int n_qbt, int n_syndr, uint8_t num_nb_max_qbt, uint8_t num_nb_max_syndr, int* nn_qbt, int* nn_syndr, uint8_t* len_nb, bool* syndrome, bool* erasure, bool* decode);
void collect_graph_and_decode_batch(int n_qbt, int n_syndr, uint8_t num_nb_max_qbt, uint8_t num_nb_max_syndr, int* nn_qbt, int* nn_syndr, uint8_t* len_nb, bool* syndrome, bool* erasure, bool* decode, int nrep);

#endif
