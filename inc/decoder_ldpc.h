#ifndef DECODER_LDPC
#define DECODER_LDPC

bool ldpc_decode_cluster(Graph* g, int root);
void ldpc_syndrome_validation_and_decode(Graph* g, int num_syndromes);
int check_correction_general(Graph* g);
void ldpc_collect_graph_and_decode(int n_qbt, int n_syndr, uint8_t num_nb_max_qbt, uint8_t num_nb_max_syndr, int* nn_qbt, int* nn_syndr, uint8_t* len_nb, bool* syndrome, bool* erasure, bool* decode);
void ldpc_collect_graph_and_decode_batch(int n_qbt, int n_syndr, uint8_t num_nb_max_qbt, uint8_t num_nb_max_syndr, int* nn_qbt, int* nn_syndr, uint8_t* len_nb, bool* syndrome, bool* erasure, bool* decode, int nrep);

#endif
