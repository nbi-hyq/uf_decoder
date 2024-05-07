#ifndef DECODER_LDPC
#define DECODER_LDPC

bool ldpc_decode_cluster(Graph* g, int root);
void ldpc_syndrome_validation_and_decode(Graph* g, int num_syndromes);
int check_correction_general(Graph* g);
void ldpc_collect_graph_and_decode(int nnode, uint8_t num_nb_max, int* nn, uint8_t* len_nb, bool* is_qbt, bool* syndrome, bool* erasure, bool* decode);
void ldpc_collect_graph_and_decode_batch(int nnode, uint8_t num_nb_max, int* nn, uint8_t* len_nb, bool* is_qbt, bool* syndrome, bool* erasure, bool* decode, int nrep);

#endif
