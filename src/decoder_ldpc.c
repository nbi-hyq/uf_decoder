#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/stabilizer_main.h"

/* get the next prime number above the next power of two (for hash table) */
static int getNextPrimePower2(int N){
  if(N <= 16) return 17;
  else if(N <= 32) return 37;
  else if(N <= 64) return 67;
  else if(N <= 128) return 131;
  else if(N <= 256) return 257;
  else if(N <= 512) return 521;
  else if(N <= 1024) return 1031;
  else if (N <= 2048) return 2053;
  else if (N <= 4096) return 4099;
  else if (N <= 8192) return 8209;
  else if (N <= 16384) return 16411;
  else if (N <= 32768) return 32771;
  else if (N <= 65536) return 65537;
  else if (N <= 131072) return 131101;
  else if (N <= 262144) return 262147;
  else if (N <= 524288) return 524309;
  else if (N <= 1048576) return 1048583;
  else if (N <= 2097152) return 2097169;
  else if (N <= 4194304) return 4194319;
  else if (N <= 8388608) return 8388617;
  else if (N <= 16777216) return 16777259;
  else if (N <= 33554432) return 33554467;
  else if (N <= 67108864) return 67108879;
  else if (N <= 134217728) return 134217757;
  else if (N <= 268435456) return 268435459;
  else if (N <= 536870912) return 536870923;
  else if (N <= 1073741824) return 1073741827;
  else return 2147483647; // integer maximum
}

/* find position of number in open addressing hash table. Wite elements into hash table if it does not exist yet
idx: number that is searched, hashTab: hashtable with -1 indicating empty spot
return: 0 when idx was not in hash table, 1 if it was */
static bool isVisited(int idx, int* hashTab, int tabSize){
  int h = idx % tabSize;
  while(hashTab[h] >= 0 && hashTab[h] != idx) h = (h + 1) % tabSize;
  if(hashTab[h] == -1){
    hashTab[h] = idx;
    return 0;
  } else return 1;
}

/* given global index idx, lookup local index in reduced H matrix */
static int idxLookup(int idx, int* hashTab, int* idxMap, int tabSize){
  int h = idx % tabSize;
  while(hashTab[h] >= 0 && hashTab[h] != idx) h = (h + 1) % tabSize; // TBD rm 1st
  return idxMap[h];
}

/* write local index into datastruture */
static void idxWrite(int idx, int* hashTab, int* idxMap, int tabSize, int idxLocal){
  int h = idx % tabSize;
  while(hashTab[h] >= 0 && hashTab[h] != idx) h = (h + 1) % tabSize;
  idxMap[h] = idxLocal;
}

/* find root and path compression */
static int findroot(Graph* g, int i){
  if (g->ptr[i]<0) return i;  // return index of root node
  return g->ptr[i] = findroot(g, g->ptr[i]);  // recursively go to root node, plus: do path-compression on the fly
}

/* decode cluster given root node of it, return if decoding is possible */
bool ldpc_decode_cluster(Graph* g, int root){
  /* go through cluster breadth first and build reduced H-matrix */
  int n_row = (- g->ptr[root]) - g->num_qbt[root]; // number of checks in cluster
  int n_col = g->num_qbt[root] + 1; // number of columns (number of data qubits + 1)
  WMat w = new_wmat(n_row, n_col); // reduced H-matrix for cluster (indexing in the order of the breadth-first traversal of the cluster)
  memset(w.mat, 0, w.num_bytes);
  int* lBf = malloc((- g->ptr[root]) * sizeof(int));
  int tableSize = getNextPrimePower2(2*(- g->ptr[root])); // size of the hash table (2x number of elements is generous)
  int* visited = malloc(sizeof(int)*tableSize); // hash table to know visited nodes of the considered cluster (without using a global array)
  int* globalToMat = malloc(sizeof(int)*tableSize); // index mapping from global index to reduced matrix (column or row index depending on whether it is syndrome or data qubit)
  for(int i=0; i<tableSize; i++) visited[i] = -1; // mark all parts of hash table as empty by -1
  int cntqbtLen = 0; // count number of qubits below index len_lBf
  int cntqbtPos = 0; // count number of qubits below index bf_pos
  int bf_pos = 0;
  int len_lBf = 1;
  lBf[0] = root;
  isVisited(root, visited, tableSize);
  idxWrite(root, visited, globalToMat, tableSize, 0);
  if(g->is_qbt[root]){
    g->decode[root] = 0; // reset previous decoding
    cntqbtLen = 1;
  }
  while(bf_pos < len_lBf){
    int idxLocal = idxLookup(lBf[bf_pos], visited, globalToMat, tableSize);
    if(!g->is_qbt[lBf[bf_pos]] && g->syndrome[lBf[bf_pos]]) write_matrix_position_bit(w.mat, idxLocal, n_col - 1, w.num_blocks_row, 1); // last column specifying right part of H*e=s
    else if (g->is_qbt[lBf[bf_pos]]) g->decode[lBf[bf_pos]] = 0; // reset previous decoding
    for(uint8_t i=0; i<g->len_nb[lBf[bf_pos]]; i++){
      int nb = g->nn[lBf[bf_pos]*g->num_nb_max + i];
      if(findroot(g, nb) == root){
        int idxLocalNb;
        if(!isVisited(nb, visited, tableSize)){
          if(g->is_qbt[nb]) idxLocalNb = cntqbtLen++;
          else idxLocalNb = len_lBf - cntqbtLen;
          idxWrite(nb, visited, globalToMat, tableSize, idxLocalNb); // store local index for new neighbor
          lBf[len_lBf++] = nb;
        } else idxLocalNb = idxLookup(nb, visited, globalToMat, tableSize);
        if(g->is_qbt[lBf[bf_pos]]) write_matrix_position_bit(w.mat, idxLocalNb, idxLocal, w.num_blocks_row, 1); // write columnwise
      }
    }
    if(g->is_qbt[bf_pos++]) cntqbtPos++;
  }

  /* decode the cluster if possible by Gaussian elimination */
  GaussElimin_bit(w.mat, n_row, n_col, n_col);
  bool* decode = malloc(g->num_qbt[root]); // store decoding at local indices
  memset(decode, 0, g->num_qbt[root]);
  bool decodable = true;
  uint64_t num_blocks_row = n_col / (8 * sizeof(dtypeBlk)) + 1;
  for(int r = n_row-1; r>=0; r--){
    if(read_matrix_position_bit(w.mat, r, n_col - 1, num_blocks_row)){ // if there is a syndrome that has to be fixed
      int c = 0; //TBD: min(r, n_col - 1)
      while(!read_matrix_position_bit(w.mat, r, c, num_blocks_row) && c < n_col - 1) c++; // find (not yet used) qubit that belongs to tagged syndrome and can be used for its correction
      if(c == n_col - 1) {decodable = false; break;} // syndome cannot be corrected
      decode[c] = 1;
      for(int r2=r-1; r2>=0; r2--){ // flip other sydromes that are affected
        if(read_matrix_position_bit(w.mat, r2, c, num_blocks_row)) flip_matrix_position_bit(w.mat, r2, n_col - 1, num_blocks_row);
      }
    }
  }

  /* apply decoding at global indices (TBD: maybe better with mapping from local to global index) */
  for(int i=0; i<tableSize; i++){
    if(visited[i] >= 0 && g->is_qbt[visited[i]] && decode[globalToMat[i]]) g->decode[visited[i]] = 1;
  }

  free_wmat(w);
  free(lBf);
  free(globalToMat);
  free(visited);
  free(decode);
  return decodable;
}

/* merge two graph fragments in g->ptr representation, return new root node index, update validity of merged cluster
   (assumes for parity update that r1 is cluster from where stuff is grown (assumed to be invalid until update at end) */
static int merge_root(Graph* g, int r1, int r2){
  if (g->parity[r2]) g->num_parity -= 1;
  if (g->ptr[r1] > g->ptr[r2]){
    g->num_qbt[r2] += g->num_qbt[r1]; // add number of data qubits in smaller cluster to bigger cluster
    g->ptr[r2] += g->ptr[r1]; // add size of smaller component to bigger one
    g->ptr[r1] = r2; // attach component with root r1 to larger component with root r2
    r1 = r2; // update root
  } else {
    g->num_qbt[r1] += g->num_qbt[r2]; // add number of data qubits in smaller cluster to bigger cluster
    g->ptr[r1] += g->ptr[r2]; // add size of smaller component to bigger one
    g->ptr[r2] = r1; // attach component with root r2 to larger component with root r1
  }
  return r1;
}

/* after growing cluster, check again if it is valid */
static void update_cluster_validity(Graph* g, int root){
  g->parity[root] = !ldpc_decode_cluster(g, root);
  if(!g->parity[root]) g->num_parity -= 1;
}

/* used for storing skipped nodes in linked list connected to cluster root */
struct nodeSk {
  int i_node;
  struct nodeSk* next;
};
typedef struct nodeSk nodeSk;

static void connect_node(nodeSk* node, int i_next){
  nodeSk* next = malloc(sizeof(nodeSk));
  node->next = next;
  next->i_node = i_next;
  next->next = NULL;
}

static void free_nodeSk_list(nodeSk* node){
  while(node->next != NULL){
    nodeSk* h = node;
    node = node->next;
    free(h);
  }
  free(node);
}

/* Algorithm 3: syndrome validation using breadth-first Tanner graph traversal to grow clusters and Gaussian elimination for validating them */
void ldpc_syndrome_validation_and_decode(Graph* g, int num_syndromes){
  memset(g->decode, 0, g->nnode * sizeof(bool));
  nodeSk** a_skipped = malloc(g->nnode * sizeof(nodeSk*));
  for(int i=0; i<g->nnode; i++) a_skipped[i] = NULL;
  nodeSk** a_skipped_last = malloc(g->nnode * sizeof(nodeSk*));
  int* bf_list = malloc(g->nnode * sizeof(int));
  int bf_next = 0; // next free position in bf_list
  memset(g->visited, 0, g->nnode * sizeof(bool));
  for(int i=0; i < g->nnode; i++){
    if (g->erasure[i]){ // erasures 1st (if only erasure errors, one is done after doing one step from here in Tanner graph)
      bf_list[bf_next++] = i; // seed
      g->visited[i] = true; // mark as visited to avoid adding to bf_list twice
    }
    g->ptr[i] = -1; // all isolated nodes in beginning
    g->num_qbt[i] = g->is_qbt[i]; // number of data qubits in cluster (0 or 1)
  }
  int num_erasure = bf_next;
  for(int i=0; i < g->nnode; i++){
    if (g->syndrome[i]){ // syndromes 2nd
      bf_list[bf_next++] = i; // seed
      g->visited[i] = true; // mark as visited to avoid adding to bf_list twice
    }
  }
  g->num_parity = num_syndromes; // num_parity has here the meaning of number of invalid clusters
  int bf_pos = 0; // current position for breadth-first graph traversal

  /* grow erasures first by one step (helps as one is done when errors happen only on erased qubits) */
  while (bf_pos < num_erasure){
    int n = bf_list[bf_pos];
    int r_n = findroot(g, n);
    if (!g->parity[r_n]) g->num_parity += 1; // if start from valid cluster, compensate g->num_parity -= 1 in update_cluster_validity
    for(uint8_t i=0; i<g->len_nb[n]; i++){
      int nb = g->nn[n*g->num_nb_max + i];
      int r_nb = findroot(g, nb);
      if(r_n != r_nb) r_n = merge_root(g, r_n, r_nb);
      if (g->visited[nb] == false) {
        bf_list[bf_next++] = nb;
        g->visited[nb] = true;
      }
    }
    update_cluster_validity(g, r_n);
    bf_pos++;
  }

  /* grow clusters with method derived from breadth-first traversal (grow here in double steps such that the cluster boundary is always syndrome checks) */
  while(g->num_parity > 0){
    int n = bf_list[bf_pos];
    int r_n = findroot(g, n);
    if(g->parity[r_n]){ // only grow invalid cluster
      for(uint8_t i=0; i<g->len_nb[n]; i++){
        int nb = g->nn[n*g->num_nb_max + i]; // this neighbor is always a data qubit (no check qubit), so it cannot have been skipped nor be part of another cluster
        int r_nb = findroot(g, nb);
        if(r_n != r_nb){
          r_n = merge_root(g, r_n, r_nb);
          for(uint8_t j=0; j<g->len_nb[nb]; j++){
            int nb2 = g->nn[nb*g->num_nb_max + j]; // this neighbor is always a check qubit, so it can have been skipped or be part of another cluster
            int r_nb2 = findroot(g, nb2);
          
            /* if r_nb2 is skipped cluster, recover it */
            if (a_skipped[r_nb2] != NULL){
              nodeSk* node = a_skipped[r_nb2];
              bf_list[bf_next] = node->i_node;
              bf_next = (bf_next + 1) % g->nnode;
              while(node->next != NULL){
                node = node->next;
                bf_list[bf_next] = node->i_node;
                bf_next = (bf_next + 1) % g->nnode;
              }
              free_nodeSk_list(a_skipped[r_nb2]);
              a_skipped[r_nb2] = NULL;
            }

            if(r_n != r_nb2) r_n = merge_root(g, r_n, r_nb2);
            if (g->visited[nb2] == false) {
              bf_list[bf_next] = nb2;
              bf_next = (bf_next + 1) % g->nnode;
              g->visited[nb2] = true;
            }
          }
        }
      }
      update_cluster_validity(g, r_n);
    } else {
      if (a_skipped[r_n] == NULL){
        a_skipped[r_n] = malloc(sizeof(nodeSk));
        a_skipped[r_n]->i_node = n;
        a_skipped[r_n]->next = NULL;
        a_skipped_last[r_n] = a_skipped[r_n];
      } else {
        connect_node(a_skipped_last[r_n], n);
        a_skipped_last[r_n] = a_skipped_last[r_n]->next;
      }
    }
    bf_pos = (bf_pos + 1) % g->nnode;
  }
  for(int i=0; i<g->nnode; i++) if (a_skipped[i] != NULL) free_nodeSk_list(a_skipped[i]);
  free(a_skipped);
  free(a_skipped_last);
  free(bf_list);
}

/* check if the decoding has worked in the sense that no syndromes are left */
int check_correction_general(Graph* g){
  memset(g->syndrome, 0, g->nnode * sizeof(bool));
  int num_syndromes = 0;
  for(int i=0; i < g->nnode; i++){
    if(g->is_qbt[i] && (g->error[i] ^ g->decode[i])){
      for(uint8_t j=0; j<g->len_nb[i]; j++){
        num_syndromes += (g->syndrome[g->nn[i*g->num_nb_max+j]] ? -1 : +1);
        g->syndrome[g->nn[i*g->num_nb_max+j]] = (g->syndrome[g->nn[i*g->num_nb_max+j]] + 1) % 2;
      }
    }
  }
  return num_syndromes;
}

/* given graph and syndrome, compute decoding for general ldpc code */
void ldpc_collect_graph_and_decode(int nnode, uint8_t num_nb_max, int* nn, uint8_t* len_nb, bool* is_qbt, bool* syndrome, bool* erasure, bool* decode){
  Graph g;
  g.ptr = malloc(nnode * sizeof(int)); // several meanings: (if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: syndrome parity of component, qubits and syndromes
  g.num_qbt = malloc(nnode * sizeof(int)); // number of data qubits in cluster
  g.nn = nn; // neighbors of a node (TBD: has a lot of zeros for tanner graph due to different vertex degrees)
  g.len_nb = len_nb; // until which index there are neighbors (255 neighbors max)
  g.is_qbt = is_qbt; // 0: syndrome, 1: qubit
  g.num_nb_max = num_nb_max; // maximum number of neighbors per node
  g.nnode = nnode; // number of nodes (qubits + syndromes)
  g.visited = malloc(nnode * sizeof(bool)); // node visited (e.g. in breadth-first graph traversal)
  g.syndrome = syndrome; // syndrome (for node type 0)
  g.erasure = erasure; // erasure (for node type 1)
  g.error = NULL; // error (for node type 1)
  g.parity = malloc(nnode * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.decode = decode; // decoder output
  g.crr_surf_x = NULL; // correlation surface 1 (for checking logical error)
  g.crr_surf_y = NULL; // correlation surface 2 (for checking logical error)
  memcpy(g.parity, g.syndrome, g.nnode * sizeof(bool)); // syndrome and parity of cluster starts as the same thing (when all nodes are isolated)

  int num_syndrome = 0;
  for(int i=0; i<g.nnode; i++) if(syndrome[i]) num_syndrome++; // no check of !is_qbt done here
  ldpc_syndrome_validation_and_decode(&g, num_syndrome);

  free(g.ptr);
  free(g.num_qbt);
  free(g.visited);
  free(g.parity);
}

/* given graph and syndrome, compute decoding in batches of nrep repetitions (for general ldpc code) */
void ldpc_collect_graph_and_decode_batch(int nnode, uint8_t num_nb_max, int* nn, uint8_t* len_nb, bool* is_qbt, bool* syndrome, bool* erasure, bool* decode, int nrep){
  Graph g;
  g.ptr = malloc(nnode * sizeof(int)); // several meanings: (if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: syndrome parity of component, qubits and syndromes
  g.num_qbt = malloc(nnode * sizeof(int)); // number of data qubits in cluster
  g.nn = nn; // neighbors of a node (TBD: has a lot of zeros for tanner graph due to different vertex degrees)
  g.len_nb = len_nb; // until which index there are neighbors (255 neighbors max)
  g.is_qbt = is_qbt; // 0: syndrome, 1: qubit
  g.num_nb_max = num_nb_max; // maximum number of neighbors per node
  g.nnode = nnode; // number of nodes (qubits + syndromes)
  g.visited = malloc(nnode * sizeof(bool)); // node visited (e.g. in breadth-first graph traversal)
  g.error = NULL; // error (for node type 1)
  g.parity = malloc(nnode * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.crr_surf_x = NULL; // correlation surface 1 (for checking logical error)
  g.crr_surf_y = NULL; // correlation surface 2 (for checking logical error)

  for(int r=0; r<nrep; r++){
    g.syndrome = syndrome + r*g.nnode; // syndrome (for node type 0)
    g.decode = decode + r*g.nnode; // decoder output
    g.erasure = erasure + r*g.nnode;
    memcpy(g.parity, g.syndrome, g.nnode * sizeof(bool)); // syndrome and parity of cluster starts as the same thing (when all nodes are isolated)
    int num_syndrome = 0;
    for(int i=0; i<g.nnode; i++) if(g.syndrome[i]) num_syndrome++; // no check of !is_qbt done here
    ldpc_syndrome_validation_and_decode(&g, num_syndrome);
  }

  free(g.ptr);
  free(g.num_qbt);
  free(g.visited);
  free(g.parity);
}

