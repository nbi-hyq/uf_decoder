#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"

/* get random number in range [0, 1] (TBD: limited resolution due to CSTD_RAND_MAX) */
static float rand_float(){
  return (float)(rand() & CSTD_RAND_MAX) / (float)CSTD_RAND_MAX;
}

/* find root and path compression */
static int findroot(Graph* g, int i){
  if (g->ptr[i]<0) return i;  // return index of root node
  return g->ptr[i] = findroot(g, g->ptr[i]);  // recursively go to root node, plus: do path-compression on the fly
}

/* merge two graph fragments in g->ptr representation, return new root node index, update cluster parity */
static int merge_root(Graph* g, int r1, int r2){
  if (g->parity[r1] && g->parity[r2]) g->num_parity -= 2;
  if (g->ptr[r1] > g->ptr[r2]){
    g->parity[r2] = (g->parity[r1] + g->parity[r2]) % 2;
    g->ptr[r2] += g->ptr[r1]; // add size of smaller component to bigger one
    g->ptr[r1] = r2; // attach component with root r1 to larger component with root r2
    r1 = r2; // update root
  } else {
    g->parity[r1] = (g->parity[r1] + g->parity[r2]) % 2;
    g->ptr[r1] += g->ptr[r2]; // add size of smaller component to bigger one
    g->ptr[r2] = r1; // attach component with root r2 to larger component with root r1
  }
  return r1;
}

/* apply erasures and Pauli errors, compute syndromes at same time
   p_erasure: probability of erasure
   p_error: probability of error */
int apply_erasure_and_error(Graph* g, float p_erasure, float p_error){
  memset(g->erasure, 0, g->n_qbt * sizeof(bool));
  memset(g->syndrome, 0, g->n_syndr * sizeof(bool));
  memset(g->error, 0, g->n_qbt * sizeof(bool));
  int num_syndromes = 0;
  for(int i=0; i < g->n_qbt; i++){
    g->erasure[i] = (rand_float() < p_erasure);
    if ((g->erasure[i] && rand_float() < 0.5) || (!g->erasure[i] && rand_float() < p_error)) { // erasure makes error with 50%
      g->error[i] = 1;
      for(uint8_t j=0; j<g->len_nb[i]; j++){
        if(g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt+j] - g->n_qbt]){
          num_syndromes--;
          g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt+j] - g->n_qbt] = 0;
        } else {
          num_syndromes++;
          g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt+j] - g->n_qbt] = 1;
        }
      }
    }
  }
  memset(g->parity, 0, g->n_qbt * sizeof(bool));
  memcpy(g->parity + g->n_qbt, g->syndrome, g->n_syndr * sizeof(bool)); // syndrome and parity of cluster starts as the same thing (when all nodes are isolated)
  return num_syndromes;
}

/* algorithm 1: get clusters with even number of syndromes by breadth-first traversal (oversimplified: also even clusters grow)
   num_syndromes: number of syndromes signalling errors */
int get_even_clusters_bfs(Graph* g, int num_syndromes){
  int* bfs_list = malloc((g->n_qbt + g->n_syndr) * sizeof(int));
  int bfs_next = 0; // next free position in bfs_list
  memset(g->visited, 0, (g->n_qbt + g->n_syndr) * sizeof(bool));
  for(int i=0; i < g->n_qbt; i++){
    if (g->erasure[i]){ // erasures 1st (if only erasure errors, one is done after doing one step from here in Tanner graph)
      bfs_list[bfs_next++] = i; // seed
      g->visited[i] = true; // mark as visited to avoid adding to bfs_list twice
    }
    g->ptr[i] = -1; // all isolated nodes in beginning
  }
  for(int i=0; i < g->n_syndr; i++){
    if (g->syndrome[i]){ // syndromes 2nd
      bfs_list[bfs_next++] = i + g->n_qbt; // seed
      g->visited[i + g->n_qbt] = true; // mark as visited to avoid adding to bfs_list twice
    }
    g->ptr[i + g->n_qbt] = -1; // all isolated nodes in beginning
  }
  g->num_parity = num_syndromes; // number of unpaired syndromes
  int bfs_pos = 0; // current position for breadth-first graph traversal
  while(g->num_parity > 0){
    int n = bfs_list[bfs_pos];
    int r_n = findroot(g, n);
    uint8_t num_nb_max;
    int* nn;
    int idx_arry; // index in syndrome or data qubit array
    if(n < g->n_qbt){
      num_nb_max = g->num_nb_max_qbt;
      nn = g->nn_qbt;
      idx_arry = n;
    } else {
      num_nb_max = g->num_nb_max_syndr;
      nn = g->nn_syndr;
      idx_arry = n - g->n_qbt;
    }
    for(uint8_t i=0; i<g->len_nb[n]; i++){
      int nb = nn[idx_arry*num_nb_max + i];
      int r_nb = findroot(g, nb);
      if(r_n != r_nb) r_n = merge_root(g, r_n, r_nb);
      if (g->visited[nb] == false) {
        bfs_list[bfs_next++] = nb;
        g->visited[nb] = true;
      }
    }
    bfs_pos++;
  }
  free(bfs_list);
  return bfs_next;
}

/* algorithm 4: get clusters with even number of syndromes by breadth-first traversal (skip even clusters)
   num_syndromes: number of syndromes signalling errors */
int get_even_clusters_bfs_skip(Graph* g, int num_syndromes){
  int* bfs_list = malloc((g->n_qbt + g->n_syndr) * sizeof(int));
  int* bfs_list_skipped = malloc((g->n_qbt + g->n_syndr) * sizeof(int));
  int bfs_next = 0; // next free position in bfs_list
  int bfs_next_skipped = 0; // next free position in bfs_list_skipped
  memset(g->visited, 0, (g->n_qbt + g->n_syndr) * sizeof(bool));
  for(int i=0; i < g->n_qbt; i++){
    if (g->erasure[i]){ // erasures 1st (if only erasure errors, one is done after doing one step from here in Tanner graph)
      bfs_list[bfs_next++] = i; // seed
      g->visited[i] = true; // mark as visited to avoid adding to bfs_list twice
    }
    g->ptr[i] = -1; // all isolated nodes in beginning
  }
  int num_erasure = bfs_next;
  for(int i=0; i < g->n_syndr; i++){
    if (g->syndrome[i]){ // syndromes 2nd
      bfs_list[bfs_next++] = i + g->n_qbt; // seed
      g->visited[i + g->n_qbt] = true; // mark as visited to avoid adding to bfs_list twice
    }
    g->ptr[i + g->n_qbt] = -1; // all isolated nodes in beginning
  }
  g->num_parity = num_syndromes; // number of unpaired syndromes
  int bfs_pos = 0; // current position for breadth-first graph traversal
  while(g->num_parity > 0){
    if(bfs_pos == bfs_next){ // rare case that a skipped element is still needed
      free(bfs_list);
      bfs_list = bfs_list_skipped;
      bfs_list_skipped = malloc((g->n_qbt + g->n_syndr) * sizeof(int));
      bfs_next = bfs_next_skipped;
      bfs_next_skipped = 0;
      bfs_pos = 0;
      num_erasure = 0; // all erasures visited at this point
    }
    int n = bfs_list[bfs_pos];
    int r_n = findroot(g, n);
    if(g->parity[r_n] || bfs_pos<num_erasure){ // only grow odd cluster || grow for initial erasure qubits
      uint8_t num_nb_max;
      int* nn;
      int idx_arry; // index in syndrome or data qubit array
      if(n < g->n_qbt){
        num_nb_max = g->num_nb_max_qbt;
        nn = g->nn_qbt;
        idx_arry = n;
      } else {
        num_nb_max = g->num_nb_max_syndr;
        nn = g->nn_syndr;
        idx_arry = n - g->n_qbt;
      }
      for(uint8_t i=0; i<g->len_nb[n]; i++){
        int nb = nn[idx_arry*num_nb_max + i];
        int r_nb = findroot(g, nb);
        if(r_n != r_nb) r_n = merge_root(g, r_n, r_nb); // TBD?: break when parity[r_n] changes?, but that would be half-skipped node
        if (g->visited[nb] == false) {
          bfs_list[bfs_next++] = nb;
          g->visited[nb] = true;
        }
      }
    } else {
      bfs_list_skipped[bfs_next_skipped++] = n; // store the skipped node
    }
    bfs_pos++;
  }
  free(bfs_list);
  free(bfs_list_skipped);
  return bfs_next;
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

/* algorithm 2: get clusters with even number of syndromes by breadth-first traversal (skip even clusters, but store skipped nodes with root for later)
   num_syndromes: number of syndromes signalling errors */
int get_even_clusters_bfs_skip_store_root(Graph* g, int num_syndromes){
  int nnode = g->n_qbt + g->n_syndr;
  nodeSk** a_skipped = malloc(nnode * sizeof(nodeSk*));
  for(int i=0; i<nnode; i++) a_skipped[i] = NULL;
  nodeSk** a_skipped_last = malloc(nnode * sizeof(nodeSk*));
  int* bfs_list = malloc(nnode * sizeof(int));
  int bfs_next = 0; // next free position in bfs_list
  memset(g->visited, 0, nnode * sizeof(bool));
  for(int i=0; i < g->n_qbt; i++){
    if (g->erasure[i]){ // erasures 1st (if only erasure errors, one is done after doing one step from here in Tanner graph)
      bfs_list[bfs_next++] = i; // seed
      g->visited[i] = true; // mark as visited to avoid adding to bfs_list twice
    }
    g->ptr[i] = -1; // all isolated nodes in beginning
  }
  int num_erasure = bfs_next;
  for(int i=0; i < g->n_syndr; i++){
    if (g->syndrome[i]){ // syndromes 2nd
      bfs_list[bfs_next++] = i + g->n_qbt; // seed
      g->visited[i + g->n_qbt] = true; // mark as visited to avoid adding to bfs_list twice
    }
    g->ptr[i + g->n_qbt] = -1; // all isolated nodes in beginning
  }
  g->num_parity = num_syndromes; // number of unpaired syndromes
  int bfs_pos = 0; // current position for breadth-first graph traversal

  /* grow erasures first by one step (helps as one is done when errors happen only on erased qubits) */
  while (bfs_pos < num_erasure){
    int n = bfs_list[bfs_pos];
    int r_n = findroot(g, n);
    uint8_t num_nb_max;
    int* nn;
    int idx_arry; // index in syndrome or data qubit array
    if(n < g->n_qbt){
      num_nb_max = g->num_nb_max_qbt;
      nn = g->nn_qbt;
      idx_arry = n;
    } else {
      num_nb_max = g->num_nb_max_syndr;
      nn = g->nn_syndr;
      idx_arry = n - g->n_qbt;
    }
    for(uint8_t i=0; i<g->len_nb[n]; i++){
      int nb = nn[idx_arry*num_nb_max + i];
      int r_nb = findroot(g, nb);
      if(r_n != r_nb) r_n = merge_root(g, r_n, r_nb);
      if (g->visited[nb] == false) {
        bfs_list[bfs_next++] = nb;
        g->visited[nb] = true;
      }
    }
    bfs_pos++;
  }

  /* grow clusters with method derived from breadth-first traversal */
  while(g->num_parity > 0){
    int n = bfs_list[bfs_pos];
    int r_n = findroot(g, n);
    if(g->parity[r_n]){ // only grow odd cluster
      uint8_t num_nb_max;
      int* nn;
      int idx_arry; // index in syndrome or data qubit array
      if(n < g->n_qbt){
        num_nb_max = g->num_nb_max_qbt;
        nn = g->nn_qbt;
        idx_arry = n;
      } else {
        num_nb_max = g->num_nb_max_syndr;
        nn = g->nn_syndr;
        idx_arry = n - g->n_qbt;
      }
      for(uint8_t i=0; i<g->len_nb[n]; i++){
        int nb = nn[idx_arry*num_nb_max + i];
        int r_nb = findroot(g, nb);

        /* if r_nb is skipped cluster, recover it */
        if (a_skipped[r_nb] != NULL){
          nodeSk* node = a_skipped[r_nb];
          bfs_list[bfs_next] = node->i_node;
          bfs_next = (bfs_next + 1) % nnode;
          while(node->next != NULL){
            node = node->next;
            bfs_list[bfs_next] = node->i_node;
            bfs_next = (bfs_next + 1) % nnode;
          }
          free_nodeSk_list(a_skipped[r_nb]);
          a_skipped[r_nb] = NULL;
        }

        if(r_n != r_nb) r_n = merge_root(g, r_n, r_nb); // TBD?: break when parity[r_n] changes?, but that would be half-skipped node
        if (g->visited[nb] == false) {
          bfs_list[bfs_next] = nb;
          bfs_next = (bfs_next + 1) % nnode;
          g->visited[nb] = true;
        }
      }
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
    bfs_pos = (bfs_pos + 1) % nnode;
  }
  for(int i=0; i<nnode; i++) if (a_skipped[i] != NULL) free_nodeSk_list(a_skipped[i]);
  free(a_skipped);
  free(a_skipped_last);
  free(bfs_list);
  return bfs_next;
}

/* get forest spanning the erasure clusters */
Forest get_forest(Graph* g){
  Forest f = new_forest(g->n_qbt + g->n_syndr);
  memset(f.visited, 0, f.nnode * sizeof(bool)); // can only change from 0->1
  memset(f.leaf, 1, f.nnode * sizeof(bool)); // can only change from 1->0
  for(int i=0; i<f.nnode; i++) f.parent[i] = -1; // -1 indicated tree root
  int* l_bfs = malloc(f.nnode * sizeof(int)); // another time breadth-first traversal to build forest
  for(int n=0; n<f.nnode; n++){
    if(g->visited[n] && !f.visited[n]){ // n is root point of tree
      int r_n = findroot(g, n);
      int bfs_next = 1; // next free position
      int bfs_pos = 0; // current position for breadth-first graph traversal
      l_bfs[0] = n;
      f.visited[n] = 1;
      while(bfs_pos < bfs_next){ // crawl cluster breadth-first to create spanning tree
        int m = l_bfs[bfs_pos];
        uint8_t num_nb_max;
        int* nn;
        int idx_arry; // index in syndrome or data qubit array
        if(m < g->n_qbt){
          num_nb_max = g->num_nb_max_qbt;
          nn = g->nn_qbt;
          idx_arry = m;
        } else {
          num_nb_max = g->num_nb_max_syndr;
          nn = g->nn_syndr;
          idx_arry = m - g->n_qbt;
        }
        for(uint8_t j=0; j<g->len_nb[m]; j++){ // look for neighbors that are not in tree yet
          int nb = nn[idx_arry*num_nb_max+j];
          if(!f.visited[nb] && r_n==findroot(g, nb)){ // neighbor not visited, same component
            l_bfs[bfs_next++] = nb;
            f.visited[nb] = 1;
            f.leaf[m] = 0;
            f.parent[nb] = m;
          }
        }
        bfs_pos++;
      }
    }
  }
  free(l_bfs);
  return f;
}

/* use the spanning forest to create the decoder output */
int peel_forest(Forest* f, Graph* g, bool print){
  memset(g->decode, 0, g->n_qbt * sizeof(bool));
  int num_leaf = 0;
  int* l_leaf = malloc(f->nnode*sizeof(int));
  for(int i=0; i<f->nnode; i++) if(f->leaf[i] && f->visited[i]) l_leaf[num_leaf++] = i;
  while(num_leaf > 0){
    int l = l_leaf[num_leaf-1];
    f->visited[l] = 0; // unvisit
    if(l < g->n_qbt){ // leaf is data qubit
      int v = f->parent[l];
      int u = (g->nn_qbt[l*g->num_nb_max_qbt+0] != v ? g->nn_qbt[l*g->num_nb_max_qbt+0] : g->nn_qbt[l*g->num_nb_max_qbt+1]); // pendent vertex
      int r = findroot(g, l);
      if(g->syndrome[u - g->n_qbt] && findroot(g, u) == r){ // flip pendent vertex if it belongs to component
        if(v >= 0) g->syndrome[v - g->n_qbt] = !g->syndrome[v - g->n_qbt]; // flip
        g->syndrome[u - g->n_qbt] = !g->syndrome[u - g->n_qbt]; // flip
        g->decode[l] = 1; // decoder output
        if(print) printf("%i \n", l); // print decoder output
      }
      if(v < 0){ // v < 0 can happen if l is erased and root
        num_leaf--;
        continue;
      }
      bool is_leaf = true;
      for(uint8_t j=0; j<g->len_nb[v]; j++){
        int nb = g->nn_syndr[(v - g->n_qbt)*g->num_nb_max_syndr+j];
        if(f->visited[nb] && f->parent[nb] == v){ // there is child node that is still visited --> v is no leaf
          is_leaf=false;
          break;
        }
      }
      if(is_leaf) l_leaf[num_leaf-1] = v;
      else num_leaf--;
    } else { // leaf is syndrome
      int p = f->parent[l];
      if(p >= 0){ // there is a parent (l is no root)
        if(f->parent[p] >= 0) l_leaf[num_leaf-1] = p; // parent is not a root
        else if(!f->visited[g->nn_qbt[p*g->num_nb_max_qbt+0]] && !f->visited[g->nn_qbt[p*g->num_nb_max_qbt+1]]) l_leaf[num_leaf-1] = p; // parent is root but last node in cluster
        else num_leaf--; // don't make root node leaf before it is last node in cluster
      }
      else num_leaf--;
    }
  }
  free(l_leaf);
  return 0;
}

/* check if the decoding has worked in the sense that no syndromes are left */
int check_correction(Graph* g){
  memset(g->syndrome, 0, g->n_syndr * sizeof(bool));
  int num_syndromes = 0;
  for(int i=0; i < g->n_qbt; i++){
    if (g->error[i] ^ g->decode[i]) {
      num_syndromes += (g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt] - g->n_qbt] ? -1 : +1);
      g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt] - g->n_qbt] = (g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt] - g->n_qbt] + 1) % 2; // 1st syndrome neighbor
      num_syndromes += (g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt + 1] - g->n_qbt] ? -1 : +1);
      g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt + 1] - g->n_qbt] = (g->syndrome[g->nn_qbt[i*g->num_nb_max_qbt + 1] - g->n_qbt] + 1) % 2; // 2nd syndrome neighbor
    }
  }
  return num_syndromes;
}

/* given graph and syndrome, compute decoding */
void collect_graph_and_decode(int n_qbt, int n_syndr, uint8_t num_nb_max_qbt, uint8_t num_nb_max_syndr, int* nn_qbt, int* nn_syndr, uint8_t* len_nb, bool* syndrome, bool* erasure, bool* decode){
  Graph g;
  g.n_qbt = n_qbt;
  g.n_syndr = n_syndr;
  g.ptr = malloc((n_qbt + n_syndr) * sizeof(int)); // if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: size of cluster, qubits and syndromes
  g.nn_qbt = nn_qbt; // neighbors of data qubit
  g.nn_syndr = nn_syndr; // neighbors of syndrome
  g.len_nb = len_nb; // until which index there are neighbors (255 neighbors max)
  g.num_nb_max_qbt = num_nb_max_qbt; // maximum number of neighbors per data qubit
  g.num_nb_max_syndr = num_nb_max_syndr; // maximum number of neighbors per syndrome
  g.visited = malloc((n_qbt + n_syndr) * sizeof(bool)); // node visited (e.g. in breadth-first traversal)
  g.syndrome = syndrome;
  g.erasure = erasure;
  g.error = NULL;
  g.parity = malloc((n_qbt + n_syndr) * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.decode = decode; // decoder output
  g.crr_surf_x = NULL; // correlation surface 1 (for checking logical error)
  g.crr_surf_y = NULL; // correlation surface 2 (for checking logical error)
  memset(g.parity, 0, g.n_qbt * sizeof(bool));
  memcpy(g.parity + g.n_qbt, g.syndrome, g.n_syndr * sizeof(bool)); // syndrome and parity of cluster starts as the same thing (when all nodes are isolated)

  int num_syndrome = 0;
  for(int i=0; i<g.n_syndr; i++) if(syndrome[i]) num_syndrome++;
  get_even_clusters_bfs_skip_store_root(&g, num_syndrome);
  Forest f = get_forest(&g);
  peel_forest(&f, &g, false);
  free_forest(&f);

  free(g.ptr);
  free(g.visited);
  free(g.parity);
}

/* given graph and syndrome, compute decoding in batches of nrep repetitions */
void collect_graph_and_decode_batch(int n_qbt, int n_syndr, uint8_t num_nb_max_qbt, uint8_t num_nb_max_syndr, int* nn_qbt, int* nn_syndr, uint8_t* len_nb, bool* syndrome, bool* erasure, bool* decode, int nrep){
  Graph g;
  g.n_qbt = n_qbt;
  g.n_syndr = n_syndr;
  g.ptr = malloc((n_qbt + n_syndr) * sizeof(int)); // if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: size of cluster, qubits and syndromes
  g.nn_qbt = nn_qbt; // neighbors of data qubit
  g.nn_syndr = nn_syndr; // neighbors of syndrome
  g.len_nb = len_nb; // until which index there are neighbors (255 neighbors max)
  g.num_nb_max_qbt = num_nb_max_qbt; // maximum number of neighbors per data qubit
  g.num_nb_max_syndr = num_nb_max_syndr; // maximum number of neighbors per syndrome
  g.visited = malloc((n_qbt + n_syndr) * sizeof(bool)); // node visited (e.g. in breadth-first traversal)
  g.error = NULL;
  g.parity = malloc((n_qbt + n_syndr) * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.crr_surf_x = NULL; // correlation surface 1 (for checking logical error)
  g.crr_surf_y = NULL; // correlation surface 2 (for checking logical error)

  for(int r=0; r<nrep; r++){
    g.syndrome = syndrome + r*g.n_syndr;
    g.decode = decode + r*g.n_qbt; // decoder output
    g.erasure = erasure + r*g.n_qbt;
    memset(g.parity, 0, g.n_qbt * sizeof(bool));
    memcpy(g.parity + g.n_qbt, g.syndrome, g.n_syndr * sizeof(bool)); // syndrome and parity of cluster starts as the same thing (when all nodes are isolated)
    int num_syndrome = 0;
    for(int i=0; i<g.n_syndr; i++) if(g.syndrome[i]) num_syndrome++;
    get_even_clusters_bfs_skip_store_root(&g, num_syndrome);
    Forest f = get_forest(&g);
    peel_forest(&f, &g, false);
    free_forest(&f);
  }

  free(g.ptr);
  free(g.visited);
  free(g.parity);
}

