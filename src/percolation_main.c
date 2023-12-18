#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"

/* get random number in range [0, 1] */
static float rand_float(){
  return (float)(rand() & CSTD_RAND_MAX) / (float)CSTD_RAND_MAX;
}

/* find root and path compression */
static int findroot(Graph* g, int i){
  if (g->ptr[i]<0) return i;  // return index of root node
  return g->ptr[i] = findroot(g, g->ptr[i]);  // recursively go to root node, plus: do path-compression on the fly
}

/* merge two graph fragments in g->ptr representation, return new root node index */
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
  memset(g->erasure, 0, g->nnode * sizeof(bool));
  memset(g->syndrome, 0, g->nnode * sizeof(bool));
  memset(g->error, 0, g->nnode * sizeof(bool));
  int num_syndromes = 0;
  for(int i=0; i < g->nnode; i++){
    if(g->is_qbt[i]){
      g->erasure[i] = (rand_float() < p_erasure);
      if ((g->erasure[i] && rand_float() < 0.5) || (rand_float() < p_error)) { // erasure makes error with 50%
        g->error[i] = 1;
        num_syndromes += (g->syndrome[g->nn[i*g->num_nb_max]] ? -1 : +1);
        g->syndrome[g->nn[i*g->num_nb_max]] = (g->syndrome[g->nn[i*g->num_nb_max]] + 1) % 2; // 1st syndrome neighbor
        num_syndromes += (g->syndrome[g->nn[i*g->num_nb_max+1]] ? -1 : +1);
        g->syndrome[g->nn[i*g->num_nb_max+1]] = (g->syndrome[g->nn[i*g->num_nb_max+1]] + 1) % 2; // 2nd syndrome neighbor
      }
    }
  }
  memcpy(g->parity, g->syndrome, g->nnode * sizeof(bool)); // syndrome and parity of cluster starts as the same thing (when all nodes are isolated)
  return num_syndromes;
}

/* get clusters with even number of syndromes by breadth-first traversal
   num_syndromes: number of syndromes signalling errors */
int get_even_clusters_bfs(Graph* g, int num_syndromes){
  int* bfs_list = malloc(g->nnode * sizeof(int));
  int* bfs_list_skipped = malloc(g->nnode * sizeof(int));
  int bfs_next = 0; // next free position in bfs_list
  int bfs_next_skipped = 0; // next free position in bfs_list_skipped
  memset(g->visited, 0, g->nnode * sizeof(bool));
  for(int i=0; i < g->nnode; i++){
    if (g->erasure[i]){ // erasures 1st (if only erasure errors one is done after BFS over this)
      bfs_list[bfs_next++] = i; // seed
      g->visited[i] = true; // mark as visited to avoid adding to bfs_list twice
    }
    g->ptr[i] = -1; // all isolated nodes in beginning
  }
  int num_erasure = bfs_next;
  for(int i=0; i < g->nnode; i++){
    if (g->syndrome[i]){ // syndromes 2nd
      bfs_list[bfs_next++] = i; // seed
      g->visited[i] = true; // mark as visited to avoid adding to bfs_list twice
    }
  }
  g->num_parity = num_syndromes; // number of unpaired syndromes
  int bfs_pos = 0; // current position for BFS
  while(g->num_parity > 0){
    if(bfs_pos == bfs_next){ // rare case that a skipped element is still needed
      free(bfs_list);
      bfs_list = bfs_list_skipped;
      bfs_list_skipped = malloc(g->nnode * sizeof(int));
      bfs_next = bfs_next_skipped;
      bfs_next_skipped = 0;
      bfs_pos = 0;
      num_erasure = 0; // all erasures visited at this point
    }
    int n = bfs_list[bfs_pos];
    int r_n = findroot(g, n);
    if(g->parity[r_n] || bfs_pos<num_erasure){ // only grow odd cluster || grow for initial erasure qubits
      for(int i=0; i<g->len_nb[n]; i++){
        int nb = g->nn[n*g->num_nb_max + i];
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

/* get forest spaning the erasure/syndrome clusters */
Forest get_forest(Graph* g){
  Forest f = new_forest(g->nnode);
  memset(f.visited, 0, f.nnode * sizeof(bool)); // can only change from 0->1
  memset(f.leaf, 1, f.nnode * sizeof(bool)); // can only change from 1->0
  for(int i=0; i<f.nnode; i++) f.parent[i] = -1; // -1 indicated tree root
  int* l_bfs = malloc(g->nnode * sizeof(int)); // another time BFS to build forest
  for(int n=0; n<g->nnode; n++){
    if(g->visited[n] && !f.visited[n]){ // n is root point of tree
      int r_n = findroot(g, n);
      int bfs_next = 1; // next free position
      int bfs_pos = 0; // current position for BFS
      l_bfs[0] = n;
      f.visited[n] = 1;
      while(bfs_pos < bfs_next){ // crawl cluster with BFS to create spanning tree
        int m = l_bfs[bfs_pos];
        for(uint8_t j=0; j<g->len_nb[m]; j++){ // look for neighbors that are not in tree yet
          int nb = g->nn[m*g->num_nb_max+j];
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
  memset(g->decode, 0, g->nnode * sizeof(bool));
  int num_leaf = 0;
  int* l_leaf = malloc(f->nnode*sizeof(int));
  for(int i=0; i<f->nnode; i++) if(f->leaf[i] && f->visited[i]) l_leaf[num_leaf++] = i;
  while(num_leaf > 0){
    int l = l_leaf[num_leaf-1];
    f->visited[l] = 0; // unvisit
    if(g->is_qbt[l]){ // leaf is qubit/edge
      int v = f->parent[l];
      int u = (g->nn[l*g->num_nb_max+0] != v ? g->nn[l*g->num_nb_max+0] : g->nn[l*g->num_nb_max+1]); // pendent vertex
      int r = findroot(g, l);
      if(g->syndrome[u] && findroot(g, u) == r){ // flip pendent vertex if it belongs to component
        if(v >= 0) g->syndrome[v] = !g->syndrome[v]; // flip
        g->syndrome[u] = !g->syndrome[u]; // flip
        g->decode[l] = 1; // decoder output
        if(print) printf("%i \n", l); // print decoder output
      }
      if(v < 0){ // v < 0 can happen is l is erased and root
        num_leaf--;
        continue;
      }
      bool is_leaf = true;
      for(uint8_t j=0; j<g->len_nb[v]; j++){
        int nb = g->nn[v*g->num_nb_max+j];
        if(f->visited[nb] && f->parent[nb] == v){ // there is child node that is still visited --> v is no leaf
          is_leaf=false;
          break;
        }
      }
      if(is_leaf) l_leaf[num_leaf-1] = v;
      else num_leaf--;
    } else { // leaf is syndrome/vertex
      int p = f->parent[l];
      if(p >= 0){ // there is a parent (l is no root)
        if(f->parent[p] >= 0) l_leaf[num_leaf-1] = p; // parent is not a root
        else if(!f->visited[g->nn[p*g->num_nb_max+0]] && !f->visited[g->nn[p*g->num_nb_max+1]]) l_leaf[num_leaf-1] = p; // parent is root but last node in cluster
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
  memset(g->syndrome, 0, g->nnode * sizeof(bool));
  int num_syndromes = 0;
  for(int i=0; i < g->nnode; i++){
    if(g->is_qbt[i]){
      if (g->error[i] ^ g->decode[i]) {
        num_syndromes += (g->syndrome[g->nn[i*g->num_nb_max]] ? -1 : +1);
        g->syndrome[g->nn[i*g->num_nb_max]] = (g->syndrome[g->nn[i*g->num_nb_max]] + 1) % 2; // 1st syndrome neighbor
        num_syndromes += (g->syndrome[g->nn[i*g->num_nb_max+1]] ? -1 : +1);
        g->syndrome[g->nn[i*g->num_nb_max+1]] = (g->syndrome[g->nn[i*g->num_nb_max+1]] + 1) % 2; // 2nd syndrome neighbor
      }
    }
  }
  return num_syndromes;
}

void collect_graph_and_decode(int nnode, int num_syndrome, uint8_t num_nb_max, int* nn, uint8_t* len_nb, bool* is_qbt, bool* syndrome, bool* erasure, bool* decode){
  Graph g;
  g.ptr = malloc(nnode * sizeof(int)); // several meanings: (if ptr[i]>0: parent index ("pointer"), elif ptr[i]<0: syndrome parity of component, qubits and syndromes
  g.nn = nn; // neighbors of a node (TBD: has a lot of zeros for tanner graph due to different vertex degrees)
  g.len_nb = len_nb; // until which index there are neighbors (255 neighbors max)
  g.is_qbt = is_qbt; // 0: syndrome, 1: qubit
  g.num_nb_max = num_nb_max; // maximum number of neighbors per node
  g.nnode = nnode; // number of nodes (qubits + syndromes)
  g.bfs_list = malloc(nnode * sizeof(int));
  g.visited = malloc(nnode * sizeof(bool)); // node visited (e.g. in BFS)
  g.syndrome = syndrome; // syndrome (for node type 0)
  g.erasure = erasure; // erasure (for node type 1)
  g.error = NULL; // error (for node type 1)
  g.parity = malloc(nnode * sizeof(bool)); // parity of syndromes in cluster (has meaning only for root node), 0: even number of syndromes
  g.decode = decode; // decoder output
  g.crr_surf_x = NULL; // coorelation surafce 1 (for checking logical error)
  g.crr_surf_y = NULL; // coorelation surafce 2 (for checking logical error)
  g.num_parity = 0; // number of unpaired syndromes

  get_even_clusters_bfs(&g, num_syndrome);
  Forest f = get_forest(&g);
  peel_forest(&f, &g, false);
  free_forest(&f);

  free(g.ptr);
  free(g.bfs_list);
  free(g.visited);
  free(g.parity);
}
