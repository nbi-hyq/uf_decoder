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
  if (-g->ptr[r1] > g->big) g->big = -g->ptr[r1]; // update largest component size
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
  int bfs_next = 0; // next free position in g->bfs_list
  memset(g->visited, 0, g->nnode * sizeof(bool));
  for(int i=0; i < g->nnode; i++){
    if (g->erasure[i] || g->syndrome[i]) g->bfs_list[bfs_next++] = i; // seed
    g->ptr[i] = -1; // all isolated nodes in beginning
  }
  g->big = 0; // size of largest connected component
  g->num_parity = num_syndromes; // number of unpaired syndromes
  int bfs_pos = 0; // current position for BFS
  while(g->num_parity > 0){
    int n = g->bfs_list[bfs_pos];
    g->visited[n] = true;
    int r_n = findroot(g, n);
    for(int i=0; i<g->len_nb[n]; i++){
      int nb = g->nn[n*g->num_nb_max + i];
      int r_nb = findroot(g, nb);
      if(r_n != r_nb) r_n = merge_root(g, r_n, r_nb);
      if (g->visited[nb] == false) {
        g->bfs_list[bfs_next++] = nb;
        g->visited[nb] = true;
      }
    }
    bfs_pos++;
  }
  return bfs_next;
}

/* get forest spaning the erasure/syndrome clusters
   num_bfs: final number in bfs list in get_even_clusters_bfs (bfs_next) */
Forest get_forest(Graph* g, int num_bfs){
  Forest f = new_forest(g->nnode);
  memset(f.visited, 0, f.nnode * sizeof(bool)); // can only change from 0->1
  memset(f.leaf, 1, f.nnode * sizeof(bool)); // can only change from 1->0
  for(int i=0; i<f.nnode; i++) f.parent[i] = -1; // -1 indicated tree root
  int* l_bfs = malloc(num_bfs * sizeof(int)); // another time BFS to build forest
  for(int i=0; i<num_bfs; i++){
    int n = g->bfs_list[i];
    if(!f.visited[n]){ // n is root point of tree
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
            //printf("%i --> %i\n",nb, m); printf("%i %i\n",bfs_next, num_bfs);
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
    if(g->is_qbt[l]){ // leaf is qubit/edge       TBD: simplify with       int p = f->parent[l];
      int u = g->nn[l*g->num_nb_max+0];
      int v = g->nn[l*g->num_nb_max+1];
      if (f->visited[u]){ // if still visited, u cannot be the pendent vertex
        u = v;
        v = g->nn[l*g->num_nb_max+0];
      }
      if(g->syndrome[u]){
        g->syndrome[v] = !g->syndrome[v]; // flip
        g->decode[l] = 1; // decoder output
        if(print) printf("%i \n", l); // print decoder output
      }
      uint8_t cnt=0;
      for(uint8_t j=0; j<g->len_nb[v]; j++) if(f->visited[g->nn[v*g->num_nb_max+j]]) cnt++; // count visited neighbors to decide if leaf or not
      if(cnt > 1) num_leaf--; // one leaf less
      else l_leaf[num_leaf-1] = v; // vertex becomes leaf
    } else { // leaf is syndrome/vertex
      /*int p = f->parent[l];
      if(p >= 0) l_leaf[num_leaf-1] = p;
      else num_leaf--;*/
      uint8_t cnt=0;
      for(uint8_t j=0; j<g->len_nb[l]; j++){
        int nb = g->nn[l*g->num_nb_max+j];
        if(f->visited[nb]){
          l_leaf[num_leaf-1] = nb; // qubit on edge becomes new leaf
          cnt++;
          break;
        }
      }
      if(cnt == 0) num_leaf--;
    }
  }
  free(l_leaf);
  return 0;
}

