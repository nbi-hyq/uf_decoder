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
  int num_seed = 0; // number of syndromes and erasures
  int* bfs_list = malloc(g->nnode * sizeof(int));
  bool* visited = malloc(g->nnode * sizeof(bool));
  memset(visited, 0, g->nnode * sizeof(bool));
  for(int i=0; i < g->nnode; i++){
    if (g->erasure[i] || g->syndrome[i]) bfs_list[num_seed++] = i;
    g->ptr[i] = -1; // all isolated nodes in beginning
  }
  g->big = 0; // size of largest connected component
  g->num_parity = num_syndromes; // number of unpaired syndromes

  int bfs_pos = 0; // current position for BFS
  int bfs_next = num_seed; // next free position in bfs_list
  while(g->num_parity > 0){ // TBD: check: || bfs_pos < num_seed
    int n = bfs_list[bfs_pos];
    visited[n] = true;
    int r_n = findroot(g, n);
    for(int i=0; i<g->len_nb[n]; i++){
      int nb = g->nn[n*g->num_nb_max + i];
      int r_nb = findroot(g, nb);
      if(r_n != r_nb) r_n = merge_root(g, r_n, r_nb);
      if (visited[nb] == false) {
        bfs_list[bfs_next++] = nb;
        visited[nb] = true;
      }
    }
    bfs_pos++;
  }
  free(bfs_list);
  free(visited);
  return bfs_next;
}

