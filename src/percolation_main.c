#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#define CSTD_RAND_MAX (0x7FFF)

/* get random number with more than 32-bit precision (uniform) */
static double rand_63_precision(void) {
  // 0x7FFF: RAND_MAX guaranteed by C standard.
  const int64_t b0 = rand() & CSTD_RAND_MAX;
  const int64_t b1 = rand() & CSTD_RAND_MAX;
  const int64_t b2 = rand() & CSTD_RAND_MAX;
  const int64_t b3 = rand() & CSTD_RAND_MAX;
  const int64_t b4 = rand() & 0x7;
  const int64_t n = (b0) | (b1 << 15) | (b2 << 30) | (b3 << 45) | (b4 << 60);
  return (double)n/(double)INT64_MAX;
}

/* get array with shuffled indices */
static int64_t* permutation(int64_t len){
  int64_t* array;
  array = malloc(len * sizeof(int64_t));
  int64_t i,j,temp;
  for (i=0;i<len;i++) array[i] = i;
  for (i=0;i<len;i++){
    j = i + (int64_t)((len-1-i) * rand_63_precision());
    temp = array[i];
    array[i] = array[j];
    array[j] = temp;
  }
  return array;
}

/* find root and path compression */
static int64_t findroot(Graph* g, int64_t i){
  if (g->ptr[i]<0) return i;  // return index of root node
  return g->ptr[i] = findroot(g, g->ptr[i]);  // recursively go to root node, plus: do path-compression on the fly
}

/* merge two graph fragments in g->ptr representation, return new root node index */
static int64_t merge_root(Graph* g, int64_t r1, int64_t r2){
  if (g->parity[r1] && g->parity[r2]) g->num_parity--;
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

/* apply erasures and Pauli errors, compute syndromes at same time */
int apply_erasure_and_error(Graph* g, double p_erasure, double p_error){
  memset(g->erasure, 0, g->nnode * sizeof(bool));
  memset(g->syndrome, 0, g->nnode * sizeof(bool));
  memset(g->error, 0, g->nnode * sizeof(bool));
  int num_syndromes = 0;
  for(int i=0; i < g->nnode; i++){
    if(g->is_qbt[i]){
      g->erasure[i] = (rand_63_precision() < p_erasure);
      if ((g->erasure[i] && rand_63_precision() < 0.5) || (rand_63_precision() < p_error)) {
        g->error[i] = 1;
        g->syndrome[g->nn[i*g->num_nb_max]] = (g->syndrome[g->nn[i*g->num_nb_max]] + 1) % 2; // 1st syndrome neighbor
        num_syndromes += (g->syndrome[g->nn[i*g->num_nb_max]] ? -1 : +1);
        g->syndrome[g->nn[i*g->num_nb_max+1]] = (g->syndrome[g->nn[i*g->num_nb_max+1]] + 1) % 2; // 2nd syndrome neighbor
        num_syndromes += (g->syndrome[g->nn[i*g->num_nb_max+1]] ? -1 : +1);
      }
    }
  }
  memcpy(g->parity, g->syndrome, g->nnode * sizeof(bool)); // syndrome and parity of cluster starts as the same thing (when all nodes are isolated)
  return num_syndromes;
}

/* get clusters with even number of syndromes by breadth-first traversal */
void get_even_clusters_bfs(Graph* g, int num_syndromes){
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
  int bfs_next = 0; // next free position in bfs_list
  while(g->num_parity > 0 || bfs_pos < num_seed){
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
}


