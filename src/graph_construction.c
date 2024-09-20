#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/decoder_main.h"

/* Define Tanner graphs of toric codes, visualization of 2d square-lattice toric code decoding. Note that data qubits must come before syndrome qubits in global indexing */

/* create half-Tanner-graph for square lattice toric code (lsize: length in one dimension) */
Graph get_2d_toric_code(int lsize){
  uint8_t num_nb_max_qbt = 2;
  uint8_t num_nb_max_syndr = 4;
  int n_qbt = 2*lsize*lsize;
  int n_syndrome = lsize*lsize;
  Graph g = new_graph(n_qbt, num_nb_max_qbt, n_syndrome, num_nb_max_syndr);
  g.num_logicals = 2;
  g.logicals = malloc(g.num_logicals * sizeof(int*));
  g.logical_weight = malloc(g.num_logicals * sizeof(int));
  for(int i=0; i<g.num_logicals; i++){
    g.logical_weight[i] = lsize;
    g.logicals[i] = malloc(g.logical_weight[i] * sizeof(int));
  }
  int cnt_x = 0;
  int cnt_y = 0;
  for(int x=0; x<lsize; x++){
    for(int y=0; y<lsize; y++){
      /* define logicals */
      if(x == 0) g.logicals[0][cnt_x++] = 2*(y*lsize+x)+1;
      if(y == 0) g.logicals[1][cnt_y++] = 2*(y*lsize+x)+0;
      /* syndromes */
      g.nn_syndr[(y*lsize+x)*num_nb_max_syndr] = 2*(y*lsize+x) + 0;
      g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 1] = 2*(y*lsize+x) + 1;
      if(x<lsize-1){
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 2] = 2*(y*lsize+x+1) + 1;
      } else {
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 2] = 2*(y*lsize+0) + 1;
      }
      if(y<lsize-1){
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 3] = 2*((y+1)*lsize+x) + 0;
      } else{
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 3] = 2*(0*lsize+x) + 0;
      }
      g.len_nb[n_qbt + y*lsize+x] = num_nb_max_syndr;

      /* data qubits */
      g.nn_qbt[(2*(y*lsize+x)+0)*num_nb_max_qbt] = y*lsize+x + n_qbt;
      if(y>0){
        g.nn_qbt[(2*(y*lsize+x)+0)*num_nb_max_qbt + 1] = (y-1)*lsize+x + n_qbt;
      } else {
        g.nn_qbt[(2*(y*lsize+x)+0)*num_nb_max_qbt + 1] = (lsize-1)*lsize+x + n_qbt;
      }
      g.len_nb[2*(y*lsize+x)+0] = num_nb_max_qbt;
      g.nn_qbt[(2*(y*lsize+x)+1)*num_nb_max_qbt] = y*lsize+x + n_qbt;
      if(x>0){
        g.nn_qbt[(2*(y*lsize+x)+1)*num_nb_max_qbt + 1] = y*lsize+(x-1) + n_qbt;
      } else {
        g.nn_qbt[(2*(y*lsize+x)+1)*num_nb_max_qbt + 1] = y*lsize+(lsize-1) + n_qbt;
      }
      g.len_nb[2*(y*lsize+x)+1] = num_nb_max_qbt;
    }
  }
  return g;
}

/* create half-Tanner-graph for triangular lattice toric code (lsize: length in one dimension) */
Graph get_2d_triangular_toric_code(int lsize){
  uint8_t num_nb_max_qbt = 2;
  uint8_t num_nb_max_syndr = 6;
  int n_qbt = 3*lsize*lsize;
  int n_syndrome = lsize*lsize;
  Graph g = new_graph(n_qbt, num_nb_max_qbt, n_syndrome, num_nb_max_syndr);
  g.num_logicals = 2;
  g.logicals = malloc(g.num_logicals * sizeof(int*));
  g.logical_weight = malloc(g.num_logicals * sizeof(int));
  for(int i=0; i<g.num_logicals; i++){
    g.logical_weight[i] = 2*lsize;
    g.logicals[i] = malloc(g.logical_weight[i] * sizeof(int));
  }
  int cnt_x = 0;
  int cnt_y = 0;
  for(int x=0; x<lsize; x++){
    for(int y=0; y<lsize; y++){
      /* define logicals */
      if(x == 0){
        g.logicals[0][cnt_x++] = 3*(y*lsize+x)+0;
        g.logicals[0][cnt_x++] = 3*(y*lsize+x)+2;
      }
      if(y == 0){
        g.logicals[1][cnt_y++] = 3*(y*lsize+x)+1;
        g.logicals[1][cnt_y++] = 3*(y*lsize+x)+2;
      }
      /* syndromes */
      g.nn_syndr[(y*lsize+x)*num_nb_max_syndr] = 3*(y*lsize+x) + 0;
      g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 1] = 3*(y*lsize+x) + 1;
      g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 2] = 3*(y*lsize+x) + 2;
      if(x>0){
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 3] = 3*(y*lsize+x-1) + 0;
      } else {
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 3] = 3*(y*lsize+lsize-1) + 0;
      }
      if(y>0){
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 4] = 3*((y-1)*lsize+x) + 1;
        if(x>0) g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 5] = 3*((y-1)*lsize+x-1) + 2;
        else g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 5] = 3*((y-1)*lsize+lsize-1) + 2;
      } else{
        g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 4] = 3*((lsize-1)*lsize+x) + 1;
        if(x>0) g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 5] = 3*((lsize-1)*lsize+x-1) + 2;
        else g.nn_syndr[(y*lsize+x)*num_nb_max_syndr + 5] = 3*((lsize-1)*lsize+lsize-1) + 2;
      }
      g.len_nb[n_qbt + y*lsize+x] = num_nb_max_syndr;

      /* data qubits */
      g.len_nb[3*(y*lsize+x)+0] = num_nb_max_qbt;
      g.len_nb[3*(y*lsize+x)+1] = num_nb_max_qbt;
      g.len_nb[3*(y*lsize+x)+2] = num_nb_max_qbt;
      g.nn_qbt[(3*(y*lsize+x)+0)*num_nb_max_qbt] = y*lsize+x + n_qbt;
      g.nn_qbt[(3*(y*lsize+x)+1)*num_nb_max_qbt] = y*lsize+x + n_qbt;
      g.nn_qbt[(3*(y*lsize+x)+2)*num_nb_max_qbt] = y*lsize+x + n_qbt;
      if(x<lsize-1){
        g.nn_qbt[(3*(y*lsize+x)+0)*num_nb_max_qbt + 1] = y*lsize+x+1 + n_qbt;
      } else {
        g.nn_qbt[(3*(y*lsize+x)+0)*num_nb_max_qbt + 1] = y*lsize+0 + n_qbt;
      }
      if(y<lsize-1){
        g.nn_qbt[(3*(y*lsize+x)+1)*num_nb_max_qbt + 1] = (y+1)*lsize+x + n_qbt;
        if(x<lsize-1) g.nn_qbt[(3*(y*lsize+x)+2)*num_nb_max_qbt + 1] = (y+1)*lsize+x+1 + n_qbt;
        else g.nn_qbt[(3*(y*lsize+x)+2)*num_nb_max_qbt + 1] = (y+1)*lsize+0 + n_qbt;
      } else {
        g.nn_qbt[(3*(y*lsize+x)+1)*num_nb_max_qbt + 1] = 0*lsize+x + n_qbt;
        if(x<lsize-1) g.nn_qbt[(3*(y*lsize+x)+2)*num_nb_max_qbt + 1] = 0*lsize+x+1 + n_qbt;
        else g.nn_qbt[(3*(y*lsize+x)+2)*num_nb_max_qbt + 1] = 0*lsize+0 + n_qbt;
      }
    }
  }
  return g;
}

/* make a few checks if Tanner graph can be correct */
int validate_graph(Graph* g){
  for(int n=0; n<g->n_qbt; n++){
    if(g->len_nb[n] != 2) return 1; // no topological code
    for(uint8_t k=0; k<g->len_nb[n]; k++){
      int nb = g->nn_qbt[n*g->num_nb_max_qbt+k];
      if(nb == n) return 2; // self-edge
      if(nb < g->n_qbt) return 3; // qubit not connected to syndrome
      bool unidir = true;
      for(uint8_t j=0; j<g->len_nb[nb]; j++){
        if(g->nn_syndr[(nb - g->n_qbt)*g->num_nb_max_syndr+j] == n){
          unidir = false;
          break;
        }
      }
      if(unidir) return 4; // connection is not bi-directional
    }
  }
  for(int n=g->n_qbt; n<g->n_qbt+g->n_syndr; n++){
    for(uint8_t k=0; k<g->len_nb[n]; k++){
      int nb = g->nn_syndr[(n - g->n_qbt)*g->num_nb_max_syndr+k];
      if(nb == n) return 2; // self-edge
      if(nb >= g->n_qbt) return 3; // syndrome not connected to data qubit
      bool unidir = true;
      for(uint8_t j=0; j<g->len_nb[nb]; j++){
        if(g->nn_qbt[nb*g->num_nb_max_qbt+j] == n){
          unidir = false;
          break;
        }
      }
      if(unidir) return 4; // connection is not bi-directional
    }
  }
  return 0;
}

/* visualize errors of square lattice 2d toric code */
int visualize_error(Graph* g, int size){
  for(int y=size-1; y>=0; y--){
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(!g->erasure[2*n+1]){
        if(g->error[2*n+1]) printf("e");
        else printf("-");
      } else {
        if(g->error[2*n+1]) printf("E");
        else printf("=");
      }
      if(g->syndrome[n]) printf("x");
      else printf("o");
    }
    printf("\n");
    for(int x=0; x<size; x++){
      printf(" ");
      int n = x + y*size;
      if(!g->erasure[2*n+0]){
        if(g->error[2*n+0]) printf("e");
        else printf("|");
      } else {
        if(g->error[2*n+0]) printf("E");
        else printf("/");
      }
    }
    printf("\n");
  }
  printf("--------------------------\n");
  return 0;
}

/* visualize decoding of square lattice 2d toric code */
int visualize_decode(Graph* g, int size){
  for(int y=size-1; y>=0; y--){
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(!g->erasure[2*n+1]){
        if(g->error[2*n+1] && !g->decode[2*n+1]) printf("e");
        else if(!g->error[2*n+1] && g->decode[2*n+1]) printf("c");
        else if(g->error[2*n+1] && g->decode[2*n+1]) printf("b");
        else printf("-");
      } else {
        if(g->error[2*n+1] && !g->decode[2*n+1]) printf("E");
        else if(!g->error[2*n+1] && g->decode[2*n+1]) printf("C");
        else if(g->error[2*n+1] && g->decode[2*n+1]) printf("B");
        else printf("=");
      }
      printf("o");
    }
    printf("\n");
    for(int x=0; x<size; x++){
      printf(" ");
      int n = x + y*size;
      if(!g->erasure[2*n]){
        if(g->error[2*n] && !g->decode[2*n]) printf("e");
        else if(!g->error[2*n] && g->decode[2*n]) printf("c");
        else if(g->error[2*n] && g->decode[2*n]) printf("b");
        else printf("|");
      } else {
        if(g->error[2*n] && !g->decode[2*n]) printf("E");
        else if(!g->error[2*n] && g->decode[2*n]) printf("C");
        else if(g->error[2*n] && g->decode[2*n]) printf("B");
        else printf("/");
      }
    }
    printf("\n");
  }
  printf("--------------------------\n");
  return 0;
}

/* visualize syndrome validation forest (for square lattice 2d toric code) */
int visualize_forest(Forest* f, int size){
  for(int y=size-1; y>=0; y--){
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(f->leaf[2*n+1] && f->visited[2*n+1]) printf("l"); // data qubit
      else if(f->visited[2*n+1]) printf("v");
      else printf(".");
      if(f->leaf[2*size*size + n] && f->visited[2*size*size + n]) printf("l"); // syndrome
      else if(f->visited[2*size*size + n]) printf("v");
      else printf(".");
    }
    printf("\n");
    for(int x=0; x<size; x++){
      printf(" ");
      int n = x + y*size;
      if(f->leaf[2*n] && f->visited[2*n]) printf("l"); // data qubit
      else if(f->visited[2*n]) printf("v");
      else printf(".");
    }
    printf("\n");
  }
  printf("--------------------------\n");
  return 0;
}

