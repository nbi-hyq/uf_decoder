#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/decoder_main.h"

/* create half-Tanner-graph for square lattice toric code (lsize: length in one dimension) */
Graph get_2d_toric_code(int lsize){
  uint8_t num_nb_max = 4;
  int n_qbt = 2*lsize*lsize;
  int n_syndrome = lsize*lsize;
  Graph g = new_graph(n_qbt + n_syndrome, num_nb_max);
  g.num_crr_x = lsize;
  g.num_crr_y = lsize;
  g.crr_surf_x = malloc(g.num_crr_x * sizeof(int));
  g.crr_surf_y = malloc(g.num_crr_y * sizeof(int));
  int cnt_crr_x = 0;
  int cnt_crr_y = 0;
  for(int x=0; x<lsize; x++){
    for(int y=0; y<lsize; y++){
      /* define correlation surfaces */
      if(x == 0) g.crr_surf_x[cnt_crr_x++] = 3*(y*lsize+x)+1;
      if(y == 0) g.crr_surf_y[cnt_crr_y++] = 3*(y*lsize+x)+2;
      /* syndrome node */
      g.nn[3*(y*lsize+x)*num_nb_max] = 3*(y*lsize+x) + 1;
      g.nn[3*(y*lsize+x)*num_nb_max + 1] = 3*(y*lsize+x) + 2;
      if(x>0){
        g.nn[3*(y*lsize+x)*num_nb_max + 2] = 3*(y*lsize+x-1) + 1;
      } else {
        g.nn[3*(y*lsize+x)*num_nb_max + 2] = 3*(y*lsize+lsize-1) + 1;
      }
      if(y>0){
        g.nn[3*(y*lsize+x)*num_nb_max + 3] = 3*((y-1)*lsize+x) + 2;
      } else{
        g.nn[3*(y*lsize+x)*num_nb_max + 3] = 3*((lsize-1)*lsize+x) + 2;
      }
      g.len_nb[3*(y*lsize+x)] = num_nb_max;
      g.is_qbt[3*(y*lsize+x)] = 0; // syndrome
      /* qubit nodes */
      g.nn[(3*(y*lsize+x)+1)*num_nb_max] = 3*(y*lsize+x);
      if(x+1 < lsize){
        g.nn[(3*(y*lsize+x)+1)*num_nb_max + 1] = 3*(y*lsize+x+1);
      } else {
        g.nn[(3*(y*lsize+x)+1)*num_nb_max + 1] = 3*(y*lsize+0);
      }
      g.len_nb[3*(y*lsize+x)+1] = 2;
      g.is_qbt[3*(y*lsize+x)+1] = 1; // qubit
      g.nn[(3*(y*lsize+x)+2)*num_nb_max] = 3*(y*lsize+x);
      if(y+1 < lsize){
        g.nn[(3*(y*lsize+x)+2)*num_nb_max + 1] = 3*((y+1)*lsize+x);
      } else {
        g.nn[(3*(y*lsize+x)+2)*num_nb_max + 1] = 3*x;
      }
      g.len_nb[3*(y*lsize+x)+2] = 2;
      g.is_qbt[3*(y*lsize+x)+2] = 1; // qubit
    }
  }
  return g;
}

/* create half-Tanner-graph for triangular lattice toric code (lsize: length in one dimension) */
Graph get_2d_triangular_toric_code(int lsize){
  uint8_t num_nb_max = 6;
  int n_qbt = 3*lsize*lsize;
  int n_syndrome = lsize*lsize;
  Graph g = new_graph(n_qbt + n_syndrome, num_nb_max);
  g.num_crr_x = 2*lsize;
  g.num_crr_y = 2*lsize;
  g.crr_surf_x = malloc(g.num_crr_x * sizeof(int));
  g.crr_surf_y = malloc(g.num_crr_y * sizeof(int));
  int cnt_crr_x = 0;
  int cnt_crr_y = 0;
  for(int x=0; x<lsize; x++){
    for(int y=0; y<lsize; y++){
      /* define correlation surfaces */
      if(x == 0){
        g.crr_surf_x[cnt_crr_x++] = 4*(y*lsize+x)+1;
        g.crr_surf_x[cnt_crr_x++] = 4*(y*lsize+x)+3; // diagonal connection
      }
      if(y == 0){
        g.crr_surf_y[cnt_crr_y++] = 4*(y*lsize+x)+2;
        g.crr_surf_y[cnt_crr_y++] = 4*(y*lsize+x)+3; // diagonal connection
      }
      /* syndrome node */
      g.nn[4*(y*lsize+x)*num_nb_max] = 4*(y*lsize+x) + 1;
      g.nn[4*(y*lsize+x)*num_nb_max + 1] = 4*(y*lsize+x) + 2;
      g.nn[4*(y*lsize+x)*num_nb_max + 2] = 4*(y*lsize+x) + 3;
      if(x>0){
        g.nn[4*(y*lsize+x)*num_nb_max + 3] = 4*(y*lsize+x-1) + 1;
      } else {
        g.nn[4*(y*lsize+x)*num_nb_max + 3] = 4*(y*lsize+lsize-1) + 1;
      }
      if(y>0){
        g.nn[4*(y*lsize+x)*num_nb_max + 4] = 4*((y-1)*lsize+x) + 2;
        if(x>0){
          g.nn[4*(y*lsize+x)*num_nb_max + 5] = 4*((y-1)*lsize+x-1) + 3;
        } else {
          g.nn[4*(y*lsize+x)*num_nb_max + 5] = 4*((y-1)*lsize+lsize-1) + 3;
        }
      } else{
        g.nn[4*(y*lsize+x)*num_nb_max + 4] = 4*((lsize-1)*lsize+x) + 2;
        if(x>0){
          g.nn[4*(y*lsize+x)*num_nb_max + 5] = 4*((lsize-1)*lsize+x-1) + 3;
        } else {
          g.nn[4*(y*lsize+x)*num_nb_max + 5] = 4*((lsize-1)*lsize+lsize-1) + 3;
        }
      }
      g.len_nb[4*(y*lsize+x)] = num_nb_max;
      g.is_qbt[4*(y*lsize+x)] = 0; // syndrome
      /* qubit nodes */
      g.nn[(4*(y*lsize+x)+1)*num_nb_max] = 4*(y*lsize+x);
      if(x+1 < lsize){
        g.nn[(4*(y*lsize+x)+1)*num_nb_max + 1] = 4*(y*lsize+x+1);
      } else {
        g.nn[(4*(y*lsize+x)+1)*num_nb_max + 1] = 4*(y*lsize+0);
      }
      g.len_nb[4*(y*lsize+x)+1] = 2;
      g.is_qbt[4*(y*lsize+x)+1] = 1; // qubit
      g.nn[(4*(y*lsize+x)+2)*num_nb_max] = 4*(y*lsize+x);
      if(y+1 < lsize){
        g.nn[(4*(y*lsize+x)+2)*num_nb_max + 1] = 4*((y+1)*lsize+x);
      } else {
        g.nn[(4*(y*lsize+x)+2)*num_nb_max + 1] = 4*x;
      }
      g.len_nb[4*(y*lsize+x)+2] = 2;
      g.is_qbt[4*(y*lsize+x)+2] = 1; // qubit

      g.nn[(4*(y*lsize+x)+3)*num_nb_max] = 4*(y*lsize+x);
      if(y+1 < lsize){
        if(x+1 < lsize){
          g.nn[(4*(y*lsize+x)+3)*num_nb_max + 1] = 4*((y+1)*lsize+x+1);
        } else {
          g.nn[(4*(y*lsize+x)+3)*num_nb_max + 1] = 4*((y+1)*lsize+0);
        }
      } else {
        if(x+1 < lsize){
          g.nn[(4*(y*lsize+x)+3)*num_nb_max + 1] = 4*(x+1);
        } else {
          g.nn[(4*(y*lsize+x)+3)*num_nb_max + 1] = 0;
        }
      }
      g.len_nb[4*(y*lsize+x)+3] = 2;
      g.is_qbt[4*(y*lsize+x)+3] = 1; // qubit
    }
  }
  return g;
}

/* make a few checks if graph can be correct */
int validate_graph(Graph* g){
  for(int n=0; n<g->nnode; n++){
    if(g->is_qbt[n] && g->len_nb[n] != 2) return 1;
    for(uint8_t k=0; k<g->len_nb[n]; k++){
      int nb = g->nn[n*g->num_nb_max+k];
      if(nb == n) return 2; // no self-edges
      if(g->is_qbt[n] == g->is_qbt[nb]) return 3; // qubit connected to syndrome (bi-partite)
      bool unidir = true;
      for(uint8_t j=0; j<g->len_nb[nb]; j++){
        if(g->nn[nb*g->num_nb_max+j] == n){
          unidir = false;
          break;
        }
      }
      if(unidir) return 4; // connection is not bi-directional
    }
  }
  return 0;
}

/* visualize errors of square lattice surface code */
int visualize_error(Graph* g, int size){
  for(int y=0; y<size; y++){
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(g->syndrome[3*n]) printf("x");
      else printf("o");
      if(g->error[3*n+1]) printf("e");
      else printf("-");
    }
    printf("\n");
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(g->error[3*n+2]) printf("e");
      else printf("|");
      printf(" ");
    }
    printf("\n");
  }
  printf("--------------------------\n");
  return 0;
}

/* visualize decoding of square lattice surface code */
int visualize_decode(Graph* g, int size){
  for(int y=0; y<size; y++){
    for(int x=0; x<size; x++){
      int n = x + y*size;
      printf("o");
      if(!g->erasure[3*n+1]){
        if(g->error[3*n+1] && !g->decode[3*n+1]) printf("e");
        else if(!g->error[3*n+1] && g->decode[3*n+1]) printf("c");
        else if(g->error[3*n+1] && g->decode[3*n+1]) printf("b");
        else printf("-");
      } else {
        if(g->error[3*n+1] && !g->decode[3*n+1]) printf("E");
        else if(!g->error[3*n+1] && g->decode[3*n+1]) printf("C");
        else if(g->error[3*n+1] && g->decode[3*n+1]) printf("B");
        else printf("=");
      }
    }
    printf("\n");
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(!g->erasure[3*n+2]){
        if(g->error[3*n+2] && !g->decode[3*n+2]) printf("e");
        else if(!g->error[3*n+2] && g->decode[3*n+2]) printf("c");
        else if(g->error[3*n+2] && g->decode[3*n+2]) printf("b");
        else printf("|");
      } else {
        if(g->error[3*n+2] && !g->decode[3*n+2]) printf("E");
        else if(!g->error[3*n+2] && g->decode[3*n+2]) printf("C");
        else if(g->error[3*n+2] && g->decode[3*n+2]) printf("B");
        else printf("/");
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("--------------------------\n");
  return 0;
}

/* visualize syndrome validation forest */
int visualize_forest(Forest* f, int size){
  for(int y=0; y<size; y++){
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(f->leaf[3*n] && f->visited[3*n]) printf("l");
      else if(f->visited[3*n]) printf("v");
      else printf(".");
      if(f->leaf[3*n+1] && f->visited[3*n+1]) printf("l");
      else if(f->visited[3*n+1]) printf("v");
      else printf(".");
    }
    printf("\n");
    for(int x=0; x<size; x++){
      int n = x + y*size;
      if(f->leaf[3*n+2] && f->visited[3*n+2]) printf("l");
      else if(f->visited[3*n+2]) printf("v");
      else printf(".");
      printf(" ");
    }
    printf("\n");
  }
  printf("--------------------------\n");
  return 0;
}

