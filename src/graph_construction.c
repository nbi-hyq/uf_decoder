#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/percolation_main.h"

/* create half-Tanner-graph for square lattice toric code (lsize: length in one dimension) */
Graph get_2d_toric_code(int lsize){
  uint8_t num_nb_max = 4;
  int n_qbt = 2*lsize*lsize;
  int n_syndrome = lsize*lsize;
  Graph g = new_graph(n_qbt + n_syndrome, num_nb_max);
  for(int x=0; x<lsize; x++){
    for(int y=0; y<lsize; y++){
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
  for(int x=0; x<lsize; x++){
    for(int y=0; y<lsize; y++){
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

