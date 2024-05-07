#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../inc/stabilizer_main.h"

/* Gaussian elimination used for decoding/validating clusters
   derived from code by: Matthias C. Loebl, Love A. Pettersson (2023-2024) */

/* create new Wmat struct */
WMat new_wmat(int n_row, int n_col){
  WMat w;
  w.num_blocks_row = n_col / (8 * sizeof(dtypeBlk)) + 1;
  w.num_bytes_row = w.num_blocks_row * sizeof(dtypeBlk);
  w.num_bytes = w.num_bytes_row * n_row;
  w.mat = (dtypeBlk*)malloc(w.num_bytes);
  w.n_row = n_row; // number of rows of active W-matrix
  w.n_col = n_col; // number of columns of active W-matrix
  w.n_row_non_zero = n_row; // (only for ker/img) the first non-zero columns
  w.idx_start = 0; // row index where the active matrix starts
  w.nqbt = n_row; // number of qubits
  return w;
}

/* free everything related to WMat struct */
void free_wmat(WMat w){
  free(w.mat);
}

/* print a row from a matrix to the terminal (with bitwise logic) */
void print_bitwise(dtypeBlk* row, uint64_t num_blocks){
  for(uint64_t i=0; i<num_blocks; i++){
    for(uint8_t byte=0; byte<sizeof(dtypeBlk); byte++){
      for(uint8_t bit=0; bit<8; bit++){
        printf("%i", (row[i] >> (byte*8+bit)) & 1ul ? 1 : 0);
      }
      printf(",");
    }
    printf("|");
  }
  printf("\n");
}

/* print matrix mat to terminal (with bitwise logic) */
void print_bitwise_matrix(dtypeBlk* mat, int n_rows, int n_cols){
  uint64_t num_blocks = n_cols / (8 * sizeof(dtypeBlk)) + 1;
  for(int r=0; r<n_rows; r++){
    print_bitwise(mat + r * num_blocks, num_blocks);
  }
  printf("------------------------\n");
}

/* in the bitwise matrix "mat" swap "row0, row1", "num_blocks_row" is the number of elements of "dtypeBlk" that represent one row, "num_bytes_row" is the corresponding number in bytes */
static void SwitchRows_bitwise(dtypeBlk* mat, dtypeBlk* dummy, uint64_t num_bytes_row, uint64_t num_blocks_row, int row0, int row1){
  dtypeBlk* p0 = mat + row0 * num_blocks_row;
  dtypeBlk* p1 = mat + row1 * num_blocks_row;
  memcpy(dummy, p0, num_bytes_row);
  memcpy(p0, p1, num_bytes_row);
  memcpy(p1, dummy, num_bytes_row);
}

/* in the bitwise matrix "mat" subtract "row0" from "row_target", "num_blocks_row" is the number of elements of "dtypeBlk" that represent one row */
static void SubtractRows_bitwise(dtypeBlk* mat, uint64_t num_blocks_row, int row0, int row_target){
  dtypeBlk* p0 = mat + row0 * num_blocks_row;
  dtypeBlk* ptarg = mat + row_target * num_blocks_row;
  for (uint64_t i = 0; i < num_blocks_row; i++) {
    *(ptarg + i) = *(ptarg + i) ^ *(p0 + i);
  }
}

/* write a single bit of matrix "mat" at position "i_r, i_c", "num_blocks_row" is the number of elements of "dtypeBlk" that represent one row */
void write_matrix_position_bit(dtypeBlk* mat, int i_r, int i_c, uint64_t num_blocks_row, bool bit){
  uint64_t pos_block = i_r * num_blocks_row + i_c / (8 * sizeof(dtypeBlk));
  uint64_t pos_bit = i_c % (8 * sizeof(dtypeBlk));
  if (bit){
    uint64_t a = 1ul << pos_bit;
    mat[pos_block] |= a;
  }
  else {
    uint64_t b = 1ul << pos_bit;
    uint64_t a = 0ul - 1;
    a -= b;
    mat[pos_block] &= a;
  }
}

/* read a single bit of matrix "mat" at position "i_r, i_c", "num_blocks_row" is the number of elements of "dtypeBlk" that represent one row */
bool read_matrix_position_bit(dtypeBlk* mat, int i_r, int i_c, uint64_t num_blocks_row){
  uint64_t pos_block = i_r * num_blocks_row + i_c / (8 * sizeof(dtypeBlk));
  uint64_t pos_bit = i_c % (8 * sizeof(dtypeBlk));
  if((mat[pos_block] >> pos_bit) & 1ul){
    return 1;
  }else{
    return 0;
  }
}

/* flip a single bit of matrix "mat" at position "i_r, i_c", "num_blocks_row" is the number of elements of "dtypeBlk" that represent one row */
void flip_matrix_position_bit(dtypeBlk* mat, int i_r, int i_c, uint64_t num_blocks_row){
  uint64_t pos_block = i_r * num_blocks_row + i_c / (8 * sizeof(dtypeBlk));
  uint64_t pos_bit = i_c % (8 * sizeof(dtypeBlk));
  uint64_t a = 1ul << pos_bit;
  mat[pos_block] ^= a;
}

/* Gaussian elimination of "mat" in F2 with a matrix element represented by on single bit
   n_rows, n_cols: dimensions of the matrix
   num_eliminated_col: terminates when that number of colums fulfill the Gaussian elimination condition */
void GaussElimin_bit(dtypeBlk* mat, int n_rows, int n_cols, int num_eliminated_col) {
  int ix_r, lead, i;
  uint64_t num_blocks_row = n_cols / (8 * sizeof(dtypeBlk)) + 1;
  uint64_t num_bytes_row = num_blocks_row * sizeof(dtypeBlk);
  dtypeBlk* dummy = (dtypeBlk*)malloc(num_bytes_row);
  lead = 0;
  for (ix_r = 0; ix_r < (n_rows - 1); ix_r++) {
    if (lead >= num_eliminated_col) {
      free(dummy);
      return;
    }
    i = ix_r;
    while (read_matrix_position_bit(mat, i, lead, num_blocks_row) == 0){
      ++i;
      if (i == n_rows) {
        i = ix_r;
        ++lead;
        if (lead == num_eliminated_col) {
          free(dummy);
          return;
        }
      }
    }
    if (i != ix_r) {
      SwitchRows_bitwise(mat, dummy, num_bytes_row, num_blocks_row, ix_r, i);
    }
    for (i = ix_r + 1; i < n_rows; i++) {
      if (read_matrix_position_bit(mat, i, lead, num_blocks_row) != 0) {
        SubtractRows_bitwise(mat, num_blocks_row, ix_r, i);
      }
    }
    lead++;
  }
  free(dummy);
  return;
}

/* Gaussian elimination of "mat" in F2 with a matrix element represented by on single bit
   n_rows, n_cols: dimensions of the matrix
   col_A: array speficying the order in which the columns are considered in the Gaussian elimination
   num_A: number of elements in col_A */
void GaussElimin_bit_col_wise(dtypeBlk* mat, int n_rows, int n_cols, int num_A, int* col_A) {
  int ix_r, lead, i;
  uint64_t num_blocks_row = n_cols / (8 * sizeof(dtypeBlk)) + 1;
  uint64_t num_bytes_row = num_blocks_row * sizeof(dtypeBlk);
  dtypeBlk* dummy = (dtypeBlk*)malloc(num_bytes_row);
  lead = 0;
  for (ix_r = 0; ix_r < (n_rows - 1); ix_r++) {
    if (lead >= num_A) {
      free(dummy);
      return;
    }
    i = ix_r;
    while (read_matrix_position_bit(mat, i, col_A[lead], num_blocks_row) == 0){
      ++i;
      if (i == n_rows) {
        i = ix_r;
        ++lead;
        if (lead == num_A) {
          free(dummy);
          return;
        }
      }
    }
    if (i != ix_r) {
      SwitchRows_bitwise(mat, dummy, num_bytes_row, num_blocks_row, ix_r, i);
    }
    for (i = ix_r + 1; i < n_rows; i++) {
      if (read_matrix_position_bit(mat, i, col_A[lead], num_blocks_row) != 0) {
        SubtractRows_bitwise(mat, num_blocks_row, ix_r, i);
      }
    }
    lead++;
  }
  free(dummy);
  return;
}

/* get W = (X|Z) stabilizer matrix (binary, bitwise representation) from a graph specified by "num_edges" edges.
   int** edges has two elements edges[0] points to an array with the first nodes in the edges and edges[1] pointing to an array with the second nodes in the edges */
void GetWmat(WMat w, int num_edges, int** edges){
memset(w.mat, 0, w.num_bytes);
for (int q = 0; q < w.nqbt; q++) write_matrix_position_bit(w.mat, q, q, w.num_blocks_row, 1);
  for (int e = 0; e < num_edges; e++){
    write_matrix_position_bit(w.mat, edges[0][e], w.nqbt + edges[1][e], w.num_blocks_row, 1);
    write_matrix_position_bit(w.mat, edges[1][e], w.nqbt + edges[0][e], w.num_blocks_row, 1);
  }
}

/* bool_array is reshaped into WMat w */
void GetWmat_from_bool_array(WMat* w, bool* bool_array, int n_row){
  w->n_row = n_row; // reset number of rows if it changed
  memset(w->mat, 0, w->num_bytes);
  for (int i_r = 0; i_r < w->n_row; i_r++){
    for (int i_c = 0; i_c < w->n_col; i_c++){
      write_matrix_position_bit(w->mat, i_r, i_c, w->num_blocks_row, bool_array[i_r * w->n_col + i_c]);
    }
  }
}

/* adjacency matrix of graph state in 1d form (adj) --> WMat w */
void GetWmat_from_adj(WMat* w, bool* adj, int n_row){
  w->n_row = n_row; // reset number of rows if it changed
  memset(w->mat, 0, w->num_bytes);
  for (int i_r = 0; i_r < w->n_row; i_r++){
    write_matrix_position_bit(w->mat, i_r, i_r, w->num_blocks_row, 1); // X-part of W-matrix
    for (int i_c = w->nqbt; i_c < w->n_col; i_c++) write_matrix_position_bit(w->mat, i_r, i_c, w->num_blocks_row, adj[i_r * w->nqbt + i_c - w->nqbt]);
  }
}

