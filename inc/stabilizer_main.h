#ifndef STABMAIN
#define STABMAIN

#include <stdint.h>
#include <stdbool.h>

#define dtypeBlk uint32_t  // other data type may be better depending on average cluster size

typedef struct {
  dtypeBlk* mat;
  uint64_t num_blocks_row;
  uint64_t num_bytes;
  uint64_t num_bytes_row;
  int n_row; // number of columns of active W-matrix
  int n_col; // number of rows of active W-matrix
  int n_row_non_zero; // the first non-zero columns
  int idx_start; // row index where the active matrix starts
  int nqbt;
} WMat;

WMat new_wmat(int n_row, int n_col);
void free_wmat(WMat w);

void GaussElimin_bit(dtypeBlk*, int, int, int);
void print_bitwise(dtypeBlk* row, uint64_t num_blocks);
void print_bitwise_matrix(dtypeBlk* mat, int n_rows, int n_cols);
void write_matrix_position_bit(dtypeBlk* mat, int i_r, int i_c, uint64_t num_blocks_row, bool bit);
bool read_matrix_position_bit(dtypeBlk* mat, int i_r, int i_c, uint64_t num_blocks_row);
void flip_matrix_position_bit(dtypeBlk* mat, int i_r, int i_c, uint64_t num_blocks_row);
void GaussElimin_bit(dtypeBlk* mat, int n_rows, int n_cols, int num_eliminated_col);
void GaussElimin_bit_col_wise(dtypeBlk* mat, int n_rows, int n_cols, int num_A, int* col_A);
void GetWmat(WMat w, int num_edges, int** edges);
void GetWmat_from_bool_array(WMat* w, bool* bool_array, int n_row);
void GetWmat_from_adj(WMat* w, bool* adj, int n_row);

#endif
