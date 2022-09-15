#ifndef SOFTWARE_PROJECT_SPKMEANS_H
#define SOFTWARE_PROJECT_SPKMEANS_H

#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h>

enum Goal{
  wam_g = 1,
  ddg_g = 2,
  lnorm_g = 3,
  jacobi_g = 4,
  spk_g = 5,
  spk_g2 = 6
};


double **adjacency_matrix(double **data_points, int dimension, int len);
double calc_euclidean_norm(double *vector_1, double *vector_2, int dimension);
double **diag_mat(double **mat, int len);
double **calc_L_mat(double **diag_mat, double **adj_mat, int len);
void calc_norm_mat(double **diag_mat, int N);
double **alloc_mat(int num_rows, int num_cols);
double **calc_mul(int len, double **mat1, double **mat2);
void sab_matrix(int len, double **mat1, double **mat2);
double **calc_id_mat(int len);
double **spk_algo(double **lnorm, int len, int *K);
double **sort_matrix_values(double **mat, int len);
void eigengap_heuristic(double *eigenvalues, int len, int *K);
double **set_T(double **U_mat, int len, int K);

/* (spkmeans.c) (C) main's functions*/
double **run_goal(enum Goal goal, double **data_input, int len, int D, int *K);
void print_result(double **mat, int num_rows, int num_cols, enum Goal goal);
void msg_and_exit(int error_type, int is_error);
int get_n_d_parameters(FILE *ifp, int situation);
void set_input(FILE *ifp, double **data_input, int num_of_rows, int num_of_cols);
void free_memory(double **ArrayToFree, int num_rows);
double **calc_jacob(int len, double **A);
void matrix_copy(int num_rows, int num_cols, double **dest_mat, double **src_mat);
void A_to_A_tag(int len, double **A, int *i_pointer, int *j_pointer);
void get_params(double **A, int i, int j, double *point_1, double *point_2);
void calc_curr_P(int len, double **P_mat, int i, int j, double d1, double d2);
void get_values_from_first_mat(double *values, int len, double **mat);
void get_mat_transe(double **mat, int N);
double **miun_of_eig(int len, double *values, double **vectors);
void calc_first_mat(int len, double **mat, double d1, double d2, int i, int j, int *return_value);
double calc_off_diag(int len, double **A);


int kMeans(int len, int K, double **points, double **cent, int dimension);
int checkTheNorm(double **fresh, double **old, int dimension, int K);
int assign_cluster(double **cent, double *Datapoint, int dimension, int K);
void idcon_cent(double **fresh, double **old, int dimension, int K);
#endif
