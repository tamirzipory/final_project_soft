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
  code_of_wam = 1,
  code_of_ddg = 2,
  code_of_lnorm = 3,
  code_of_jacobi = 4,
  code_of_spk = 5,
  code_of_kmeans = 6
};


double **calc_mul(int , double **, double **);
void sab_matrix(int , double **, double **);
double **calc_id_mat(int );
double **calc_spk_method(double **, int , int *);
double **mat_sorting(double **, int );
void msg_and_exit(int , int );
int get_n_d_parameters(FILE *, int );
void func_heruicsic(double *, int , int *);
double **calc_the_T(double **, int , int );
double **run_goal(enum Goal , double **, int , int , int *);
void print_result(double **, int , int , enum Goal );
void set_input(FILE *, double **, int , int );
void free_memory(double **, int );
double **calc_jacob(int , double **);
void matrix_copy(int , int , double **, double **);
int kMeans(int , int , double **, double **, int );
int checkTheNorm(double **, double **, int , int );
int assign_cluster(double **, double *, int , int );
void idcon_cent(double **, double **, int , int );
double **mat_adj(double **, int , int );
double euc_norm_calc(double *, double *, int );
double **diag_mat(double **, int );
double **calc_L_mat(double **, double **, int );
void A_to_A_tag(int , double **, int *, int *);
void get_params(double **, int , int , double *, double *);
void calc_curr_P(int , double **, int , int , double , double );
void get_values_from_first_mat(double *, int , double **);
void get_mat_transe(double **, int );
double **miun_of_eig(int , double *, double **);
void calc_first_mat(int , double **, double , double , int , int , int *);
double calc_off_diag(int , double **);
void calc_norm_mat(double **, int );
double **alloc_mat(int , int );
#endif
