#ifndef SPKMEANS_H_
#define SPKMEANS_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*constants*/
#define CONVERGENCE_EPSILON pow(10, -5)
#define MAX_ITER_JACOBI 100
#define NUMBER_OF_VALID_ENUM_OPTIONS 4

/*structs*/
/*struct that represents a matrix of doubles, including it's dimensions*/
typedef struct
{
    double **values;
    int rows;
    int columns;
} matrix;

/*struct for the output of the jacobi algorithm.
the struct will conatin:
1. eigenvalues - the diagonal of A', represented as one dimeantional matrix (1 row).
2. eigenvectors - V matrix.*/
typedef struct
{
    matrix *eigenvalues;
    matrix *eigenvectors_mat;
} eigen_data;

/*struct for connection between eigenvalue to it's eigenvector (which is represented as one dimeantional matrix (1 row))*/
typedef struct
{
    matrix *eigenvector;
    double eigenvalue;
} eigen_pair;

/*strcut cluster, representing a cluster of the kmeans++ algorithm.*/
typedef struct
{
    int *points_indexes;
    int num_points;
} cluster;

/*enum representing the output options*/
typedef enum
{
    wam,
    ddg,
    lnorm,
    jacobi,
    invalid
} GOAL;

/*for conversion function, from https://stackoverflow.com/questions/16844728/converting-from-string-to-enum-in-c*/
static const struct
{
    GOAL val;
    const char *str;
} conversion[] = {{wam, "wam"}, {ddg, "ddg"}, {lnorm, "lnorm"}, {jacobi, "jacobi"}};

/*declaretion - general*/
GOAL string_to_goal(const char *str);

/*declarations - The Normalized Spectral Clustering Algorithm*/
double euclidean_norm_squared(double *point1, double *point2, int d);
matrix *create_matrix(int rows, int columns);
void print_matrix(matrix *mat);
void print_matrix_formatted(matrix *mat, int format_zero);
matrix *dot_product(matrix *mat1, matrix *mat2);
matrix *calculate_w(matrix *data_points);
matrix *calculate_d_from_w(matrix *w_mat);
matrix *calculate_d(matrix *data_points);
void sqrt_diagonal(matrix *mat);
matrix *calculate_l(matrix *data_points);
void normalize_matrix(matrix *mat);
matrix *create_identity_matrix(int n);
eigen_data *calculate_jacobi(matrix *a);
matrix *calculate_p(matrix *a, int i, int j);
void find_pivot(matrix *a, int *i, int *j);
void rotate_matrix(matrix *a, matrix *p, int i, int j);
double off(matrix *a);
int compare_eigen_pair(const void *p1, const void *p2);
eigen_pair *extract_sorted_eigen_pairs(eigen_data *jacobi_data);
int find_k_heuristic(eigen_pair *eigen_pairs, int n);
matrix *calculate_t(eigen_pair *eigen_pairs, int n, int k);
matrix *fit_nsc_c(matrix *data_points, int k);
void free_matrix(matrix *mat);
void free_eigen_data(eigen_data *jacobi_data);
void free_eigen_pairs(eigen_pair *eigen_pairs, int length);

/*declarations - kmeans++ algorithm*/
matrix *read_input_file(char *input_file_name);
matrix *fit_kmeanspp_c(matrix *data_points, matrix *centroids, int max_iter, double epsilon);
void reset_clusters(cluster *clusters, int k);
cluster *create_clusters(int k, int n);
int find_closest_centroid(double *point, matrix *centroids);
void add_point(cluster *cluster, int point_index);
int update_centroid(cluster *current_cluster, double *current_centroid, double *new_centroid, matrix *data_points, double epsilon);
void free_clusters(cluster *clusters, int k);
int invalid_input();
int error();

#endif
