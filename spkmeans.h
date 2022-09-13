#ifndef SPKMEANSMODULE_C_SPKMEANS_H
#define SPKMEANSMODULE_C_SPKMEANS_H
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#define EPSILON 0
#define max_iter 300

/*
 * a struct for the use of the 'spk' function.
 * in order to obtain the k eigenvectors that correspond to the k-smallest eigenvalues
 */
typedef struct{
    double val;
    int index;
}PAIR;

/* Global variables */
double C, S;    /* for the jacobi algorithm step 4 (Obtain c, t) */
int index_i, index_j; /* for the jacobi algorithm step 3 (Rotation Matrix P) */
PAIR *eigen_values = NULL;  /* for the spk algorithm step 3 (Determine k and obtain the first k eigenvectors) */

/* Free allocated memory of a given 2D array*/
void free_Mem(double **array, int rows){
    int i;
    for (i = 0; i < rows; i++)
        free(array[i]);
    free(array);
}

/* Create a new double 2D matrix by dynamic memory allocation.
 * args : rows, cols - the wanted matrix dimensions.
 *        z != 0  uninitialized matrix
 *        z == 0 - initializes matrix to zeros. */
double** Array_2D(int rows, int cols, int z) {
    int i;
    double **array = (double **) malloc(rows * sizeof(double *));
    if (array == NULL)
        return NULL;

    if (z == 0) {
        for (i = 0; i < rows; i++) {
            array[i] = (double *) calloc(cols, sizeof(double));
            if (array[i] == NULL)
                return NULL;
        }
    }
    else {
        for (i = 0; i < rows; i++) {
            array[i] = (double *) malloc(cols * sizeof(double));
            if (array[i] == NULL)
                return NULL;
        }
    }
    return array;
}

/* Computes the euclidean distance between two given vectors v1, v2.
 * n - the dimension of the given vectors*/
static double euclidean_norm(double *v1, double *v2, int n) {
    double temp, totalSum = 0;
    int i;
    for(i = 0; i < n; i++){
        temp = v1[i] - v2[i];
        temp = temp * temp;
        totalSum += temp;
    }
    return sqrt(totalSum);
}

/* Util function for the k-means algorithm step 4 (Assign x_i to the closest cluster)
 * args: clusters - matrix of the clusters as rows
 *       dist - matrix to contain the min distances from the clusters
 *       v - data vector to assign
 */
static void minDistance(double **clusters, double **dist, double *v, int cols, int K) {
    int i, j, indexOfMin = 0;
    double temp = euclidean_norm(clusters[0], v, cols);
    double min = temp;
    for(i = 1; i < K; i++){
        temp = euclidean_norm(clusters[i], v, cols);
        if(temp < min){
            min = temp;
            indexOfMin = i;
        }
    }
    for(j = 0; j < cols; j++)
        dist[indexOfMin][j] += v[j];
    dist[indexOfMin][j]++;
}

/* Util function for the k-means algorithm step 5 (Update the centroids)
 * args: clusters - matrix of the clusters as rows
 *       dist - matrix to contain the min distances from the clusters
 *       K, cols - the dimension of the matrix
 */
static double update (double **dist, double **clusters, int K, int cols) {
    int i,j;
    double temp, sum,  max = 1;
    for(i = 0; i < K; i++){
        sum = dist[i][cols];
        for (j = 0; j < cols; j++)
            dist[i][j] = (dist[i][j] / sum);
        temp = euclidean_norm(clusters[i],dist[i],cols);
        if (temp > max)
            max = temp;
        for (j = 0; j < cols; j++){
            clusters[i][j] = (dist[i][j]);
            dist[i][j] = 0;
        }
        dist[i][j] = 0;
    }
    return max;
}

/*  Implantation of the k-means algorithm, as detailed in HW1.
 *   return: clusters - is a 2D matrix of the final centroids.
 */
int fit_kmeans(double **vectors, double **clusters, int K, int rows, int cols){
    int i, count = 0;
    double max = EPSILON + 1;
    double **dist = Array_2D(K, cols+1, 0);
    if (dist == NULL)
        return 1;

    while ((max > EPSILON) && (count < max_iter)){
        for (i = 0; i < rows; i++)
            minDistance(clusters, dist, vectors[i], cols, K);
        max = update(dist, clusters, K, cols);
        count++;
    }
    free_Mem(dist, K);
    return 0;
}

/*
 * Computes the Weighted Adjacency Matrix of the given data vectors.
 * args: rows, cols - dimension on the data matrix
 *       vectors - a 2D matrix of the data vectors as rows.
 * return: weighted adjacency matrix of the data vectors
 */
double** wam(double** vectors, int rows, int cols) {
    double w_ij, x;
    int i, j;
    double** W = Array_2D(rows, rows, 0);
    if (W == NULL)
        return NULL;
    for(i = 0; i < rows; i++){
        for(j = i+1; j < rows; j++){
            x = euclidean_norm(vectors[i], vectors[j], cols);
            w_ij = exp(-0.5 * x);
            W[i][j] = w_ij;
            W[j][i] = w_ij;
        }
    }
    return W;
}

/* Util function to compute the sum of a given vector element wise */
static double sumRow(double* weightedVector, int length) {
    double result = 0;
    int i;
    for (i = 0; i < length; i++)
        result += weightedVector[i];
    return result;
}

/*
 * Computes the Diagonal Degree Matrix of the given data vectors
 * args: rows, cols - dimension on the data matrix
 *       vectors - a 2D matrix of the data vectors as rows.
 * return: diagonal degree matrix of the data vectors.
 */
double** ddg(double** weight_Adj_Mat, int rows) {
    int i;
    double** ddg_mat = Array_2D(rows, rows, 0);
    if (weight_Adj_Mat == NULL || ddg_mat == NULL){
        free_Mem(ddg_mat, rows);
        return NULL;
    }
    for(i = 0; i < rows; i++) {
        ddg_mat[i][i] = sumRow(weight_Adj_Mat[i], rows);
    }
    return ddg_mat;
}

/* Util function to multiply two N*N matrix.
 * NOTE: matrix B will be overwritten with the result.
 */
static double** mul_matrix(double** A, double** B, int N){
    int i, j, k;
    double** res = Array_2D(N, N, 0);
    if (res == NULL || A == NULL || B == NULL)
        return NULL;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
                res[i][j] += A[i][k] * B[k][j];

    for (i = 0; i < N; i ++)
        for (j = 0; j < N; j++)
            B[i][j] = res[i][j];
    free_Mem(res, N);
    return  B;
}


double** L_norm(double** vectors, int rows, int cols) {
    int i, j;
    double** L_norm_Mat;
    double** W = wam(vectors, rows, cols);
    double** D = ddg(W, rows);
    if (W == NULL || D == NULL)
        return NULL;

    for (i = 0; i < rows; i++)
        D[i][i] = pow(D[i][i], -0.5);

    L_norm_Mat = mul_matrix(mul_matrix(D, W,rows), D, rows);
    if (L_norm_Mat == NULL)
        return NULL;

    for (i = 0; i < rows; i++)
        for (j = 0; j < rows; j++) {
            if(L_norm_Mat[i][j] != 0)
                L_norm_Mat[i][j] = -1 * L_norm_Mat[i][j];
            if (i == j)
                L_norm_Mat[i][j]++;
        }

    free_Mem(W, rows);
    return L_norm_Mat;
}



static double** identity_mat(int n) {
    int i;
    double** I = Array_2D(n,n, 0);
    for (i = 0; i< n; i++) { I[i][i] = 1; }
    return I;
}

static double computeConvergence(double** matrix, double** newMatrix, int n){
    int i, j;
    double offA = 0.0;
    double offATag = 0.0;

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(i == j)
                continue;
            offA += pow(matrix[i][j], 2.0);
            offATag += pow(newMatrix[i][j], 2.0);
        }
    }
    return offA - offATag;
}

static int isDiaganol(double** matrix, int n) {
    int i, j;
    for(i = 0; i < n; i++) {
        for (j = 0;j < n; j++) {
            if (i == j)
                continue;
            if (matrix[i][j] != 0)
                return 0;
        }
    }
    return 1;
}

static void findLargestOffDiaganolElement(double** matrix, int n) {
    int i, j;
    double maxValue = (-1.0) * DBL_MAX;
    index_i = -1, index_j = -1;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if(i == j)
                continue;
            if(fabs(matrix[i][j]) > maxValue) {
                maxValue = fabs(matrix[i][j]);
                index_i = i, index_j = j;
            }
        }
    }
}

static void computeCAndS(double** matrix) {
    int i = index_i, j = index_j;
    double signTeta, t;
    double teta = matrix[j][j] - matrix[i][i];

    teta = teta / (2.0 * matrix[i][j]);
    if(teta >= 0)
        signTeta = 1.0;
    else
        signTeta = -1.0;

    t = signTeta / (fabs(teta) + sqrt(pow(teta, 2.0) + 1.0));
    C = 1.0 / sqrt(pow(t, 2.0) + 1.0);
    S = C * t;
}

static double** buildRotationMatrix(int n) {
    double** rotationMatrix = identity_mat(n);
    int i = index_i, j = index_j;

    rotationMatrix[i][i] = C;
    rotationMatrix[j][j] = C;
    rotationMatrix[i][j] = S;
    rotationMatrix[j][i] = S * (-1.0);
    return rotationMatrix;
}

static double** computeNewMatrix(double** matrix, int n) {
    int r, j;
    double** newMatrix = Array_2D(n,n,0);
    int I = index_i, J = index_j;
    for(r = 0; r < n; r++)
        for(j = 0; j < n; j++)
            newMatrix[r][j] = matrix[r][j];

    for(r = 0; r < n; r++) {
        if (r == I || r == J) {continue;}
        newMatrix[r][I] = (C * matrix[r][I]) - (S * matrix[r][J]);
        newMatrix[I][r] = newMatrix[r][I];
        newMatrix[r][J] = (C * matrix[r][J]) + (S * matrix[r][I]);
        newMatrix[J][r] = newMatrix[r][J];
    }
    newMatrix[I][I] = (pow(C, 2.0) * matrix[I][I]) + (pow(S, 2.0) * matrix[J][J]) - (2 * C * S * matrix[I][J]);
    newMatrix[J][J] = (pow(S, 2.0) * matrix[I][I]) + (pow(C, 2.0) * matrix[J][J]) + (2 * C * S * matrix[I][J]);
    newMatrix[I][J] = ((pow(C, 2.0) - pow(S, 2.0)) * matrix[I][J]) + (C * S) * (matrix[I][I] - matrix[J][J]);
    newMatrix[J][I] = newMatrix[I][J];
    return newMatrix;
}

static double** jacobi(double** matrix, int n) {
    const double epsilon = 0.00001;
    int i, rotationCount = 0;
    double convergence = epsilon + 1.0;
    double **newMatrix, **temp, **rotationMatrix = identity_mat(n);
    while(convergence > epsilon && rotationCount < 100 && isDiaganol(matrix, n) == 0) {
        findLargestOffDiaganolElement(matrix, n);
        computeCAndS(matrix);
        temp = rotationMatrix;
        rotationMatrix = mul_matrix(rotationMatrix, buildRotationMatrix(n), n);
        free_Mem(temp, n);
        newMatrix = computeNewMatrix(matrix, n);
        convergence = computeConvergence(matrix, newMatrix, n);
        free_Mem(matrix, n);
        matrix = newMatrix;
        rotationCount++;
    }
    eigen_values = (PAIR *) malloc(n * sizeof(PAIR));
    if (eigen_values == NULL)
        return NULL;
    for(i = 0; i < n; i++){
        eigen_values[i].val = matrix[i][i];
        eigen_values[i].index = i;
    }
    free_Mem(matrix, n);
    return rotationMatrix;
}


static void compute_T(double** U, double** T, int rows, int cols) {
    int i, j;
    double sum;
    for (i = 0; i < rows; i++) {
        sum = 0;
        for (j = 0; j < cols; j++)
            sum += pow(U[i][j], 2.0);

        sum = sqrt(sum);
        for (j = 0; j < cols; j++)
            if (sum != 0)
                T[i][j] = U[i][j] / sum;
    }
}

static void compute_U(double** matrix, double** U, int rows, int cols){
    int i, j, q;
     for (j = 0; j < cols; j++){
        q = eigen_values[j].index;
        for (i = 0; i < rows; i++){
            U[i][j] = matrix[i][q];
        }
    }
}

int compare (const void * a, const void * b)
{
    const PAIR *p1 = (PAIR *)a;
    const PAIR *p2 = (PAIR *)b;

    if (p1->val > p2->val) return 1;
    else if (p1->val < p2->val) return -1;
    else return 0;
}

static int eigen_gap(int n) {
    int i, argmax = 0;
    double max = 0;

    qsort(eigen_values, n, sizeof(PAIR), compare);
    n = n / 2;
    for (i = 0; i < n; i++){
        if (fabs(eigen_values[i].val - eigen_values[i+1].val) > max){
            max = fabs(eigen_values[i].val - eigen_values[i+1].val);
            argmax = i;
        }
    }
    return argmax+1;
}

double** spk(double** vectors, int rows, int cols, int* k) {
    double **T, **U;
    double** lnorm_mat = L_norm(vectors, rows, cols);
    double** mat = jacobi(lnorm_mat, rows);

    if (mat == NULL)
        return NULL;
    *k = eigen_gap(rows);
    if (*k == -1)
        return NULL;

    U = Array_2D(rows, *k, 0);
    T = Array_2D(rows, *k, 0);
    if (T == NULL || U == NULL)
        return NULL;

    compute_U(mat, U, rows, *k);
    compute_T(U, T, rows, *k);

    free_Mem(U, rows);
    free_Mem(mat, rows);
    return T;
}

#endif
