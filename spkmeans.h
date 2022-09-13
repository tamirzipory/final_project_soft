#ifndef SPKMEANS_H_
#define SPKMEANS_H_


typedef struct Matrix {
    int N;
    int M;
    double * nodes;
    double ** edges;
} Matrix;



/* declareing variables */
int i, j, k, l, ind, d, n, maxRotations, rotation, m_i, m_j, argmax, goal_int, cc;
Matrix graph, lnorm, eigenVectors, p, A_t, matrix, U, T, wam, ddg, result, lnorm, jacobi;
FILE * ifp;
double offA, offA_tag, maxValue, s, **A, e, theta, abs_theta, t, c, max, tmp;
double * tmpList;
char *goal, *filename, *file_name;

/* declareing functions */
Matrix MatOutput(int goal, char *file_name, int n, int d);
double * ReadObservations(char* file_name, int n, int d);
Matrix Wam(char* file_name, int n, int d);
Matrix Ddg(char* file_name, int n, int d);
Matrix Lnorm(char* file_name, int n, int d);
Matrix allocate_matrix(int n, int m);
void free_matrix(Matrix A);
void fillWithZeros(Matrix *A);
void copyMatrix(Matrix *A, Matrix *B);
void multiply(Matrix *A, Matrix *B);
Matrix transpose(Matrix A);
Matrix Jacobi(Matrix lnorm);
Matrix Eigengap(char* file_name, int n, int d);
double Norm(double number);
void freeObservations(void);

double * observations;

#endif
