#define PY_SSIZE_T_CLEAN

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "spkmeans.h"
#define eps pow(10, -5)

/* GOLOBAL ELEMENTS: */

typedef struct eigenVector
{
    double eigenValue;
    int index;
} eigenVector;

double **clusters;
double **centroids;
double **vector_list;
int *clustersindexes;
int vector_len;
int vector_num;
int k;
int clusters_num;
double **weightedAdjMatrix;
double **diagDegMatrix;
double **LnormMatrix;
eigenVector *eigenVectors;
double *eigenValues;
double **VMat;
double **UMat;
double **TMat;
int max_iter;
float inputK;
char *goal;
double **python_centroids;

/* -------------------- FOR SPK ----------------- */

int eigenComperator(const void *a, const void *b)
{
    struct eigenVector *A = (struct eigenVector *)a;
    struct eigenVector *B = (struct eigenVector *)b;
    A = (eigenVector *)a;
    B = (eigenVector *)b;
    if (A->eigenValue == B->eigenValue)
    {
        return A->index - B->index;
    }
    else
    {
        if (A->eigenValue > B->eigenValue)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
}

void sortEigenVectors()
{
    int i;
    eigenVectors = (eigenVector *)calloc(vector_num, vector_num * sizeof(eigenVector));
    assertmsg(eigenVectors != NULL);
    for (i = 0; i < vector_num; i++)
    {
        eigenVectors[i].index = i;
        eigenVectors[i].eigenValue = eigenValues[i];
    }

    qsort(eigenVectors, vector_num, sizeof(eigenVector), eigenComperator);
    for (i = 0; i < vector_num; i++)
    {
        eigenValues[i] = eigenVectors[i].eigenValue;
    }
}

int eigengapHeuristic()
{
    int i, maxIndex, k = 0;
    double maxGap = -1.0;
    double *eigenGaps;
    calcJacobi(Lnorm());
    sortEigenVectors();
    eigenGaps = (double *)calloc(vector_num - 1, sizeof(double));
    for (i = 0; i < vector_num - 1; i++)
    {
        eigenGaps[i] = fabs(eigenValues[i] - eigenValues[i + 1]);
    }
    maxIndex = (int)floor(vector_num / 2);
    for (i = 0; i < maxIndex; i++)
    {
        if (eigenGaps[i] > maxGap)
        {
            maxGap = eigenGaps[i];
            k = i;
        }
    }
    free(eigenGaps);
    return k + 1;
}

void createTMat()
{ /*using Vmat*/
    int i, j;
    double sum;
    /* Transpose VMat */
    transpMat(VMat);
    UMat = (double **)calloc(vector_num, k * sizeof(double));
    assertmsg(UMat != NULL);
    for (i = 0; i < vector_num; i++)
    {
        UMat[i] = (double *)calloc(k, sizeof(double));
        assertmsg(UMat[i] != NULL);
        for (j = 0; j < k; j++)
        {
            UMat[i][j] = VMat[eigenVectors[j].index][i];
        }
    }
    /*normalizing UMat*/
    for (i = 0; i < vector_num; i++)
    {
        sum = 0;
        for (j = 0; j < k; j++)
        {
            sum += (UMat[i][j]) * (UMat[i][j]);
        }
        sum = sqrt(sum);
        if (sum != 0)
        {
            for (j = 0; j < k; j++)
            {
                UMat[i][j] = (UMat[i][j]) / sum;
            }
        }
    }
    TMat = UMat; /* TMat is a normalized Umat*/
}

void init_centroids(int k)
{
    int i, j, x;
    assertmsg(k < vector_num);
    centroids = (double **)calloc(k, vector_len * sizeof(double));
    assertmsg(centroids != NULL);
    for (i = 0; i < k; i++)
    {
        centroids[i] = (double *)calloc(vector_len, sizeof(double));
        assertmsg(centroids[i] != NULL);
    }
    for (x = 0; x < k; x++)
    {
        for (j = 0; j < vector_len; j++)
        {
            centroids[x][j] = vector_list[x][j];
        }
    }
    clusters = (double **)calloc(k, sizeof(double *));
}

double distance(double *v1, double *v2)
{
    double res;
    int length;
    int i;
    res = 0;
    length = vector_len;
    for (i = 0; i < length; i++)
    {
        res += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return res;
}

int min_dist_centroid(double *v)
{
    double minimum;
    int ind;
    double dist;
    int i;
    minimum = distance(v, centroids[0]);
    ind = 0;
    for (i = 0; i < k; i++)
    {
        dist = distance(v, centroids[i]);
        if (dist < minimum)
        {
            minimum = dist;
            ind = i;
        }
    }
    return ind;
}

void vector_to_cluster(int k)
{
    int *clusterssizes; /*for realloc*/
    int ind;
    int i;
    free(clustersindexes);
    clustersindexes = (int *)calloc(k, sizeof(int));
    assertmsg(clustersindexes != NULL);
    clusterssizes = (int *)calloc(k, sizeof(int));
    assertmsg(clusterssizes != NULL);
    for (i = 0; i < k; i++)
    { /*initialize each cluster's size to 100*/
        clusterssizes[i] = 100;
    }
    for (i = 0; i < k; i++)
    {
        free(clusters[i]);
        clusters[i] = (double *)calloc(100, sizeof(double));
        assertmsg(clusters[i] != NULL);
    }
    for (i = 0; i < vector_num; i++)
    {
        ind = min_dist_centroid(vector_list[i]);
        if (clustersindexes[ind] > ((clusterssizes[ind]) / 2))
        { /*Increase if necessary*/
            clusters[ind] = (double *)realloc(clusters[ind], 2 * clusterssizes[ind] * sizeof(double));
            clusterssizes[ind] *= 2;
        }
        clusters[ind][clustersindexes[ind]] = i;
        clustersindexes[ind]++; /*increase number of vectors in specified cluster*/
    }
    free(clusterssizes);
}

double *cluster_to_centroid(int index)
{
    int i;
    int j;
    int vector_index;
    int num = clustersindexes[index]; /* number of vectors in given cluster */
    double *res = (double *)calloc(vector_len, sizeof(double));
    assertmsg(res != NULL);
    if (num != 0)
    {
        for (i = 0; i < vector_len; i++)
        {
            for (j = 0; j < num; j++)
            {
                vector_index = (int)clusters[index][j]; /* not actual vector but index in vector_list */
                res[i] += vector_list[vector_index][i]; /*relevant cluster*/
            }
        }

        for (i = 0; i < vector_len; i++)
        {
            res[i] = res[i] / (num);
        }
    }
    else
    {
        for (i = 0; i < vector_len; i++)
        {
            res[i] = centroids[index][i];
        }
    }
    return res;
}

int areequal(double *arr1, double *arr2)
{
    int length;
    int y;
    length = vector_len;
    for (y = 0; y < length; y++)
    {
        if (arr1[y] != arr2[y])
        {
            return 0;
        }
    }
    return 1;
}

int update_centroids()
{
    int changed = 0;
    int x, i, j;
    for (i = 0; i < k; i++)
    {
        double *newcentroid;
        double *res;
        newcentroid = (double *)calloc(vector_len, sizeof(double));
        assertmsg(newcentroid != NULL);
        res = cluster_to_centroid(i);
        for (j = 0; j < vector_len; j++)
        {
            newcentroid[j] = res[j];
        }
        if (areequal(centroids[i], newcentroid) == 0)
        {
            changed++;
        }
        for (x = 0; x < vector_len; x++)
        {
            centroids[i][x] = newcentroid[x];
        }
        free(newcentroid);
        free(res);
    }
    return (changed != 0);
}

double **calccentroids(int max_iter)
{
    int counter;
    int isequal;
    assertmsg(clusters_num == 0 || clusters_num > 0);
    counter = 0;
    isequal = 1;
    clusters = (double **)calloc(clusters_num, sizeof(double *)); /* originally in init_centroids */
    while (counter < max_iter && isequal == 1)
    {
        vector_to_cluster(clusters_num);
        isequal = update_centroids();
        counter++;
    }
    freearray(clusters, clusters_num);
    return centroids;
}
void fullSpectral()
{
    int counter;
    int isequal;
    int max_iter;
    int newK;
    max_iter = 300;
    newK = eigengapHeuristic();
    if (k == 0)
    {
        k = newK;
    }
    vector_len = k;
    createTMat();
    freearray(vector_list, vector_num);
    vector_list = UMat;

    init_centroids(k);

    counter = 0;
    isequal = 1;

    while (counter < max_iter && isequal == 1)
    {
        vector_to_cluster(k);
        isequal = update_centroids();
        counter++;
    }
}

/* ------------------------- END SPK FUNC ------------------------ */
/* ASSERT MES*/
void assertmsg(int x)
{
    if (!x)
    {
        printf("An Error Has Occured\n");
        assert(x);
    }
}

/* NUMBER OF LINES IN FILE - THE N*/
int find_vector_num(FILE *file)
{
    int lines;
    char ch;
    char prev_char;
    lines = 0;
    ch = fgetc(file);
    prev_char = ch;
    while (ch != EOF)
    {
        if (ch == '\n')
        {
            lines++;
        }
        prev_char = ch;
        ch = fgetc(file);
    }
    if (prev_char != '\n')
    {
        lines++;
    }
    vector_num = lines;
    return lines;
}

/* NUMBER OF LINES IN FILE - THE d*/
int find_vector_len(FILE *file)
{

    int vectorcomp;
    char ch;
    vectorcomp = 1;
    ch = fgetc(file);
    while (ch != '\n')
    {
        if (ch == ',')
        {
            vectorcomp++;
        }
        ch = fgetc(file);
    }
    vector_len = vectorcomp;
    return vectorcomp;
}

void readfile(FILE *file)
{
    char *split_line;
    int i, j;
    char line[1000 * sizeof(char)];
    find_vector_len(file); /* initail vector_len */
    rewind(file);
    find_vector_num(file); /* initail vector_num */
    rewind(file);
    vector_list = (double **)malloc(vector_num * sizeof(double)); /* vector_num (n) tabs that will contain double*/
    assertmsg(vector_list != NULL);
    for (i = 0; i < vector_num; i++)
    {
        vector_list[i] = (double *)calloc(vector_len, sizeof(double));
        assertmsg(vector_list[i] != NULL);
    }
    j = 0;
    while (fgets(line, 1000, file) != NULL)
    {                                   /* reads line by line as a string*/
        split_line = strtok(line, ","); /* split line by commas*/
        for (i = 0; i < vector_len; i++)
        {
            vector_list[j][i] = atof(split_line); /* convert str to float*/
            split_line = strtok(NULL, ",");
        }
        j++;
    }
    fclose(file);
}

double **matMult(double **A, double **B)
{
    int i, j, k, l;
    double **res = (double **)calloc(vector_num, vector_num * sizeof(double));
    assertmsg(res != NULL);
    for (i = 0; i < vector_num; i++)
    {
        res[i] = (double *)calloc(vector_num, sizeof(double));
        assertmsg(res[i] != NULL);
    }
    for (j = 0; j < vector_num; j++)
    {
        for (k = 0; k < vector_num; k++)
        {
            res[j][k] = 0;
            for (l = 0; l < vector_num; l++)
            {
                res[j][k] += A[j][l] * B[l][k];
            }
        }
    }
    return res;
}

/* RETURN THE WAM MATRIX */
double **weightedAdjMat()
{
    int i, j, k;
    weightedAdjMatrix = (double **)calloc(vector_num, vector_num * sizeof(double));
    assertmsg(weightedAdjMatrix != NULL);
    for (i = 0; i < vector_num; i++)
    {
        weightedAdjMatrix[i] = (double *)calloc(vector_num, sizeof(double));
        assertmsg(weightedAdjMatrix[i] != NULL);
    }
    for (j = 0; j < vector_num; j++)
    {
        double *v1 = vector_list[j];
        for (k = j + 1; k < vector_num; k++)
        {
            double *v2 = vector_list[k];
            weightedAdjMatrix[j][k] = weightedAdjMat_calc(v1, v2);
            weightedAdjMatrix[k][j] = weightedAdjMatrix[j][k]; /*symmetry*/
        }
    }
    return weightedAdjMatrix;
}

double weightedAdjMat_calc(double *v1, double *v2)
{
    double res;
    double dist = 0;
    int i;
    for (i = 0; i < vector_len; i++)
    {
        dist += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    res = sqrt(dist);
    res = -((res) / 2);
    res = exp(res);
    return res;
}

double **Lnorm()
{
    int i, j, k;
    double **res;
    weightedAdjMat();
    diagDegMat();
    for (i = 0; i < vector_num; i++)
    { /* D^(-1/2) */
        diagDegMatrix[i][i] = 1 / (sqrt(diagDegMatrix[i][i]));
    }
    res = matMult(diagDegMatrix, weightedAdjMatrix);
    LnormMatrix = matMult(res, diagDegMatrix);
    for (j = 0; j < vector_num; j++)
    {
        for (k = 0; k < vector_num; k++)
        {
            if (j == k)
            { /*minus id mat*/
                LnormMatrix[j][k] = 1 - LnormMatrix[j][k];
            }
            else
            {
                LnormMatrix[j][k] = -(LnormMatrix[j][k]);
            }
        }
    }
    freearray(res, vector_num);
    return LnormMatrix;
}

/* PRINTING THE MAT BY GIVIEN ROW AND COL */
void printMat(double **Mat, int rowNum, int colNum)
{
    int i, j;
    for (i = 0; i < rowNum; i++)
    {
        for (j = 0; j < colNum; j++)
        {
            if ((Mat[i][j] < 0) && (Mat[i][j] > -0.00005))
            {
                Mat[i][j] = 0;
            }
            if (j == colNum - 1)
            {
                printf("%0.4f", Mat[i][j]);
                if (i < rowNum - 1)
                {
                    printf("\n");
                }
            }
            else
            {
                printf("%0.4f,", Mat[i][j]);
            }
        }
    }
}

double **diagDegMat()
{
    int i, j, k;
    double sum;
    diagDegMatrix = (double **)calloc(vector_num, vector_num * sizeof(double));
    assertmsg(diagDegMatrix != NULL);
    for (i = 0; i < vector_num; i++)
    {
        diagDegMatrix[i] = (double *)calloc(vector_num, sizeof(double));
        assertmsg(diagDegMatrix[i] != NULL);
    }
    for (j = 0; j < vector_num; j++)
    {
        sum = 0;
        for (k = 0; k < vector_num; k++)
        {
            sum += weightedAdjMatrix[j][k];
        }
        diagDegMatrix[j][j] = sum;
    }
    return diagDegMatrix;
}

void calcRotationMat(double **P, double c, double s, int row, int column)
{
    int i, j;
    for (i = 0; i < vector_num; i++)
    {
        for (j = 0; j < vector_num; j++)
        { /*id mat*/
            if (i == j)
            {
                P[i][j] = 1;
            }
            else
            {
                P[i][j] = 0;
            }
        }
    }
    P[row][column] = s;
    P[column][row] = -s;
    P[row][row] = c;
    P[column][column] = c;
}

double calc_theta(double **Mat, int i, int j)
{
    double res, up, down;
    up = Mat[j][j] - Mat[i][i];
    down = 2 * Mat[i][j];
    res = up / down;
    return res;
}

double calc_t(double theta)
{
    int sign;
    double down;
    if (theta >= 0)
    {
        sign = 1;
    }
    else
    {
        sign = -1;
    }
    down = fabs(theta) + sqrt(theta * theta + 1);
    return sign / down;
}

double calc_c(double t)
{
    double res;
    double down = sqrt(t * t + 1);
    res = 1 / down;
    return res;
}

double calc_s(double t, double c) { return t * c; }

int isConverged(double **A, double **Aprime)
{
    double epsilon = pow(10, -5);
    int i, j, k, l;
    double res;
    double sumA = 0;
    double sumAprime = 0;
    for (i = 0; i < vector_num; i++)
    {
        for (j = 0; j < vector_num; j++)
        {
            if (i != j)
            {
                sumA += (A[i][j]) * (A[i][j]);
            }
        }
    }
    for (k = 0; k < vector_num; k++)
    {
        for (l = 0; l < vector_num; l++)
        {
            if (k != l)
            {
                sumAprime += (Aprime[k][l]) * (Aprime[k][l]);
            }
        }
    }
    res = sumA - sumAprime;
    if (res <= epsilon)
    {
        return 1;
    }
    return 0;
}

void calcAprime(double **A, double **Aprime, int i, int j, double c, double s)
{
    int k;
    for (k = 0; k < vector_num; k++)
    {
        if ((k != i) && (k != j))
        {
            Aprime[k][i] = c * A[k][i] - s * A[k][j];
            Aprime[i][k] = Aprime[k][i];
            Aprime[k][j] = c * A[k][j] + s * A[k][i];
            Aprime[j][k] = Aprime[k][j];
        }
    }
    Aprime[i][i] = c * c * A[i][i] + s * s * A[j][j] - 2 * s * c * A[i][j];
    Aprime[j][j] = s * s * A[i][i] + c * c * A[j][j] + 2 * s * c * A[i][j];
    Aprime[i][j] = 0;
    Aprime[j][i] = Aprime[i][j];
}

double **calcJacobi(double **A)
{
    double **Aprime;
    double **P;
    double **tmp;
    int i, j, k, l;
    int row = 0;
    int column = 1;
    double c, t, s, theta;
    int isConvergedBool = 0;
    int counter = 0;
    Aprime = (double **)calloc(vector_num, vector_num * sizeof(double));
    assertmsg(Aprime != NULL);
    for (i = 0; i < vector_num; i++)
    {
        Aprime[i] = (double *)calloc(vector_num, sizeof(double));
        assertmsg(Aprime[i] != NULL);
    }
    P = (double **)calloc(vector_num, vector_num * sizeof(double));
    assertmsg(P != NULL);
    for (j = 0; j < vector_num; j++)
    {
        P[j] = (double *)calloc(vector_num, sizeof(double));
        assertmsg(P[j] != NULL);
    }
    VMat = (double **)calloc(vector_num, vector_num * sizeof(double));
    assertmsg(VMat != NULL);
    for (k = 0; k < vector_num; k++)
    {
        VMat[k] = (double *)calloc(vector_num, sizeof(double));
        assertmsg(VMat[k] != NULL);
    }
    for (l = 0; l < vector_num; l++)
    { /* init to id mat*/
        VMat[l][l] = 1;
    }

    for (i = 0; i < vector_num; i++)
    { /*copy A to Aprime*/
        for (j = 0; j < vector_num; j++)
        {
            Aprime[i][j] = A[i][j];
        }
    }
    while ((isConvergedBool == 0) && (counter < 100))
    {
        for (i = 0; i < vector_num; i++)
        { /*finds max off-diagonal indices*/
            for (j = i + 1; j < vector_num; j++)
            {
                if (fabs(A[i][j]) > fabs(A[row][column]))
                {
                    row = i;
                    column = j;
                }
            }
        }

        if (A[row][column] == 0)
        {
            break;
        }

        theta = calc_theta(A, row, column);
        t = calc_t(theta);
        c = calc_c(t);
        s = calc_s(t, c);

        calcRotationMat(P, c, s, row, column);
        tmp = VMat;
        VMat = matMult(VMat, P);
        freearray(tmp, vector_num);
        calcAprime(A, Aprime, row, column, c, s);
        isConvergedBool = isConverged(A, Aprime);

        for (i = 0; i < vector_num; i++)
        { /*copy Aprime to A*/
            for (j = 0; j < vector_num; j++)
            {
                A[i][j] = Aprime[i][j];
            }
        }
        counter++;
    }
    eigenValues = (double *)calloc(vector_num, sizeof(double));
    assertmsg(eigenValues != NULL);
    for (i = 0; i < vector_num; i++)
    {
        eigenValues[i] = Aprime[i][i];
    }
    freearray(Aprime, vector_num);
    freearray(P, vector_num);
    return A;
}

double **transpMat(double **Mat)
{
    int i, j;
    double temp;
    for (i = 1; i < vector_num; i++)
    {
        for (j = 0; j < i; j++)
        {
            temp = VMat[i][j];
            Mat[i][j] = Mat[j][i];
            Mat[j][i] = temp;
        }
    }
    return Mat;
}

void freearray(double **array, int length)
{
    int i;
    for (i = 0; i < length; i++)
    {
        free(array[i]);
    }
    free(array);
}
void assert_double_mat(double **mat)
{
    /*Replaces assert*/
    if (mat == NULL)
    {
        printf("An Error Has Occured");
        exit(0);
    }
}

void assert_double_arr(const double *arr)
{
    /*Replaces assert*/
    if (arr == NULL)
    {
        printf("An Error Has Occured");
        exit(0);
    }
}

double **gen_id_mat(int N)
{
    /*Generates ID matrix*/
    int i, j;
    double **I;
    double *block;

    block = calloc(N * N, sizeof(double));
    assert_double_arr(block);
    I = calloc(N, sizeof(double *));
    assert_double_mat(I);
    for (i = 0; i < N; i++)
    {
        I[i] = block + i * N;
        for (j = 0; j < N; j++)
        {
            if (i == j)
            {
                I[i][j] = 1.0;
            }
        }
    }
    return I;
}
double off(double **A, int N)
{
    /*Calculates Off(A)^2*/
    int i, j;
    double res;

    res = 0;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i != j)
            {
                res += pow(A[i][j], 2);
            }
        }
    }
    return res;
}
int sign(double num)
{
    if (num < 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}
void V_multi_P(double **V, double s, double c, int N, int i, int j)
{
    /*V = VP; since P is almost-diagonal, no need to perform full matrix multiplication*/
    int r;
    double V_ri, V_rj;

    for (r = 0; r < N; r++)
    {
        V_ri = V[r][i];
        V_rj = V[r][j];

        V[r][i] = (c * V_ri) - (s * V_rj);
        V[r][j] = (s * V_ri) + (c * V_rj);
    }
}
double *get_ith_column(double **mat, int col_ind, int N)
{
    /*returns [mat[0][col_ind],mat[1][col_ind],...,mat[N-1][col_ind]]*/
    int i;
    double *col;
    col = calloc(N, sizeof(double));
    assert_double_arr(col);
    for (i = 0; i < N; i++)
    {
        col[i] = mat[i][col_ind];
    }
    return col;
}
void assert_int_arr(const int *arr)
{
    /*Replaces assert*/
    if (arr == NULL)
    {
        printf("An Error Has Occured");
        exit(0);
    }
}

int *max_indices_off_diag(double **A, int N)
{
    /*returns [i,j] so that A[i][j] is the off-diagonal largest absolute element in A*/
    double val;
    int i, j, max_i, max_j;
    int *arr;

    val = -1;
    max_i = 0;
    max_j = 0;

    arr = calloc(2, sizeof(int));
    assert_int_arr(arr);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i != j)
            {
                if (fabs(A[i][j]) > val)
                {
                    val = fabs(A[i][j]);
                    max_i = i;
                    max_j = j;
                }
            }
        }
    }
    arr[0] = max_i;
    arr[1] = max_j;
    return arr;
}

void A_to_A_tag(double **A, double **V, int N)
{
    /*Calculates A' from A using the relation between them as explained in 1.2.6*/
    int i, j, r;
    int *arr_max;
    double theta, s, t, c, a_ri, a_rj, a_ii, a_jj;

    arr_max = max_indices_off_diag(A, N); /*Finding pivot - 1.2.3*/
    i = arr_max[0];
    j = arr_max[1];
    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    t = sign(theta) / (fabs(theta) + sqrt((pow(theta, 2)) + 1));
    c = 1 / sqrt((pow(t, 2)) + 1);
    s = t * c;
    V_multi_P(V, s, c, N, i, j); /*V = VP*/

    for (r = 0; r < N; r++)
    {
        if ((r != j) && (r != i))
        {
            a_ri = c * A[r][i] - s * A[r][j];
            a_rj = c * A[r][j] + s * A[r][i];
            A[r][i] = a_ri;
            A[r][j] = a_rj;
            A[j][r] = a_rj;
            A[i][r] = a_ri;
        }
    }
    /*1.2.3*/
    a_ii = pow(c, 2) * A[i][i] + pow(s, 2) * A[j][j] - 2 * s * c * A[i][j];
    a_jj = pow(c, 2) * A[j][j] + pow(s, 2) * A[i][i] + 2 * s * c * A[i][j];
    A[j][j] = a_jj;
    A[i][i] = a_ii;
    A[i][j] = 0.0;
    A[j][i] = 0.0;

    free(arr_max);
}

double **jacobi(double **A, int N)
{
    /*Calculates and prints the eigenvalues and eigenvectors as described in 1.2.1*/
    int iter, max_iter;
    double diff;
    double **V;

    iter = 0;
    max_iter = 100;
    V = gen_id_mat(N);
    diff = 1;

    while (diff > eps && iter < max_iter)
    { /*Stops when off(A)^2- - off(A')^2 <= epsilon OR 100 iterations*/
        double off_A = off(A, N);
        iter++;
        A_to_A_tag(A, V, N); /*Also change V during A->A' iteration so that eventually V = P1*P2*P3*...*/
        diff = off_A - off(A, N);
    }
    return V;
}
void print_double(double num)
{
    /*Prevents "-0.0000" situation*/
    if ((num < 0) && (num > -0.00005))
    {
        printf("0.0000");
    }
    else
    {
        printf("%.4f", num);
    }
}

void print_row(double *row, int len)
{
    /*Prints array (row) according to given dimension: len == number of items*/
    int i;
    for (i = 0; i < len; i++)
    {
        print_double(row[i]);
        if (i != len - 1)
        {
            printf(",");
        }
    }
}
double *get_diag(double **mat, int N)
{
    /*returns [mat[0][0],mat[1][1],...,mat[N-1][N-1]]*/
    int i;
    double *diag;
    diag = calloc(N, sizeof(double));
    assert_double_arr(diag);

    for (i = 0; i < N; i++)
    {
        diag[i] = mat[i][i];
    }
    return diag;
}
void free_mat(double **mat)
{
    /*Since the matrices are allocated as a contigous block, we first free the block of N*k size and then the pointer*/
    free(mat[0]);
    free(mat);
}

/* ----------------------------------------------MAIN--------------------------------------------------- */
int main(int argc, char *argv[])
{
    FILE *file;
    if (!(argc == 3))
    {
        printf("Invalid Input!\n");
        assert(argc == 3);
    }
    goal = argv[1];

    file = fopen(argv[2], "r");
    readfile(file);

    /* JACOBI CASE
    if (strcmp(goal, "jacobi") == 0)
    {
        double **W, **V;
        double *eigenvalues;
        int i, N;

        W = vector_list;
        N = vector_num;
        V = jacobi(W, N);
        eigenvalues = get_diag(W, N);
        print_row(eigenvalues, N);
        printf("\n");
        for (i = 0; i < N; i++)
        {
            double *eigenvector = get_ith_column(V, i, N);
            print_row(eigenvector, N);
            if (i < N - 1)
            {
                printf("\n");
            }
            free(eigenvector);
        }
        printf("============================= SIGN 1 ================== \n");
        free_mat(W);
        printf("============================= SIGN 2 ================== \n");

        free_mat(V);
        printf("============================= SIGN 3 ================== \n");

        free(eigenvalues);
        printf("============================= SIGN 4 ================== \n");
    }*/
    if (strcmp(goal, "wam") == 0)
    {
        printMat(weightedAdjMat(), vector_num, vector_num);
        freearray(weightedAdjMatrix, vector_num);
        freearray(vector_list, vector_num);
    }
    /* DDG CASE */
    else if (strcmp(goal, "ddg") == 0)
    {
        weightedAdjMat();
        printMat(diagDegMat(), vector_num, vector_num);
        freearray(weightedAdjMatrix, vector_num);
        freearray(diagDegMatrix, vector_num);
        freearray(vector_list, vector_num);
    }
    /* LNORM CASE */
    else if (strcmp(goal, "lnorm") == 0)
    {
        printMat(Lnorm(), vector_num, vector_num);
        freearray(weightedAdjMatrix, vector_num);
        freearray(diagDegMatrix, vector_num);
        freearray(LnormMatrix, vector_num);
        freearray(vector_list, vector_num);
    }
    /* JACOBI CASE*/
    else if (strcmp(goal, "jacobi") == 0)
    {
        int i;
        double **A;
        printf("--------------------------END JACOBI- before 1 =========");
        A = calcJacobi(vector_list);
        printf("--------------------------END JACOBI- after 1 =========");
        for (i = 0; i < vector_num; i++)
        {
            if ((A[i][i] < 0) && (A[i][i] > -0.00005))
            {
                A[i][i] = 0;
            }
            if (i == vector_num - 1)
            {
                printf("%0.4f \n", A[i][i]);
            }
            else
            {
                printf("%0.4f,", A[i][i]);
            }
        }
        printf("--------------------------END JACOBI-1 =========");
        transpMat(VMat);
        printf("--------------------------END JACOBI -2 ======");
        printMat(VMat, vector_num, vector_num);
        printf("--------------------------END JACOBI -3 ========");
        freearray(VMat, vector_num);
        printf("--------------------------END JACOBI -4 -=======");
        free(eigenValues);
        printf("--------------------------END JACOBI --5 ========");
        /* segmantion fault ??? */
        freearray(vector_list, vector_num);
        printf("--------------------------END JACOBI --6 ========");
    }
    else
    {
        assertmsg(0 != 0);
    }
    return 0;
}
