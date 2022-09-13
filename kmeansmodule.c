#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"


Matrix Spkmeans(char *file_name, int n, int d){
    return Eigengap(file_name, n, d);
}


Matrix MatOutput(int goal, char *file_name, int n, int d){
    switch (goal)
    { 
    case 0: return Eigengap(file_name, n, d);
        break;
    case 1: return Wam(file_name, n, d);
        break;
    case 2:  return Ddg(file_name, n, d);
        break;
    case 3:  return Lnorm(file_name, n, d);
        break;
    case 4:
            
        /* Checking matrix's size */
        if (n != d)
        {
            printf("Invalid Input!\n");
            exit(1);
        }

        observations = ReadObservations(file_name, n, n);
        lnorm = allocate_matrix(n,n);

        for (i=0; i<graph.N; i++) {
            graph.nodes[i]= observations[i];
        }   
        ind = 0;
        for (i = 0; i < lnorm.N; i++){
            for (j = 0; j < lnorm.N; j++){
                lnorm.edges[i][j] = observations[ind++];
            }
        }
        return Jacobi(lnorm);
        break;
    }
    exit(1);
}


double* ReadObservations(char* file_name, int n, int d){
    observations = (double *)malloc(n*d*sizeof(double));
    if (observations == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    ifp = fopen(file_name, "r");
    if (ifp == NULL)
    {
        freeObservations();
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n*d; i++){
        fscanf(ifp,"%lf", &observations[i]);
        fgetc(ifp); /* skip ',' */
    }
    fclose(ifp);

    return observations;
}



Matrix Wam(char* file_name, int n, int d) {
    observations = ReadObservations(file_name,n,d);
    graph = allocate_matrix(n,n);
    for (i = 0; i < graph.N; i++) {
        for (j = 0; j < graph.N; j++) {
            if (i == j) { /* self edge or no edge at all */
                 graph.edges[i][j] = 0;
            } else {  /* edge exist */
                double norm = 0.0;
                for (l = 0; l < d; l++){
                    norm += pow(observations[d*i+l] - observations[d*j+l], 2);
                }
                norm = sqrt(norm);
                
                graph.edges[i][j] = exp(-norm / 2);
            }
        }
    }
    return graph;
}


Matrix Ddg(char* file_name, int n, int d) {
    wam = Wam(file_name, n, d);
    for (i = 0; i < wam.N; i++) { /* wam[i][i] */
        for (j = 0; j < wam.N; j++) { /* sums and resets rows */
            if (i != j) {
                wam.edges[i][i] += wam.edges[i][j];
                wam.edges[i][j] = 0;
            }
        }
    }
    return wam;
}


Matrix Lnorm(char* file_name, int n, int d) {
    ddg = Ddg(file_name,n,d);
    freeObservations();
    wam = Wam(file_name,n,d);
    for (i = 0; i < wam.N; i++) {
        for (j = 0; j < wam.N; j++) {
            wam.edges[i][j] = - sqrt(1/ddg.edges[i][i]) * wam.edges[i][j] * sqrt(1/ddg.edges[j][j]);;
            if (i==j)
                wam.edges[i][j] = 1 + wam.edges[i][j];
        }
    }
    free_matrix(ddg);
    return wam;
}


Matrix allocate_matrix(int n, int m) {
    result.N = n;
    result.M = m;
    result.nodes = (double *)malloc(n*sizeof(double));
    if (result.nodes == NULL)
    {
        printf("An Error Has Occurred\n");
        freeObservations();
        exit(1);
    }
    result.edges = (double **)malloc(n*sizeof(double*));
    if (result.edges == NULL)
    {
        free_matrix(result);
        printf("An Error Has Occurred\n");
        freeObservations();
        exit(1);
    }
    /* allocates memory for edges */
    for (i = 0; i < n; i++){
        result.edges[i] = (double *)malloc(m*sizeof(double));
        if (result.edges[i] == NULL)
        {
            free_matrix(result);
            printf("An Error Has Occurred\n");
            freeObservations();
            exit(1);
        }
    }
    return result;
}



void free_matrix(Matrix A) {
    
    if (A.nodes != NULL)
        free(A.nodes);

    if (A.edges != NULL) {
        for (i = 0; i<A.N; i++) {
            if (A.edges[i] != NULL)
                free(A.edges[i]);
        }
        free(A.edges);

    }
}


/* Fills matrix with zeros */
void fillWithZeros(Matrix *A) {
    for (i = 0; i < A->N; i++)
        for (j = 0; j < A->N; j++)
            A->edges[i][j] = 0;
}

/* copy B to A */
void copyMatrix(Matrix *A, Matrix *B) {
    for (i = 0; i < A->N; i++)
        for (j = 0; j < A->M; j++)
            A->edges[i][j] = B->edges[i][j];
}


/* Multiply 2 matrices NxN in order AXB, put value in A*/
void multiply(Matrix *A, Matrix *B) {
    result = allocate_matrix(A->N, A->N);
   
    fillWithZeros(&result);

    for (i = 0; i < A->N; i++)
        for (j = 0; j < B->N; j++)
            for (k = 0; k < A->N; k++)
                result.edges[i][j] += A->edges[i][k] * B->edges[k][j];

    copyMatrix(A,&result);
    free_matrix(result);
}


Matrix transpose(Matrix A) {
    result = allocate_matrix(A.M, A.N);
    for (i=0; i<A.N;i++) {
        result.nodes[i] = A.nodes[i];
    }
    for (i = 0; i < A.M; i++)
        for (j = 0; j < A.N; j++) 
            result.edges[i][j] = A.edges[j][i];
    free_matrix(A);
    return result;
}


Matrix Jacobi(Matrix lnorm) {

    maxRotations = 100;
    e = pow(10, -5);
    rotation = 0;

    /* initial offA_tag */
    offA=0;
    offA_tag = 0;
        for (i = 0; i < lnorm.N; i++) {
            for (j = 0; j < lnorm.N; j++) {
                if(i!=j) {
                    offA_tag += pow(lnorm.edges[i][j], 2);
                }
            }
        }
    eigenVectors = allocate_matrix(lnorm.N, lnorm.N);
    fillWithZeros(&eigenVectors);
    for (i = 0; i<eigenVectors.N; i++) {
        eigenVectors.edges[i][i] = 1;
    }

    p = allocate_matrix(lnorm.N, lnorm.N);
    A_t = allocate_matrix(lnorm.N, lnorm.N);

    do {
        /* >>> find max abs value */
        /* coordinates for max value */
        m_i = 0;
        m_j = 0;
        maxValue = 0;

        for (i = 0; i < lnorm.N; i++) {
            for (j = i + 1; j < lnorm.N; j++) {
                if (fabs(lnorm.edges[i][j]) > maxValue) {
                    maxValue = fabs(lnorm.edges[i][j]);
                    m_i = i;
                    m_j = j;
                }
            }
        }
        /* <<< */

        /* >>> c, t */
        theta = (lnorm.edges[m_j][m_j] - lnorm.edges[m_i][m_i]) / (2.0*lnorm.edges[m_i][m_j]);
        abs_theta = theta;
        if (theta<0)
            abs_theta= - theta;
        t = 1.0 / (abs_theta + sqrt(theta*theta + 1));
        if (theta<0)
            t = -t;
        c = 1.0 / sqrt(t*t + 1);
        s = t * c;
        /* <<< */

        /* >>> build A_tag */
        fillWithZeros(&A_t);
        A = lnorm.edges;
        i = m_i, j = m_j;

        A_t.edges[i][i] = c*c*A[i][i] - 2*s*c*A[i][j] + s*s*A[j][j];
        A_t.edges[j][j] = s*s*A[i][i] + 2*s*c*A[i][j] + c*c*A[j][j];
        A_t.edges[i][j] = (c*c-s*s)*A[i][j] + s*c*(A[i][i]-A[j][j]);
        A_t.edges[j][i] = (c*c-s*s)*A[i][j] + s*c*(A[i][i]-A[j][j]);

        for (k = 0; k < A_t.N; k++ ){
            for (l = 0; l < A_t.N; l++) {
                if (k != i && k != j) {
                    if (l==i) {
                        A_t.edges[i][k]=c*A[i][k] - s*A[j][k];
                        A_t.edges[k][i]=c*A[i][k] - s*A[j][k];
                    } else if (l==j) {
                        A_t.edges[j][k]=s*A[i][k] + c*A[j][k];
                        A_t.edges[k][j]=s*A[i][k] + c*A[j][k];
                    } else {
                        A_t.edges[k][l] = A[k][l];
                    }
                }
                
            }
        }
        
        copyMatrix(&lnorm, &A_t);


        /* >>> calculate eigencectors */
        fillWithZeros(&p);

        /* Diagonal = 1 */
        for (i = 0; i<p.N; i++)
            p.edges[i][i]=1;

        p.edges[m_i][m_j] = s;
        p.edges[m_j][m_i] = -s;
        p.edges[m_i][m_i] = c;
        p.edges[m_j][m_j] = c;

        multiply(&eigenVectors,&p);

        /* <<< */

        /* >>> next */
        rotation++;
        offA = offA_tag;
        offA_tag = 0;
        for (i = 0; i < lnorm.N; i++) {
            for (j = 0; j<lnorm.N; j++) {
                if(i!=j) {
                    offA_tag += pow(lnorm.edges[i][j], 2);
                }
            }
        }
    }  while (rotation<maxRotations &&  offA - offA_tag > e);

    
    free_matrix(p); 
    free_matrix(A_t);  

    /* assign eigenValues */
    for (i = 0; i < lnorm.N; i++){
            eigenVectors.nodes[i] = lnorm.edges[i][i];
    }

    free_matrix(lnorm);
    return eigenVectors;
}

Matrix Eigengap(char* file_name, int n, int d) {
    lnorm = Lnorm(file_name, n, d);
    jacobi = Jacobi(lnorm);
    argmax = 0;
    max = 0.0;


    jacobi = transpose(jacobi);

    /* sort array */
    for (i = 0; i < n; ++i){
        for (j = i + 1; j < n; ++j){
            if (jacobi.nodes[i] > jacobi.nodes[j]){
                tmp = jacobi.nodes[i];
                tmpList = jacobi.edges[i];
                jacobi.nodes[i] = jacobi.nodes[j];
                jacobi.edges[i] = jacobi.edges[j];
                jacobi.nodes[j] = tmp;
                jacobi.edges[j] = tmpList;
            }
        }
    }
    jacobi = transpose(jacobi);

    /* find eigengap */
    for (i = 0; i < jacobi.N / 2; i++) {
        tmp = fabs(jacobi.nodes[i]-jacobi.nodes[i+1]);
        if (tmp > max){
            max = tmp;
            argmax=i;
        }
        
    }
    argmax++;

    if (argmax==1)
    {
        printf("An Error Has Occurred\n");
        free_matrix(jacobi);
        exit(1);
    }

    U = allocate_matrix(jacobi.N, argmax);
    for (i = 0; i < jacobi.N; i++)
    {
        for (j = 0; j < argmax; j++)
        {
            U.edges[i][j] = jacobi.edges[i][j];
        }
        
    }

    T = allocate_matrix(jacobi.N, argmax);
    for (i = 0; i < jacobi.N; i++)
    {
        for (j = 0; j < argmax; j++)
        {
            tmp=0;
            for (l = 0; l < argmax; l++)
            {
                tmp+=pow(U.edges[i][l],2);
            }
            if(tmp==0) {
                T.edges[i][j]=0;
            } else {
                T.edges[i][j] = U.edges[i][j] / pow(tmp,0.5);
            }

        }
    }    

    free_matrix(U);
    free_matrix(jacobi);

    return T;
}

double Norm(double number){
    double norm = number;
    return norm;
}


void freeObservations() {
    if (observations!=NULL)
        free(observations);
}


int main(int argc, char *argv[])
{
    FILE *ifp;
    if(argc == 3)
    {
        goal=argv[1];
        filename = argv[2];
    } 
    else if (argc == 4)
    {
        goal=argv[2];
        filename = argv[3];
    }
    else
    {
        printf("Invalid Input!\n");
        exit(1);
    }


   /* ----- Counting d ----- */
    ifp = fopen(filename, "r");
    if (ifp == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    d = 1;
    /* counting d */
    while ((cc = fgetc(ifp)) != '\n')
    {
        if (cc == ',')
            ++d;

    }
    fclose(ifp);
    /* counting n */
    ifp = fopen(filename, "r");
    if (ifp == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    n = 0;
    while ((cc = fgetc(ifp)) != EOF)
    {
        if (cc == '\n') 
            ++n;
    }
    fclose(ifp);

    goal_int = 0;

    if (strcmp(goal, "wam") ==0 )
        goal_int=1;
    else if (strcmp(goal, "ddg") == 0)
        goal_int=2;
    else if (strcmp(goal, "lnorm") == 0)
        goal_int=3;
    else if (strcmp(goal, "jacobi") == 0)
        goal_int=4;
    else
        {
        printf("Invalid Input!\n");
        return 1;
        }
    matrix = MatOutput(goal_int, filename, n, d);

    /* Prints the desire matrix */
    if(goal_int==4) {
        for (i = 0; i < matrix.N; i++){
            if(fabs(matrix.nodes[i])<0.00005){
                printf("0.0000");
            } else {
                printf("%.4f", matrix.nodes[i]);
            }
            if (i + 1 < matrix.N)
                printf("%c", ',');
         }
        printf("\n");

    }
    for (i = 0; i < matrix.N; i++){
        for (j = 0; j < matrix.N; j++){
            if(fabs(matrix.edges[i][j])<0.00005){
                printf("0.0000");
            } else {
                printf("%.4f", matrix.edges[i][j]);
            }
            if (j + 1 < matrix.N)
                printf("%c", ',');
        }
        printf("\n");
    }
    freeObservations();
    free_matrix(matrix);

    return 0;
}
