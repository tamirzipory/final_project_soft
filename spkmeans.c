#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "spkmeans.h"

/* for kmeans */
int fitToKmeans(double *points, int n, int dim, int k, int maxIter, double *centroids, double *results);
double euclid(double x[], double y[], int dim);
int areCentroidsEqual(double *oldArr, double *newArr, int k, int dim);
void reArrangePointsToClusters(double *points, int *pointToWhichCluster, double *centroidsArr, int n, int k, int dim);
void reSelectCentroids(double *points, int *clustersCounters, int *pointToWhichCluster, double *centroidsArr, double *centroidsSum, int n, int k, int dim);
void copyResults(double *newArr, int k, int dim, double *results);
void freeArraysMemory(double *points, int *pointToWhichCluster, int *clustersCounters,
                      double *centroidsOdd, double *centroidsEven,
                      double *newCentroidAverageBeforeDivision);
int kmeans(int n, int k, int dim, int maxIter, double *points, double *centroids, double *results);

/* for SPkmeans */
int atoi(const char *nptr);
double atof(const char *nptr);
void rewind(FILE *stream);
char *strtok(char *str, const char *delim);
double pow(double x, double y);
int commasCount(char *p);
double vectorNorm(double *vector, int dim);
double euclideanDistance(double x[], double y[], int dim);
void subTwoMatrices(double *matrix1, double *matrix2, int rowDim1, int colDim1, double *subMatrix);
void multiTwoMatrices(double *matrix1, double *matrix2, int rowDim1, int colDim1, int colDim2, double *multiMatrix);
void transpose(double *matrix, int rowDim, int colDim, double *transposeMatrix);
void negativeSqrtMatrix(double *matrixD, int dim, double *negativeSqrtMatrix);
void WeightedAdjacencyMatrixFromPoints(double *points, int n, int dim, double *weightedAdjancy);
void diagonalMatrixFromW(double *weightedMatrix, int len, double *diagonal);
void normalizedGraphLaplacian(double *Dsqrt, double *W, int dim, double *Lnorm);
void renormalizingMatrixRow(double *matrixU, int rowDim, int colDim, double *renormalizedRowMatrix);
double sumSquareOffDiagonalElem(double *matrix, int rowDim, int colDim);
void IndexOfMaxValInMatrix(double *matrix, int rowDim, int colDim, int *indicesIandJ);
void computingCandSforRotationMatrix(double *A, int rowDim, int *indices, double *valsCandS);
void eigenvaluesSortingAndHeuristic(double *Atag, double *V, double **resultArray, int dim);
int eigengapHeuristic(double *sortedEigenValues, int len);
double **jacobiAlgorithm(double *A, int rowDim, int colDim, double **resultArray, const char *goal);
void freeJacobiAlgorithm(int *indicesIandJ, double *valsCandS, double *Ptranspose, double *Vmulti, double *firstMulti, double *P, double *Atag);
int callingJacobiFromGoal(double *symmetricMatrix, double *matrixU, const char *goal, int k, int n);
void printFromGoal(double *matrix, int rowDim, int colDim);
void stableSelectionSort(double *Array, int *indicesArray, int len);

int g_np_context = 0;
/*--------------- auxiliary functions ---------------*/

/* return euclidean distance between points x, y in d dimension */
double euclid(double x[], double y[], int dim)
{
    double sumUp = 0;
    double res;
    int coordinate;

    double *coordinateSub = (double *)malloc(dim * sizeof(double));
    if (!coordinateSub)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    for (coordinate = 0; coordinate < dim; coordinate++)
    {
        coordinateSub[coordinate] = x[coordinate] - y[coordinate];
        /* pow(coordinateSub[coordinate], 2) */
        coordinateSub[coordinate] = coordinateSub[coordinate] * coordinateSub[coordinate];
        sumUp += coordinateSub[coordinate];
    }
    /*res = sqrt(sumUp);*/
    res = sumUp;

    free(coordinateSub);

    return res;
}

/* checking whether vectors are equal */
int areCentroidsEqual(double *oldArr, double *newArr, int k, int dim)
{
    int centEqIndex;
    int cord;

    for (centEqIndex = 0; centEqIndex < k; centEqIndex++)
    {
        for (cord = 0; cord < dim; cord++)
        {
            double sub = oldArr[centEqIndex * dim + cord] - newArr[centEqIndex * dim + cord];
            if (sub > 0.00001 || sub < -0.00001)
            {
                return 0; /* false */
            }
        }
    }
	
    return 1; /* true */
}

/* step 3 - re-arrange points to clusters, after we already have cluster coordinates
    clustersCounters seems not needed */
void reArrangePointsToClusters(double *points, int *pointToWhichCluster, double *centroidsArr, int n, int k, int dim)
{
    double distance;
    double closestCentroid;
    double minDistance;
    int CurrPointIndex;
    int centIndex1;
    double disSubMin;

    for (CurrPointIndex = 0; CurrPointIndex < n; CurrPointIndex++)
    {
        minDistance = euclid(&points[CurrPointIndex * dim], &centroidsArr[0], dim);
        closestCentroid = 0;
        for (centIndex1 = 0; centIndex1 < k; centIndex1++)
        {
            distance = euclid(&points[CurrPointIndex * dim], &centroidsArr[centIndex1 * dim], dim);
            disSubMin = distance - minDistance;
            if (disSubMin < 0)
            {
                minDistance = distance;
                closestCentroid = centIndex1;
            }
        }
        pointToWhichCluster[CurrPointIndex] = closestCentroid;
    }
}

/* step 4 - re-select centroid coordinates, after we already have a full pointToWhich array
    centroidsSum should be empty
    clustersCounters should be empty */
void reSelectCentroids(double *points, int *clustersCounters, int *pointToWhichCluster, double *clusterCoordinates, double *sumCoordinates, int n, int k, int dim)
{
    int pointIndex, coordinate, clusterIndex, index;

    memset(sumCoordinates, 0, k * dim * sizeof(double));
    memset(clustersCounters, 0, k * sizeof(int));

    for (pointIndex = 0; pointIndex < n; pointIndex++)
    {
        index = pointToWhichCluster[pointIndex];

        clustersCounters[index] = clustersCounters[index] + 1;

        for (coordinate = 0; coordinate < dim; coordinate++)
        {

            sumCoordinates[index * dim + coordinate] =
                sumCoordinates[index * dim + coordinate] + points[pointIndex * dim + coordinate];
        }
    }

    for (clusterIndex = 0; clusterIndex < k; clusterIndex++)
    {
        for (coordinate = 0; coordinate < dim; coordinate++)
        {

            clusterCoordinates[clusterIndex * dim + coordinate] =
                sumCoordinates[clusterIndex * dim + coordinate] / clustersCounters[clusterIndex];
        }
    }
}

/* copy centroids from newArr to results */
void copyResults(double *newArr, int k, int dim, double *results)
{
    int resultsIndex = 0;
    int cordIndex1;
    int centIndex1;

    for (centIndex1 = 0; centIndex1 < k; centIndex1++)
    {
        for (cordIndex1 = 0; cordIndex1 < dim; cordIndex1++)
        {
            results[resultsIndex] = newArr[centIndex1 * dim + cordIndex1];
            resultsIndex++;
        }
    }
}

/* free arrays memory */
void freeArraysMemory(double *points, int *pointToWhichCluster, int *clustersCounters,
                      double *centroidsOdd, double *centroidsEven,
                      double *newCentroidAverageBeforeDivision)
{
    points++;
    points--;

    free(pointToWhichCluster);
    free(clustersCounters);
    free(centroidsOdd);
    free(centroidsEven);
    free(newCentroidAverageBeforeDivision);
}

int kmeans(int n, int k, int dim, int maxIter, double *points, double *centroids, double *results)
{
    /* indices for algorithm functions */
    int iter;
    int centroidIndex;

    /* Arrays Allocation */
    double *centroidsEven = (double *)0;
    double *centroidsOdd = (double *)0;
    int *pointToWhichCluster = (int *)0;
    int *clustersCounters = (int *)0;
    double *newCentroidAverageBeforeDivision = (double *)0;

    /* Arrays malloc & Assert */
    centroidsEven = (double *)malloc(k * dim * sizeof(double));
    if (!centroidsEven)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    centroidsOdd = (double *)malloc(k * dim * sizeof(double));
    if (!centroidsOdd)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    pointToWhichCluster = (int *)malloc(n * sizeof(int));
    if (!pointToWhichCluster)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    clustersCounters = (int *)malloc(k * sizeof(int));
    if (!clustersCounters)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    newCentroidAverageBeforeDivision = (double *)malloc(k * dim * sizeof(double));
    if (!newCentroidAverageBeforeDivision)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* --------------- Kmeans algorithm --------------- */

    /* step 1 - copy centroids to centroidsOdd */
    for (centroidIndex = 0; centroidIndex < k * dim; centroidIndex++)
    {
        centroidsOdd[centroidIndex] = centroids[centroidIndex];
    }


    /* step 3 - first arrangement of points to clusters */
    reArrangePointsToClusters(points, pointToWhichCluster, centroidsOdd, n, k, dim);

    /* step 4 - first update of centroids & check whether reached convergence */
    reSelectCentroids(points, clustersCounters, pointToWhichCluster, centroidsEven,
                      newCentroidAverageBeforeDivision, n, k, dim);

    /* copy centroids to results if reached convergence */
    if (areCentroidsEqual(centroidsOdd, centroidsEven, k, dim) == 1)
    {
        copyResults(centroidsOdd, k, dim, results);
        freeArraysMemory(points, pointToWhichCluster, clustersCounters,
                         centroidsOdd, centroidsEven, newCentroidAverageBeforeDivision);
        return 0;
    }

    /* iterate until convergence - steps 3&4 */
    for (iter = 0; iter < maxIter - 1; iter++)
    {
        if (iter % 2 == 0)
        {
            /* step 3 - arrangement of points to clusters */
            reArrangePointsToClusters(points, pointToWhichCluster, centroidsEven, n, k, dim);

            /* step 4 - update of centroids & check whether reached convergence */
            reSelectCentroids(points, clustersCounters, pointToWhichCluster, centroidsOdd,
                              newCentroidAverageBeforeDivision, n, k, dim);

            /* copy centroids to results if reached convergence */
            if (areCentroidsEqual(centroidsEven, centroidsOdd, k, dim) == 1)
            {
                copyResults(centroidsOdd, k, dim, results);
                freeArraysMemory(points, pointToWhichCluster, clustersCounters,
                                 centroidsOdd, centroidsEven, newCentroidAverageBeforeDivision);
                return 0;
            }
        }
        else
        {
            /* step 3 - arrangement of points to clusters */
            reArrangePointsToClusters(points, pointToWhichCluster, centroidsOdd, n, k, dim);

            /* step 4 - update of centroids & check whether reached convergence */
            reSelectCentroids(points, clustersCounters, pointToWhichCluster, centroidsEven,
                              newCentroidAverageBeforeDivision, n, k, dim);

            /* copy centroids to results if reached convergence */
            if (areCentroidsEqual(centroidsOdd, centroidsEven, k, dim) == 1)
            {
                copyResults(centroidsEven, k, dim, results);
                freeArraysMemory(points, pointToWhichCluster, clustersCounters,
                                 centroidsOdd, centroidsEven, newCentroidAverageBeforeDivision);
                return 0;
            }
        }
    }

    /* copy centroidsOdd if reached iter == maxIter */
    if (maxIter % 2 == 0)
    {
        copyResults(centroidsOdd, k, dim, results);
        freeArraysMemory(points, pointToWhichCluster, clustersCounters,
                         centroidsOdd, centroidsEven, newCentroidAverageBeforeDivision);
    }
    else
    { /* copy centroidsEven if reached iter == maxIter */
        copyResults(centroidsEven, k, dim, results);
        freeArraysMemory(points, pointToWhichCluster, clustersCounters,
                         centroidsOdd, centroidsEven, newCentroidAverageBeforeDivision);
    }
    return 0;
}

/* used when python calls to kmeans in c at kmeansPP algorithm */
int fitToKmeans(double *x, int n, int dim, int k, int maxIter, double *centroids, double *results)
{
    int ret = 0;

	kmeans(n, k, dim, maxIter, x, centroids, results);

    if (x == NULL || centroids == NULL || results == NULL) {
        ret = -2;
        goto out;
    }

out:
    return ret;
}

/* ######################## SP-kmeans Algorithm ######################## */

/*--------------- auxiliary functions ---------------*/

/* count commas at each line when reading file */
int commasCount(char *p)
{
    int count = 0;
    while (*p)
    {
        if (*p == ',')
            count++;
        p++;
    }
    return count;
}

/* return euclidean distance between points x, y in dim dimension */
double euclideanDistance(double x[], double y[], int dim)
{
    double sumUp = 0;
    double res;
    int coordinate;

    double *coordinateSub = (double *)malloc(dim * sizeof(double));
    if (!coordinateSub)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    for (coordinate = 0; coordinate < dim; coordinate++)
    {
        coordinateSub[coordinate] = x[coordinate] - y[coordinate];

        /* pow(coordinateSub[coordinate], 2) */
        coordinateSub[coordinate] = coordinateSub[coordinate] * coordinateSub[coordinate];
        sumUp += coordinateSub[coordinate];
    }

    /* res = sqrt(sumUp) */
    res = pow(sumUp, 0.5);

    free(coordinateSub);

    return res;
}

/* return norma of vector in dim dimension */
double vectorNorm(double *vector, int dim)
{
    int i;
    double sumOfRowsPower = 0.0;
    double norma;

    for (i = 0; i < dim; i++)
    {
        sumOfRowsPower += vector[i] * vector[i];
    }
    norma = pow(sumOfRowsPower, 0.5);

    return norma;
}

/* sub of two matrices
   matrix1 of dimension rowDim1, colDim1
   matrix2 of dimension rowDim2, colDim2
   return new matrix equal to matrix1 - matrix2 of dimension rowDim1, colDim1
   not changing  matrix1, matrix2 itself
   pre: rowDim1 = rowDim2 and colDim1 = colDim2 */
void subTwoMatrices(double *matrix1, double *matrix2, int rowDim1, int colDim1, double *subMatrix)
{
    int i;
    int j;

    for (i = 0; i < rowDim1; i++)
    {
        for (j = 0; j < colDim1; j++)
        {
            subMatrix[i * rowDim1 + j] = matrix1[i * rowDim1 + j] - matrix2[i * rowDim1 + j];
        }
    }
}

/* multiplication of two matrices
   matrix1 of dimension rowDim1, colDim1
   matrix2 of dimension rowDim2, colDim2
   return new matrix equal to matrix1 * matrix2 of dimension rowDim1, colDim2
   not changing  matrix1, matrix2 itself
   pre: colDim1 = rowDim2 */
void multiTwoMatrices(double *matrix1, double *matrix2, int rowDim1, int colDim1, int colDim2, double *multiMatrix)
{
    int i;
    int j;
    int k;
    double sumOfRowByColumnMulti = 0.0;

    for (i = 0; i < rowDim1; i++)
    {
        for (j = 0; j < colDim2; j++)
        {
            for (k = 0; k < colDim1; k++)
            {
                sumOfRowByColumnMulti += matrix1[i * rowDim1 + k] * matrix2[k * colDim1 + j];
            }

            multiMatrix[i * rowDim1 + j] = sumOfRowByColumnMulti;
            sumOfRowByColumnMulti = 0.0;
        }
    }
}

/* recieve matrix of dimension rowCol, dimRow
   return the transpose matrix of dimension dimRow, rowCol
   not changing the matrix itself */
void transpose(double *matrix, int rowDim, int colDim, double *transposeMatrix)
{
    int i;
    int j;

    for (i = 0; i < rowDim; i++)
    {
        for (j = 0; j < colDim; j++)
        {
            transposeMatrix[i * rowDim + j] = matrix[j * colDim + i];
        }
    }
}

/* creating diagonal degree matrix D from Weight matrix W
   diagonal[i,j] = sum of weight[i,z] when z goes from 1 to n
   recieves empty matrix of size n*n and fill it
   not changing weightedMatrix itself */
void diagonalMatrixFromW(double *weightedMatrix, int len, double *diagonal)
{
    int i;
    int z;
    double rowSum = 0.0;

#if 0
   for (i = 0; i < len; i++)
		for (z = 0; z < len; z++)
			diagonal[i*len+z] = 0.0;
#endif
	
    for (i = 0; i < len; i++)
    {
        for (z = 0; z < len; z++)
        {
            rowSum += weightedMatrix[i * len + z];
        }
        diagonal[i * len + i] = rowSum;
        rowSum = 0.0;
    }
}

/* pow(-1/2) of a diagonal matrix D
   return new matrix, not changing D itself
   pre: matrix D must be diagonal dim * dim order */
void negativeSqrtMatrix(double *matrixD, int dim, double *negativeSqrtMatrix)
{
    int i;

    for (i = 0; i < dim; i++)
    {
        /* in case eleme is zero, avoid dividing by zero */
        if (matrixD[i * dim + i] == 0.0)
        {
            negativeSqrtMatrix[i * dim + i] = 0.0;
        }
        /* otherwise */
        else
        {
            negativeSqrtMatrix[i * dim + i] = pow(matrixD[i * dim + i], -0.5);
        }
    }
}

/* creating Weighted Adjacency Matrix matrix W.
   recieve points array, n, dim and empty matrix weightedAdjancy
   return the matrix weightedAdjancy of size n*n after filling it
   not changing points itself */
void WeightedAdjacencyMatrixFromPoints(double *points, int n, int dim, double *weightedAdjancy)
{
    int i;
    int j;
    double euclidDistance;

    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            /* weights are symmetric - so euclidean distance calculates once */
            euclidDistance = euclideanDistance(&points[i * dim], &points[j * dim], dim);
            euclidDistance = exp(euclidDistance * -0.5);
            weightedAdjancy[i * n + j] = euclidDistance;
            weightedAdjancy[j * n + i] = euclidDistance;
        }
    }
}

/* form new matrix from matrixU by renormalizing each of U’s rows to have unit length
   return pointer to new matrix renormalizedRowMatrix, not changing U itself
   renormalizedRowMatrix is empty of size rowDim * colDim */
void renormalizingMatrixRow(double *matrixU, int rowDim, int colDim, double *renormalizedRowMatrix)
{
    int i;
    int j;
    int l;
    double iRowNorma;
    double *iRow;
	
    iRow = (double *)calloc(colDim, sizeof(double));
    if (!iRow)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }
	
	
    for (i = 0; i < rowDim; i++)
    {
        /* find the norma of row i */
        for (l = 0; l < colDim; l++)
        {
            iRow[l] = matrixU[i * colDim + l];
        }
        iRowNorma = vectorNorm(iRow, colDim);
		
		/* update T by taking U and renormalize */
        for (j = 0; j < colDim; j++)
        {
            if (iRowNorma != 0.0)
            {
                renormalizedRowMatrix[i * colDim + j] = matrixU[i * colDim + j] / iRowNorma;
            }
            else
            {
                renormalizedRowMatrix[i * colDim + j] = 0.0;
            }
        }
    }
    free(iRow);
}

/* creating the normalized graph Laplacian
   recieve D^(-0.5) and W matrices
   calculate Lnorm = I - D^(-0.5) * W * D^(-0.5)
   return the matrix Lnorm, not changing W, D^(-0.5)
   Lnorm is of size dim * dim */
void normalizedGraphLaplacian(double *Dsqrt, double *W, int dim, double *Lnorm)
{
    int i;
    double *identityMatrix;
    double *firstMulti;
    double *secondMulti;

    /* identity matrix of dimension dim, dim */
    identityMatrix = (double *)calloc(dim * dim, sizeof(double));
    if (!identityMatrix)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }
    for (i = 0; i < dim; i++)
    {
        identityMatrix[i * dim + i] = 1;
    }

    /* matrix for multiplication */
    firstMulti = (double *)calloc(dim * dim, sizeof(double));
    if (!firstMulti)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* matrix for multiplication */
    secondMulti = (double *)calloc(dim * dim, sizeof(double));
    if (!secondMulti)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    multiTwoMatrices(Dsqrt, W, dim, dim, dim, firstMulti);
    multiTwoMatrices(firstMulti, Dsqrt, dim, dim, dim, secondMulti);
    subTwoMatrices(identityMatrix, secondMulti, dim, dim, Lnorm);

    free(identityMatrix);
    free(firstMulti);
    free(secondMulti);
}

/* receive matrix and return indices i, j.
   where Aij is largest absolute off-diagonal value in matrix */
void IndexOfMaxValInMatrix(double *matrix, int rowDim, int colDim, int *resultIndices)
{
    int i;
    int j;
    int maxI = 0;
    int maxJ = 0;
    double maxVal = -1.0;

    for (i = 0; i < rowDim; i++)
    {
        for (j = 0; j < colDim; j++)
        {
            if (i != j)
            {

                if (maxVal < fabs(matrix[i * colDim + j]))
                {
                    maxVal = fabs(matrix[i * colDim + j]);
                    maxI = i;
                    maxJ = j;
                }
            }
        }
    }
    resultIndices[0] = maxI;
    resultIndices[1] = maxJ;
}

/* receive matrix A, dim of matrix A and indices i, j
   return c, s values according to formula in project notes
   c and s needed for buliding the rotation matrix P */
void computingCandSforRotationMatrix(double *A, int rowDim, int *indices, double *valsCandS)
{
    int i = indices[0];
    int j = indices[1];
    int sign;
    double theta;
    double t;
    double c;
    double s;

    /* computing theta & theta sign according to project notes */
    theta = (A[j * rowDim + j] - A[i * rowDim + i]) / (2 * A[i * rowDim + j]);
    if (theta > 0.0 || theta == 0.0)
    {
        sign = 1;
    }
    else
    {
        sign = -1;
    }

    /* computing t according to project notes */
    t = (sign) / (fabs(theta) + pow(pow(theta, 2) + 1, 0.5));

    /* computing c according to project notes */
    c = 1 / (pow(pow(t, 2) + 1, 0.5));

    /* computing s according to project notes */
    s = t * c;

    valsCandS[0] = c;
    valsCandS[1] = s;
}

/* calculate sum of squares of all off-diagonal elements in matrix */
double sumSquareOffDiagonalElem(double *matrix, int rowDim, int colDim)
{
    int i;
    int j;
    double sumOfSquares = 0;

    for (i = 0; i < rowDim; i++)
    {
        for (j = 0; j < colDim; j++)
        {
            if (i != j)
            {
                sumOfSquares += pow(matrix[i * colDim + j], 2);
            }
        }
    }

    return sumOfSquares;
}

/* find eigenvalues.
   sort eigenvalues & eigenVectors indices.
   if k=0 find k with Eigengap Heuristic.
   update eigenvalues in return array of jacobiAlgorithm.
   update k in return array of jacobiAlgorithm.
   update matrix U of eigenvectors in return array of jacobiAlgorithm. */
void eigenvaluesSortingAndHeuristic(double *Atag, double *V, double **resultArray, int dim)
{
    int i;
    int j;
    int k;
    int *eigenVectorsIndices;

    eigenVectorsIndices = (int *)calloc(dim, sizeof(int));
    if (!eigenVectorsIndices)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* find eigenvalues */
    for (i = 0; i < dim; i++)
    {
        resultArray[1][i] = Atag[i * dim + i];
    }

    /* sort eigenvalues & eigenVectors indices */
    /* ------- this function update eigenvalues in return array of jacobiAlgorithm ------- */
    stableSelectionSort(resultArray[1], eigenVectorsIndices, dim);

    k = resultArray[2][0];

    /* if k=0 find k with Eigengap Heuristic */
    if (k == 0)
    {
        k = eigengapHeuristic(resultArray[1], dim);
    }

    /* ------- update k in return array of jacobiAlgorithm ------- */
    resultArray[2][0] = (double)k;

    /* changes U by taking k first eigenvectors columns from matrix V */
    /* ------- it updates matrix U in return array of jacobiAlgorithm ------- */
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < k; j++)
        {
            resultArray[0][i * dim + j] = V[i * dim + eigenVectorsIndices[j]];
        }
    }

    free(eigenVectorsIndices);
}

/* stable selection sort.
    used for eigenvalues sort, meaning Array = eigenvalues array.
    at same time sort eigenVectors indices array. */
void stableSelectionSort(double *Array, int *indicesArray, int len)
{
    int i;
    int j;
    int minIndex;
    double minValue;
    double *originEigenValue;

    originEigenValue = (double *)calloc(len, sizeof(double));
    if (!originEigenValue)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }
    for (i = 0; i < len; i++)
    {
        originEigenValue[i] = Array[i];
    }

    /* loop invariant: elements till Array[i - 1] are already sorted */
    for (i = 0; i < len; i++)
    {
        /* find minimum element from Array[i] to Array[n - 1] */
        minIndex = i;
        for (j = i + 1; j < len; j++)
            if (Array[minIndex] > Array[j])
                minIndex = j;

        /* move elements from i to minIndex one step right */
        minValue = Array[minIndex];
        while (minIndex > i)
        {
            Array[minIndex] = Array[minIndex - 1];
            minIndex--;
        }
        /* put minimum element at current index i */
        Array[i] = minValue;
    }

    for (i = 0; i < len; i++)
    {
        for (j = 0; j < len; j++)
        {
            if (Array[i] == originEigenValue[j])
            {
                indicesArray[i] = j;
                originEigenValue[j] = -1.0; /* in spk all eigenvalues of Lnorm is non-negative*/
                break;
            }
        }
    }
    free(originEigenValue);
}

/* computing k in case k=0 */
int eigengapHeuristic(double *sortedEigenValues, int len)
{
    int i;
    double tmpMax;
    int k = 0;
    double max = fabs(sortedEigenValues[0] - sortedEigenValues[1]);

    for (i = 0; i < (len / 2); i++)
    {
        tmpMax = fabs(sortedEigenValues[i] - sortedEigenValues[i + 1]);
        if (max < tmpMax)
        {
            max = tmpMax;
            k = i + 1;
        }
    }

    return k;
}

/* free jacobiAlgorithm function arrays */
void freeJacobiAlgorithm(int *indicesIandJ, double *valsCandS, double *Ptranspose, double *Vmulti,
                         double *firstMulti, double *P, double *Atag)
{
    free(indicesIandJ);
    free(valsCandS);
    free(Ptranspose);
    free(Vmulti);
    free(firstMulti);
    free(P);
    free(Atag);
}

/* The Jacobi eigenvalue algorithm is an iterative method,
   calculate the eigenvalues and eigenvectors of a real symmetric matrix (diagonalization) */
double **jacobiAlgorithm(double *A, int rowDim, int colDim, double **resultArray, const char *goal)
{
    /* function variables */
    int i;
    int j;
    int iterationCnt;
    double epsilon = pow(10, -15);
    double offA;
    double offAtag;

    /* array for the function that calculate i, j */
    int *indicesIandJ;
    int maxI;
    int maxJ;

    /* array for the function that calculate c, s */
    double *valsCandS;
    double c;
    double s;

    /* arrays for the function */
    double *P;
    double *Ptranspose;
    double *Atag;
    double *Vmulti;
    double *firstMulti;

    indicesIandJ = (int *)calloc(2, sizeof(int));
    if (!indicesIandJ)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    valsCandS = (double *)calloc(2, sizeof(double));
    if (!valsCandS)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* creating P as diagonal matrix of dimension rowCol, dimRow */
    P = (double *)calloc(rowDim * colDim, sizeof(double));
    if (!P)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* at first P is identity matrix */
    for (i = 0; i < rowDim; i++)
    {
        P[i * rowDim + i] = 1.0;
    }

    /* creating matrix for transpose function */
    Ptranspose = (double *)calloc(rowDim * colDim, sizeof(double));
    if (!Ptranspose)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* creating Atag as matrix of dimension rowCol, dimRow */
    Atag = (double *)calloc(rowDim * colDim, sizeof(double));
    if (!Atag)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* creating matrix for multiplication function */
    firstMulti = (double *)calloc(rowDim * colDim, sizeof(double));
    if (!firstMulti)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* creating matrix for multiplication function */
    Vmulti = (double *)calloc(rowDim * colDim, sizeof(double));
    if (!Vmulti)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    /* at first Vmulti is identity matrix */
    for (i = 0; i < rowDim; i++)
    {
        Vmulti[i * rowDim + i] = 1.0;
    }

    /* ----- repeat a, b steps until A' is diagonal matrix ----- */

    /* stop if iteration amount is bigger than 100 */
    iterationCnt = 0;
    while (iterationCnt < 100)
    {
        /* geting indices of absolute off-diagonal value in matrix A */
        IndexOfMaxValInMatrix(A, rowDim, colDim, indicesIandJ);
        maxI = indicesIandJ[0];
        maxJ = indicesIandJ[1];

        /* getting c, s according to the formula in project notes */
        computingCandSforRotationMatrix(A, rowDim, indicesIandJ, valsCandS);
        c = valsCandS[0];
        s = valsCandS[1];

        /* reset P to identity matrix */
        memset(P, 0, sizeof(double) * rowDim * colDim);
        for (i = 0; i < rowDim; i++)
        {
            P[i * rowDim + i] = 1.0;
        }

        /* build rotation matrix P - put vals in relevent indices */
        P[maxI * rowDim + maxJ] = s;
        P[maxI * rowDim + maxI] = c;
        P[maxJ * rowDim + maxI] = -s;
        P[maxJ * rowDim + maxJ] = c;

        /* Transform the matrix A: Atag = P-transpose * A * P */
        transpose(P, rowDim, colDim, Ptranspose);
        multiTwoMatrices(Ptranspose, A, rowDim, colDim, colDim, firstMulti);
        multiTwoMatrices(firstMulti, P, rowDim, colDim, colDim, Atag);

        /* Calculate the eigenvectors of A by multiplying all the rotation matrices */
        /* update eigenvectors matrix V in return array */
        multiTwoMatrices(Vmulti, P, rowDim, colDim, colDim, resultArray[0]);
        for (i = 0; i < rowDim; i++)
        {
            for (j = 0; j < colDim; j++)
            {
                Vmulti[i * rowDim + j] = resultArray[0][i * rowDim + j];
            }
        }

        /* stop if reach convergence */
        offA = sumSquareOffDiagonalElem(A, rowDim, colDim);
        offAtag = sumSquareOffDiagonalElem(Atag, rowDim, colDim);
        if (fabs(offA - offAtag) < epsilon || fabs(offA - offAtag) == epsilon)
        {
            break;
        }
        /* otherwise A = Atag */
        memcpy(A, Atag, sizeof(double) * rowDim * colDim);
        iterationCnt++;
    }

    /* in case A is diagonal or there was a convergence or iterationCnt is 100 */

    if (!strcmp(goal, "spk"))
    {
        /* in case of spk - sort eigenvalues, doing heuristic
        in case k=0 and update result array */
        eigenvaluesSortingAndHeuristic(Atag, Vmulti, resultArray, rowDim);
    }
    else if (!strcmp(goal, "jacobi"))
    {
        for (i = 0; i < rowDim; i++)
        {
            /* update eigenvectors in result array */
            for (j = 0; j < colDim; j++)
            {
                resultArray[0][i * rowDim + j] = Vmulti[i * rowDim + j];
            }

            /* update eigenvalues in resultArray */
            resultArray[1][i] = Atag[i * rowDim + i];
        }
    }
	
    freeJacobiAlgorithm(indicesIandJ, valsCandS, Ptranspose, Vmulti, firstMulti, P, Atag);

    return resultArray;
}

/* Calling to jacobiAlgorithm function,
   in case goal is jacobi the symmetricMatrix is points.
   in case goal is spk the symmetricMatrix is Lnorm.
   Creating return Array for jacobiAlgorithm,
   jacobiResArray[0] = eigenVectors matrix U of size n*n.
   jacobiResArray[1] = eigenValues matrix of size n.
   jacobiResArray[2] = matrix of size 1, will contain k. */
int callingJacobiFromGoal(double *symmetricMatrix, double *matrixU, const char *goal, int k, int n)
{
    int i;
    int j;
    double **jacobiResArray;
    double *V;
    double *Vtranspose;
    double *eigenvalues;
    double *kArray;

    jacobiResArray = (double **)calloc(3, sizeof(double *));
    if (!jacobiResArray)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    V = (double *)calloc(n * n, sizeof(double));
    if (!V)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    Vtranspose = (double *)calloc(n * n, sizeof(double));
    if (!Vtranspose)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    eigenvalues = (double *)calloc(n, sizeof(double));
    if (!eigenvalues)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    kArray = (double *)calloc(1, sizeof(double));
    if (!kArray)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }

    kArray[0] = (double)k;

    jacobiResArray[0] = V;
    jacobiResArray[1] = eigenvalues;
    jacobiResArray[2] = kArray;

    jacobiAlgorithm(symmetricMatrix, n, n, jacobiResArray, goal);

    if (!strcmp(goal, "jacobi"))
    {
        /* we print by rows and eigenvectors are columns of matrix */
        transpose(V, n, n, Vtranspose);

        /* fill matrix U with k first eigenvectors
           we run till n+1 for placing eigenvalues on first row */
        for (i = 0; i < n + 1; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (i == 0)
                {
                    matrixU[i * n + j] = eigenvalues[i * n + j];
                }
                else
                {
                    matrixU[i * n + j] = Vtranspose[(i - 1) * n + j];
                }
            }
        }
    }
    else if (!strcmp(goal, "spk"))
    {
        if (k == 0)
        {
            /* taking real k - needed in case k was 0 */
            k = (int)jacobiResArray[2][0];
        }

        /* fill matrix U with k first eigenvectors */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < k; j++)
            {
                matrixU[i * k + j] = V[i * n + j];
            }
        }
    }

    free(V);
    free(Vtranspose);
    free(eigenvalues);
    free(kArray);
    free(jacobiResArray);

    return k;
}

/* print matrix - using in main. */
void printFromGoal(double *matrix, int rowDim, int colDim)
{
    int i;
    int j;
	
	if (g_np_context)
		return;

    for (i = 0; i < rowDim; i++)
    {
        for (j = 0; j < colDim; j++)
        {
            if (j != colDim - 1)
            {
                if (matrix[i * colDim + j] < 0 && matrix[i * colDim + j] > -0.0001)
                { /* in case number is -0.0000 */
                    printf("%.4f", 0.0000);
                    printf(",");
                }
                else
                {
                    printf("%.4f", matrix[i * colDim + j]);
                    printf(",");
                }
            }
            else
            {
                if (matrix[i * colDim + j] < 0 && matrix[i * colDim + j] > -0.0001)
                { /* in case number is -0.0000 */
                    printf("%.4f", 0.0000);
                }
                else
                {
                    printf("%.4f", matrix[i * colDim + j]);
                }
            }
        }
        printf("\n");
    }
}


/*--------------- main program of SPkmeans ---------------*/

int do_main(const char *arg, const char *goal, const char *fileName, double *results_array, double *k_array)
{
	    /* using in atoi when scanning args */
    double kAtof;

    /* in use at scan of input file for points array */
    FILE *file;
    char line[1000]; /* up to 1000 points and up to 10 features */
    int n;
    int dim;
    int k = 0;
    int nCounter;
    int dimCounter;
    char *token;

    /* algorithm variables */
    int i;

    /* arrays pointers */
    double *points = (double *)0;
    double *weighted = (double *)0;
    double *diagonal = (double *)0;
    double *negativeSqrtDiag = (double *)0;
    double *Lnorm = (double *)0;
    double *U = (double *)0;
    double *T = (double *)0;
    double *kmeansResult = (double *)0;
    double *initCentroids = (double *)0;

    /*--------------- saving arguments ---------------*/
	if (results_array)
		g_np_context = 1;

	if (strcmp(goal, "spk") && strcmp(goal, "wam") && strcmp(goal, "ddg") &&
		strcmp(goal, "lnorm") && strcmp(goal, "jacobi"))
	{ /* goal is not valid */
		printf("%s\n", "Invalid Input!");
		assert(0);
	}

	k = atoi(arg);
	kAtof = atof(arg);

	if (!(strcmp(goal, "spk")))
	{
		if ((double)(k) != kAtof)
		{ /* Invalid input. K must be integer */
			printf("%s\n", "Invalid Input!");
			assert(0);
		}

		if (k < 0)
		{ /* Invalid input. K must be non negative integer */
			printf("%s\n", "Invalid Input!");
			assert(0);
		}
	}
 

    /*--------------- scanning file ---------------*/

	n = 0;
	dim = 0;
	
    file = fopen(fileName, "r");
    if (file == NULL)
    { /* Error while opening the file */
        printf("%s\n", "An Error Has Occured");
    }

    /* first scan for getting n, dim */
    while (fgets(line, sizeof(line), file))
    {
        n++;
        if (dim == 0)
			dim = commasCount(line) + 1;
    }


    /* check correctness of n, dim  */

    if (n < 1)
    { /* Invalid input. At least one point is needed. */
        printf("%s\n", "Invalid Input!");
        assert(0);
    }

    if ((n <= k) && !(strcmp(goal, "spk")))
    { /* Invalid input. K must be smaller than N and not equal to N. */
        printf("%s\n", "Invalid Input!");
        assert(0);
    }

    if (dim == 0)
    { /* Invalid input. Dim of points is zero. */
        printf("%s\n", "Invalid Input!");
        assert(0);
    }

    /* create points array at exact length of points amount  */
    points = (double *)malloc(n * dim * sizeof(double));
    if (!points)
    { /* Allocation memory problem. */
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }


    /* fill points array */
    nCounter = 0;
    rewind(file); /* getting to head of file for re-scanning */
	while (fgets(line, sizeof(line), file))
    {
        /* get first token */
        token = strtok(line, ",");

        /* walk through other tokens */
		dimCounter = 0;
        while (token != NULL)
        {
			
			points[nCounter*dim + dimCounter] = atof(token);
			
            dimCounter += 1;
            token = strtok(NULL, ",");
        }
        nCounter += 1;
    }
    fclose(file);

    /*--------------- divided into goals cases ---------------*/

    if (!strcmp(goal, "jacobi"))
    {
        /* creating matrix of first k eigenvectors */
        U = (double *)calloc((n + 1) * n, sizeof(double));
        if (!U)
        {
            printf("%s\n", "An Error Has Occured");
            assert(0);
        }

        /* k matters only in spk case */
        callingJacobiFromGoal(points, U, goal, k, n);
        printFromGoal(U, n + 1, n);

		if (results_array) {
			memcpy(results_array, U, sizeof(double) * (n+1)*n);
		}


        free(points);
        free(U);

        return 0;
    }

    /* creating matrix for weighted adjacency function */
    weighted = (double *)calloc(n * n, sizeof(double));
    if (!weighted)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }
    WeightedAdjacencyMatrixFromPoints(points, n, dim, weighted);

    if (!strcmp(goal, "wam"))
    {
        printFromGoal(weighted, n, n);

		if (results_array) {
			memcpy(results_array, weighted, sizeof(double) * (n)*n);
		}


        free(points);
        free(weighted);

        return 0;
    }

    /* creating matrix for diagonalization function */
    diagonal = (double *)calloc(n * n, sizeof(double));
    if (!diagonal)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }
    diagonalMatrixFromW(weighted, n, diagonal);

    if (!strcmp(goal, "ddg"))
    {
        printFromGoal(diagonal, n, n);

		if (results_array) {
			memcpy(results_array, diagonal, sizeof(double) * n*n);
		}

        free(points);
        free(weighted);
        free(diagonal);

        return 0;
    }

    /* creating matrix for negative sqrt of diagonal function */
    negativeSqrtDiag = (double *)calloc(n * n, sizeof(double));
    if (!negativeSqrtDiag)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }
    negativeSqrtMatrix(diagonal, n, negativeSqrtDiag);

    /* creating matrix for normalized graph laplacian function */
    Lnorm = (double *)calloc(n * n, sizeof(double));
    if (!Lnorm)
    {
        printf("%s\n", "An Error Has Occured");
        assert(0);
    }
    normalizedGraphLaplacian(negativeSqrtDiag, weighted, n, Lnorm);

    if (!strcmp(goal, "lnorm"))
    {
        printFromGoal(Lnorm, n, n);

		if (results_array) {
			memcpy(results_array, Lnorm, sizeof(double) * n * n);
		}
		
        free(points);
        free(weighted);
        free(diagonal);
        free(negativeSqrtDiag);
        free(Lnorm);

        return 0;
    }

    if (!strcmp(goal, "spk"))
    {
        /* creating matrix of first k eigenvectors */
        U = (double *)calloc(n * n, sizeof(double));
        if (!U)
        {
            printf("%s\n", "An Error Has Occured");
            assert(0);
        }

        k = callingJacobiFromGoal(Lnorm, U, goal, k, n);

        /* creating matrix for normalize eigenvectors matrix */
        T = (double *)calloc(n * k, sizeof(double));
        if (!T)
        {
            printf("%s\n", "An Error Has Occured");
            assert(0);
        }

        renormalizingMatrixRow(U, n, k, T);
		
		/*
		*result_to_save = (double*)malloc(sizeof(double) * n * k);
		*/
		if (results_array) {
			memcpy(results_array, T, sizeof(double) * n * k);
		}
		
		if (k_array) {
			k_array[0] = k;
		}


        /* each row of T is a point in Rk, cluster them into k clusters
        via the K-means algorithm. maxIter = 300. */

        /* creating matrix for first k centroids */
        initCentroids = (double *)calloc(k * k, sizeof(double));
        if (!initCentroids)
        {
            printf("%s\n", "An Error Has Occured");
            assert(0);
        }

        /* init first k centroids for kmeans */
        for (i = 0; i < k*k; i++)
        {
            initCentroids[i] = T[i];
        }

        /* creating matrix for kmeans result */
        kmeansResult = (double *)calloc(k * k, sizeof(double));
        if (!kmeansResult)
        {
            printf("%s\n", "An Error Has Occured");
            assert(0);
        }

        kmeans(n, k, k, 300, T, initCentroids, kmeansResult);

        /* printing  */
		printFromGoal(kmeansResult, k, k);

        /* free & return */
        free(points);
        free(weighted);
        free(diagonal);
        free(negativeSqrtDiag);
        free(Lnorm);
        free(U);
        free(T);
        free(initCentroids);
        free(kmeansResult);

        return 0;
    }
    return 0;
}

int main(int argc, char *argv[])
{
	/* using in atoi when scanning args */
    const char *arg = 0;
    const char *fileName = 0;
    const char *goal = 0;
	int n=0, dim=0;
    FILE *file;
    char line[1000]; /* up to 1000 points and up to 10 features */

    /*--------------- saving arguments ---------------*/

    if (argc == 4)
    {
        goal = argv[2];
        arg = argv[1];
        fileName = argv[3];
    }
    else
    { /* Amount of arguments is invalid */
        printf("%s\n", "Invalid Input!");
        assert(0);
    }
	
	
	file = fopen(fileName, "r");
    while (fgets(line, sizeof(line), file))
    {
        n++;
        if (dim == 0)
			dim = commasCount(line) + 1;
    }
	fclose(file);

	return do_main(arg, goal, fileName, 0, 0);
}


