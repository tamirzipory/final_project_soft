#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

int fitToKmeans(double *points, int n, int dim, int k, int maxIter, double *centroids, double *results);
int do_main(const char *arg, const char *goal, const char *fileName, double *results_array, double *k_array);

