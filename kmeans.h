#ifndef KMEANS_H_
#define KMEANS_H_


/* declaring variables */
int k, max_iter, iterates, *to_avg, points_counter, cluster_length, assign_to_cluster;
int i, j, m, p, t, is_norm;
double *centroids, *datapoints, *cluster_points, temp_point, min_dist, distance, normDelta, eps;

/* declaring functions */

void freePointers();
double* Kmeans(int k, int n, int d, int max_iter, double epss, double* dp, double *c);

#endif
