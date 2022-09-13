#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "kmeans.h"

double* Kmeans(int kk, int n, int d, int max_it, double epss, double* dp, double *cent)
{
    k = kk;
    max_iter = max_it;
    cluster_length = d;
    points_counter = n*d;
    datapoints = dp;
    centroids = cent;
    eps = epss;
    cluster_points = (double *)malloc(k * cluster_length * sizeof(double));
    to_avg = (int *)malloc(k * sizeof(int));
    if (cluster_points == NULL || to_avg == NULL)
    {
        printf("An Error Has Occurred");
        free(centroids);
        freePointers();
        exit(1);
    }

    /* ----- The actual algorithm ----- */
    iterates = max_iter;
    is_norm = 1;
    do
    { /* Reseting to_avg */
        for (i = 0; i < k; i++)
            to_avg[i] = 0;

        /* Reseting cluster_points */
        for (i = 0; i < k * cluster_length; i++)
            cluster_points[i] = 0;

        /* ----- Assigning Xi to the closest Sj ----- */
        for (i = 0; i < points_counter; i += cluster_length)
        { /* Xi loop */

            /* ----- Finding the closest centroid ----- */
            min_dist = DBL_MAX;
            assign_to_cluster = 0;

            for (j = 0, p = 0; j < k * cluster_length; j += cluster_length, p++)
            {
                distance = 0;
                for (m = 0; m < cluster_length; m++)
                {
                    /* Cell m in dist = (Xi_m-Cj_m)^2 ====== The distance from centroid Cj */
                    distance += (datapoints[i + m] - centroids[j + m]) * (datapoints[i + m] - centroids[j + m]);
                }
                if (distance < min_dist)
                {
                    min_dist = distance;
                    assign_to_cluster = p;
                }
            }
            ++to_avg[assign_to_cluster];

            /* Assignin Xi to the closest centroid */
            for (t = 0; t < cluster_length; t++)
            {
                cluster_points[(assign_to_cluster * cluster_length) + t] += datapoints[i + t];
            }
        }
        /* ----- Updating the centroids ----- */
        for (i = 0, p = 0; i < k * cluster_length; i += cluster_length, p++)
        {
            for (j = 0; j < cluster_length; j++)
            {
                /* Calculating the avarage */
                if (to_avg[p] > 0)
                    cluster_points[i + j] /= to_avg[p];
            }
        }

        /* Now cluster_points is the new centroids array */
        normDelta = 0;
        /* Calculating the norms */
        for (i = 0; i < k; i++)
        {
            for (j = 0; j < k * cluster_length; j += cluster_length)
            {
                normDelta += pow((cluster_points[i + j] - centroids[i + j]), 2);
            }
        }
        /* Copying cluster_points which are the new centroids to centroids */
        for (i = 0; i < k * cluster_length; i++)
        {
            centroids[i] = cluster_points[i];
        }

        /* Checking if norm is breaking the condition */
        if (normDelta < eps)
            is_norm = 0;

    } while (iterates-- > 0 && is_norm == 1);

    freePointers();

    return centroids;
}



void freePointers() {
    free(cluster_points);
    free(datapoints);
    free(to_avg);
}
