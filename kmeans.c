#include "spkmeans.h"


int kMeans(int N, int K, double **Datapoints, double **Centroids, int dimension){

    int i, j, counter, cluster;
    double **used_cent;
    int max_iter = 300;
    counter = 0;

    used_cent = malloc((sizeof(double *)) * K);
    if (used_cent == NULL)
        return -1;
    
    for(i=0;i<K;i++){
        used_cent[i] = malloc((sizeof(double)) * (dimension));
        if (used_cent[i] == NULL)
            return -1;
        for(j=0; j<dimension; j++)
            used_cent[i][j]=Centroids[i][j];
        
    }

    /* Calculate k-means:
    Stop iteration when number of iteration is more then man_iter
    or when all of the centroids have changed less then epsilon*/
    
    while (max_iter > counter){
        for (i = 0; i < N; i++){
            cluster = find_cluster(Centroids, Datapoints[i], dimension, K);

            /* Update for each datapoint- what number of cluster it belongs
            and for each centroids update the number of datapoints that belong to it*/
            Datapoints[i][dimension] = cluster;
            Centroids[cluster - 1][dimension] += 1;
        }

        /* Set all centroids to zero so that updated centroids can be calculated next*/
        for (i = 0; i < K; i++){
            for (j = 0; j < dimension; j++)
                Centroids[i][j] = 0;
        }

        /* Update centroids according to the calculations*/
        for (i = 0; i < N; i++){
            cluster = Datapoints[i][dimension];
            for (j = 0; j < dimension; j++)
                Centroids[cluster - 1][j] += Datapoints[i][j];  
        }
        for (i = 0; i < K; i++){
            for (j = 0; j < dimension; j++){
                /*Centroids[i][j]= sum of datapoints that belong to cluster i+1
                Centroids[i][dimension]= number of datapoints that belong to cluster i*/
                Centroids[i][j] = Centroids[i][j] / Centroids[i][dimension];
            }
            Centroids[i][dimension] = 0;
        }

        /* if all centroids changed less then epsilon -done, else- countinue*/
        if (check_euclidean_norm(Centroids, used_cent, dimension, K))
            break;

        /* make oldcentroids be the new ones for next iteration*/
        update_old_centroids(Centroids, used_cent, dimension, K);

        counter++;
    }

    free_memory(used_cent, K);
    return SUCCESS;
}

/* Receives the new and old centroids
 * Returns 1 if all of the centroids didn't change more then epsilon,else-0*/
int check_euclidean_norm(double **newCentroids, double **used_cent, int dimension, int K){
    int i;
    double e_norm;

    /*Calculate euclidean norm for each centroid*/
    for (i = 0; i < K; i++){
        e_norm = calc_euclidean_norm(newCentroids[i],used_cent[i],dimension);
        /* One centroid changed more then epsilon*/
        if (e_norm >= 0)
            return 0;
    }
    /* Every centroids changed less then epsilon */
    return 1;
}

/* Receives the centroids and one datapoint
 * Returns datapoint's cluster*/
int find_cluster(double **Centroids, double *Datapoint, int dimension, int K){
    int i, j, cluster;
    double sum, min_clust;

    cluster=0; /*Default*/

    min_clust = DBL_MAX;
    for (i = 0; i < K; i++){
        sum = 0;
        for (j = 0; j < dimension; j++)
            sum += pow((Datapoint[j] - Centroids[i][j]), 2);
        if (min_clust > sum){
            min_clust = sum;
            cluster = i + 1;
        }
    }

    return cluster;
}

/* Receives the updated centroids and old ones- update the old centroids*/
void update_old_centroids(double **newCentroids, double **used_cent, int dimension, int K){
    int i, j;
    for (i = 0; i < K; i++){
        for (j = 0; j < dimension; j++)
            used_cent[i][j] = newCentroids[i][j];
        
    }
}
