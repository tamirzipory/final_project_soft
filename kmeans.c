#include "spkmeans.h"

int kMeans(int len, int in1, double **points, double **cent, int dim){
    int i, j, count, clust;
    double **cent_that_used;
    int max_iter = 300;
    count = 0;
    cent_that_used = calloc(in1 , sizeof(double *));
    if (cent_that_used == NULL)
        return -1;
    i = 0;
    while(i < in1){
        cent_that_used[i] = calloc(dim , sizeof(double));
        if (cent_that_used[i] == NULL)
            return -1;
        j = 0;
        while(j < dim){
           cent_that_used[i][j]=cent[i][j];
           j++;
        }
        i++;
    }

    /* Calculate k-means:
    Stop iteration when number of iteration is more then man_iter
    or when all of the centroids have changed less then epsilon*/
    
    while (max_iter > count){
        i = 0;
        while(i < len){
            clust = find_cluster(cent, points[i], dim, in1);
            points[i][dim] = clust;
            cent[clust - 1][dim]++;
            i++;
        }
        i = 0;
        while(i < in1){
            for (j = 0; j < dim; j++)
                cent[i][j] = 0;
            i++;
        }
        i = 0;
        while(i < len){
            clust = points[i][dim];
            for (j = 0; j < dim; j++)
                cent[clust - 1][j] = points[i][j] + cent[clust - 1][j];  
            i++;
        }
        i = 0;
        while(i < in1){
            for (j = 0; j < dim; j++)
                cent[i][j] = cent[i][j] / cent[i][dim];
            cent[i][dim] = 0;
            i++;
        }
        if (check_euclidean_norm(cent, cent_that_used, dim, in1))
            break;
        update_old_centroids(cent, cent_that_used, dim, in1);
        count++;
    }

    free_memory(cent_that_used, in1);
    return 0;
}

int check_euclidean_norm(double **newcent, double **cent_that_used, int dim, int in1){
    int i;
    double e_norm;
    for (i = 0; i < in1; i++){
        e_norm = calc_euclidean_norm(newcent[i],cent_that_used[i],dim);
        if (e_norm >= 0)
            return 0;
    }
    return 1;
}


int find_cluster(double **cent, double *Datapoint, int dim, int K){
    int i, j, cluster;
    double sum, min_clust;

    cluster=0; 

    min_clust = DBL_MAX;
    for (i = 0; i < K; i++){
        sum = 0;
        for (j = 0; j < dim; j++)
            sum += pow((Datapoint[j] - cent[i][j]), 2);
        if (min_clust > sum){
            min_clust = sum;
            cluster = i + 1;
        }
    }

    return cluster;
}

void update_old_centroids(double **newcent, double **cent_that_used, int dim, int K){
    int i, j;
    for (i = 0; i < K; i++){
        for (j = 0; j < dim; j++)
            cent_that_used[i][j] = newcent[i][j];
        
    }
}
