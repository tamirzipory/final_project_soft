#include "spkmeans.h"



int checkTheNorm(double **centNew, double **cent_that_used, int dim, int in1){
    int i;
    double e;
    i = 0;
    while (i < in1){
        e = calc_euclidean_norm(centNew[i] , cent_that_used[i], dim);
        if (e >= 0)
            return 0;
        i++;
    }
    return 1;
}


int assign_cluster(double **cent, double *data, int dim, int K){
    int i, j, ret;
    double sum, min;
    ret=0; 
    min = DBL_MAX;
    i = 0;
    while (i < K)
    {
        sum = 0;
        for (j = 0; j < dim; j++)
            sum = pow((data[j] - cent[i][j]), 2) + sum;
        if (min > sum){
            ret = i + 1;
            min = sum;
        }
        i++;
    }
    return ret;
}


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

    /* calculates kmeans like in ex1*/
    
    while (max_iter > count){
        i = 0;
        while(i < len){
            clust = assign_cluster(cent, points[i], dim, in1);
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
        if (checkTheNorm(cent, cent_that_used, dim, in1))
            break;
        idkun_the_cents(cent, cent_that_used, dim, in1);
        count++;
    }

    free_memory(cent_that_used, in1);
    return 0;
}


void idkun_the_cents(double **centNew, double **cent_that_used, int dim, int K){
    int i, j;
    i = 0;
    while (i < K){
        j = 0;
        while (j < dim){
            cent_that_used[i][j] = centNew[i][j];
            j++;
        }
        i++;
    }
}
