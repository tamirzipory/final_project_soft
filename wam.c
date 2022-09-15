double calc_euclidean_norm(double *arr1, double *arr2, int dim){
    double sum, ret;   
    int j;
    sum = 0;
    j = 0;
    while(j < dim){
        
        sum = ((arr1[j]-arr2[j]) * (arr1[j]-arr2[j])) + sum;
        /*sum = pow(arr1[j] - arr2[j], 2) + sum;*/
        j++;
    }
    ret = sqrt(sum);
    return ret;
}

double **adjacency_matrix(double **points, int dim, int len){
    int i, j;
    double **a_mat = alloc_for_mat(len, len);
    if (a_mat == NULL)
        return NULL;
    i = 0;
    while ( i< len) {
        j = 0;
        while (j < len)
        {
            if(i == j)
                a_mat[i][j] = 0;
            else a_mat[i][j] = (exp((calc_euclidean_norm(points[i], points[j], dim)) / (-2)));
            j++;
        }
        i++;
    }
    return a_mat;
}

