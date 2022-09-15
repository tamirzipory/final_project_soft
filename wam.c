double **adjacency_matrix(double **mat_of_vectors, int dimension, int len){
    int i, j;

    double **adj_mat = alloc_for_mat(len, len);
    if (adj_mat == NULL)
        return NULL;

    for (i = 0; i < len; i++){
        for (j = i; j < len; j++){
            adj_mat[i][j] = (i == j) ? 0 : (exp((euc_norm_calc(mat_of_vectors[i], mat_of_vectors[j], dimension)) / (-2)));
            adj_mat[j][i] = adj_mat[i][j];
        }
    }
    return adj_mat;
}

double euc_norm_calc(double *arr1, double *arr2, int dimension){
    int j;
    double sum, ret;
    sum = 0;
    for (j = 0; j < dimension; j++)
        sum = sum+((arr1[j]-arr2[j]) * (arr1[j]-arr2[j]));    
    sum = sqrt(sum);
    ret = sum;
    return ret;
}
