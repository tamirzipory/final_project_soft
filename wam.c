double **adjacency_matrix(double **data_points, int dimension, int N){
    int i, j;

    double **adj_mat = alloc_for_mat(N, N);
    if (adj_mat == NULL)
        return NULL;

    for (i = 0; i < N; i++){
        for (j = i; j < N; j++){
            adj_mat[i][j] = (i == j) ? 0 : (exp((euc_norm_calc(data_points[i], data_points[j], dimension)) / (-2)));
            adj_mat[j][i] = adj_mat[i][j];
        }
    }
    return adj_mat;
}


double euc_norm_calc(double *arr1, double *arr2, int dimension){
    int j;
    double sum;
    sum = 0;
    for (j = 0; j < dimension; j++)
        sum += pow(arr1[j] - arr2[j], 2);
    
    sum = sqrt(sum);
    return sum;
}
