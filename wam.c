double **adjacency_matrix(double **data_points, int dimension, int N){
    int i, j;

    double **adj_mat = matrix_allocation(N, N);
    if (adj_mat == NULL)
        return NULL;

    for (i = 0; i < N; i++){
        for (j = i; j < N; j++){
            adj_mat[i][j] = (i == j) ? 0 : (exp((calc_euclidean_norm(data_points[i], data_points[j], dimension)) / (-2)));
            adj_mat[j][i] = adj_mat[i][j];
        }
    }
    return adj_mat;
}

/* Receives 2 vectors- x,y and thier dimension
 * Returns their distance: ||x-y||2*/
double calc_euclidean_norm(double *x, double *y, int dimension){
    int j;
    double sum;
    sum = 0;
    for (j = 0; j < dimension; j++)
        sum += pow(x[j] - y[j], 2);
    
    sum = sqrt(sum);
    return sum;
}
