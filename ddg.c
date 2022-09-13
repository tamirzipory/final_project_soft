double **diagonal_matrix(double **adj_mat, int N){
    int i, j;
    double sum;
    double **diag_mat = matrix_allocation(N, N);
    if (diag_mat == NULL)
        return NULL;
    
    sum=0;
    for (i = 0; i < N; i++){
        sum = 0;
        for (j = 0; j < N; j++){
            sum += adj_mat[i][j];
            diag_mat[i][j] = 0;
        }
        diag_mat[i][i] = sum;
    }
    return diag_mat;
}