double **diagonal_matrix(double **mat, int len){
    int i, j;
    double sum;
    double **diag = mat_alloc_by_row_col(len, len);
    if (diag == NULL)
        return NULL;
    sum = 0, i = 0;
    while(i < len){
        sum = 0;
        for (j = 0; j < len; j++){
            sum = mat[i][j] + sum;
            diag[i][j] = 0;
        }
        diag[i][i] = sum;
        i++;
    }
    return diag;
}
