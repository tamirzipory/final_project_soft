double **diag_mat(double **mat, int len){
    int i, j;
    double sum;
    double **diag = matrix_allocation(len, len);
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
