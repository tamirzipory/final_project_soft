double **diag_mat(double **mat, int len){
    int i, j;
    double sum;
    double **diag = alloc_for_mat(len, len);
    if (diag == NULL)
        return NULL;
    sum = 0, i = 0;
    while(i < len){
        sum = 0, j = 0;
        while (j < len){
            sum = mat[i][j] + sum;
            diag[i][j] = 0;
            j++;
        }
        diag[i][i] = sum;
        i++;
    }
    return diag;
}
