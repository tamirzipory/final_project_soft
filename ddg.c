double **diag_mat(double **mat, int len){
    int i, j;
    double sum;
    double **ret = alloc_mat(len, len);
    if (ret == NULL)
        return NULL;
    sum = 0, i = 0;
    while(i < len){
        sum = 0;
        j = 0;
        while (j < len)
        {
            sum = mat[i][j] + sum;
            ret[i][j] = 0;
            j++;
        }
        ret[i][i] = sum;
        i++;
    }
    return ret;
}
