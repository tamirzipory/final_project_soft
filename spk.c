double **calc_spk_method(double **ln, int len, int *K){ 
    double **output_j, **vectors_e, **ret, **sort_t;
    output_j = (double **)calc_jacob(len, ln);
    if (output_j == NULL)
        return NULL;
    vectors_e = output_j + 1; 
    get_mat_transe(vectors_e, len);
    sort_t = mat_sorting(output_j, len);
    if (sort_t == NULL)
        return NULL;
    if (*K == 0)
        func_heruicsic(sort_t[0], len, K);
    vectors_e = sort_t + 1; 
    get_mat_transe(vectors_e, len);
    ret = calc_the_T(vectors_e, len, *K);
    free_memory(sort_t, len + 1);
    return ret;
}

double **mat_sorting(double **mat, int len){
    int i, j, max_index;
    double max_value;
    double **ret = alloc_mat(len + 1, len);
    if (ret == NULL){
        free_memory(mat, len + 1);
        return NULL;
    }
    i = 0;
    while (i < len)
    {
        max_index = -1;
        max_value = -1;
        j=0;
        while(j < len){
            if (max_value < mat[0][j]){
                max_index = j;
                max_value = mat[0][j];
            }
            j++;
        }
        ret[0][i] = max_value;
        ret[i + 1] = mat[max_index + 1];
        mat[0][max_index] = -1;
        i++;
    }
    free(mat[0]);
    free(mat);
    return ret;
}

double **calc_the_T(double **u_mat, int len, int K){
    int i, j, m;
    double sum = 0;
    double **ret = alloc_mat(len, K);
    if (ret == NULL)
        return NULL;
    i = 0;
    while (i < len)
    {
        j = 0;
        while(j < K){
            if (j == 0){
                sum = 0;
                for (m = 0; m < K; m++)
                    sum = pow(u_mat[i][m], 2) + sum;
            }
            if(sum != 0)
                ret[i][j] = (u_mat[i][j] / sqrt(sum));
            else ret[i][j] = 0;
            j++;
        }
        i++;
    }
    return ret;
}


void func_heruicsic(double *eigenvalues, int len, int *K){ 
    int i;
    double curr_max_gap = DBL_MIN;

    int max_iter = (int)(len / 2);
    for (i = 1; i <= max_iter; i++){
        if (curr_max_gap < fabs(eigenvalues[i - 1] - eigenvalues[i])){
            curr_max_gap = fabs(eigenvalues[i - 1] - eigenvalues[i]);
            *K = i;
        }
    }
}
