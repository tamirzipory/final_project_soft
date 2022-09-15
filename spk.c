double **calc_spk_method(double **ln, int len, int *K){ 
    double **output_j, **vectors_e, **ret, **sort_t;
    output_j = (double **)calc_jacob(len, ln);
    if (output_j == NULL)
        return NULL;
    vectors_e = output_j + 1; 
    get_mat_transe(vectors_e, len);
    sort_t = sort_matrix_values(output_j, len);
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

double **sort_matrix_values(double **mat, int N){
    int i, j, max_index;
    double max_value;
    double **sort_mat = alloc_mat(N + 1, N);
    if (sort_mat == NULL){
        free_memory(mat, N + 1);
        return NULL;
    }
    for (i = 0; i < N; i++){
        max_index = -1;
        max_value = -1;
        for (j = 0; j < N; j++){
            if (max_value < mat[0][j]){
                max_index = j;
                max_value = mat[0][j];
            }
        }
        sort_mat[0][i] = max_value;
        sort_mat[i + 1] = mat[max_index + 1];
        mat[0][max_index] = -1;
    }
    free(mat[0]);
    free(mat);
    return sort_mat;
}

double **calc_the_T(double **U, int N, int K){
    int i, j, q;
    double sum = 0;

    double **T = alloc_mat(N, K);
    if (T == NULL)
        return NULL;

    for (i = 0; i < N; i++){
        for (j = 0; j < K; j++){
            if (j == 0){
                sum = 0;
                for (q = 0; q < K; q++)
                    sum += pow(U[i][q], 2);
            }
            T[i][j] = (sum != 0) ? (U[i][j] / sqrt(sum)) : 0;
        }
    }
    return T;
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
