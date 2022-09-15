double **spk_algo(double **lnorm, int len, int *K){ /* Called after steps 1-2 have been made*/
    double **jacobi_pattern, **eigenvectors, **the_T_mat, **miun_T_mat;

    jacobi_pattern = (double **)calc_jacob(len, lnorm);
    if (jacobi_pattern == NULL)
        return NULL;

  
    eigenvectors = jacobi_pattern + 1; 
    get_mat_transe(eigenvectors, len);  
    miun_T_mat = sort_matrix_values(jacobi_pattern, len);
    if (miun_T_mat == NULL)
        return NULL;


    if (*K == 0)
        eigengap_heuristic(miun_T_mat[0], len, K);

    eigenvectors = miun_T_mat + 1; 
    get_mat_transe(eigenvectors, len);
    the_T_mat = set_T(eigenvectors, len, *K);
    free_memory(miun_T_mat, len + 1);
    return the_T_mat;
}


double **sort_matrix_values(double **mat, int N){
    int i, j, max_index;
    double max_value;
    double **sort_mat = alloc_for_mat(N + 1, N);
    if (sort_mat == NULL){
        free_memory(mat, N + 1);
        return NULL;
    }
    for (i = 0; i < N; i++){
        max_index = -1;
        max_value = -1;
        for (j = 0; j < N; j++){
            /* Found new max*/
            if (max_value < mat[0][j]){
                max_index = j;
                max_value = mat[0][j];
            }
        }
       
        sort_mat[0][i] = max_value;
        sort_mat[i + 1] = mat[max_index + 1];
        mat[0][max_index] = -1;
    }
   =
    free(mat[0]);
    free(mat);
    return sort_mat;
}


double **set_T(double **the_u_mat, int len, int K){
    int i, j, in;
    double sum = 0;
    double **the_T_mat = alloc_for_mat(len, K);
    if (the_T_mat == NULL)
        return NULL;

    for (i = 0; i < len; i++){
        for (j = 0; j < K; j++){
           
            if (j == 0){
                sum = 0;
                for (in = 0; in < K; in++)
                    sum += pow(the_u_mat[i][in], 2);
            }
            the_T_mat[i][j] = (sum != 0) ? (the_u_mat[i][j] / sqrt(sum)) : 0;
        }
    }
    return the_T_mat;
}


void eigengap_heuristic(double *eigenvalues, int len, int *K){ 
    int i;
    double curr_max_gap = DBL_MIN;


    for (i = 1; i <= (int)(len / 2); i++){
        if (curr_max_gap < fabs(eigenvalues[i - 1] - eigenvalues[i])){
            curr_max_gap = fabs(eigenvalues[i - 1] - eigenvalues[i]);
            *K = i;
        }
    }
}
