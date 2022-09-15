double **cal_spk(double **lnorm, int len, int *K){ 
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

void free_sort(double** mat){
    free(mat[0]);
    free(mat);
}

double **sort_matrix_values(double **mat, int dim){
    int i, j, max_of_index;
    double the_max_value;
    double **sort_mat = alloc_for_mat(dim + 1, dim);
    if (sort_mat == NULL){
        free_memory(mat, dim + 1);
        return NULL;
    }

    i  = 0;
    while(i < dim){
        max_of_index = -1, the_max_value = -1;
        j = 0;
        while(j < dim){
            if (the_max_value < mat[0][j]){
                max_of_index = j;
                the_max_value = mat[0][j];
            }
            j++;
        }
        sort_mat[0][i] = the_max_value;
        sort_mat[i + 1] = mat[max_of_index + 1];
        mat[0][max_of_index] = -1;
        i++;
    }
   
    free_sort(mat);
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
