double **set_T(double **the_u_mat, int len, int K){
    int i, j, in;
    double sum = 0;
    double **the_T_mat = alloc_for_mat(len, K);
    if (NULL == the_T_mat)
         return NULL;
    i = 0;
    while (i < len){
        j = 0;
        while (j < K){
             if (j == 0){
                sum = 0;
                for (in = 0; in < K; in++)
                    sum = sum + (the_u_mat[i][in] * the_u_mat[i][in]);
            }
            if(sum != 0)
                the_T_mat[i][j] = the_u_mat[i][j] / sqrt(sum);
            else the_T_mat[i][j] = 0;
        j++;
        }
        i++;
    }
    return the_T_mat;
}

void heru_eigen(double *values, int len, int *K){ 
    int i;
    double curr_max_gap = DBL_MIN;
    i =0;
    while((int)(len / 2) >= i){
         if (fabs(values[i-1]-values[i]) > curr_max_gap){
            curr_max_gap = fabs(values[i-1]-values[i]);
            *K = i;
        }
        i++;
    }
}
double **cal_spk(double **lNormat, int len, int *K){ 
    double **patt_jacobi;
    double **the_e_vectors;
    double **the_T_mat;
    double **miun_T_mat;
    patt_jacobi = (double **)calc_jacob(len, lNormat);
    if (NULL == patt_jacobi)
        return NULL;
    the_e_vectors = 1 + patt_jacobi;
    get_mat_transe(the_e_vectors, len);  
    miun_T_mat = sortMatValues(patt_jacobi, len);
    if(miun_T_mat == NULL || *K == 0){
        if (miun_T_mat == NULL)
             return NULL;
        if (*K == 0) heru_eigen(miun_T_mat[0], len, K);
    }
    the_e_vectors = 1 + miun_T_mat; 
    get_mat_transe(the_e_vectors, len);
    the_T_mat = set_T(the_e_vectors, len, *K);
    free_memory(miun_T_mat, 1 + len);
    return the_T_mat;
}

void free_sort(double** mat){
    free(mat[0]);
    free(mat);
}

double **sortMatValues(double **mat, int dim){
    int i, j, max_of_index;
    double the_max_value, temp;
    double **sort_mat = alloc_for_mat(dim + 1, dim);
    if (NULL == sort_mat){
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
        temp = the_max_value;
        sort_mat[0][i] = temp;
        sort_mat[i + 1] = mat[max_of_index + 1];
        mat[0][max_of_index] = -1;
        i++;
    }
    free_sort(mat);
    return sort_mat;
}
