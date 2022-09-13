double **calc_jacob(int len_mat, double **A){
    int count, i, j, ret;
    double d1, d2, first_A_mat, second_A_mat; 
    double **V_mat, **the_p_mat, **mat_of_jacobi, **finish_V_mat;
    double *values; 
    V_mat = calc_id_mat(len_mat);
    if (V_mat == NULL)
        return NULL;

    the_p_mat = matrix_allocation(len_mat, len_mat);
    if (the_p_mat == NULL){
        free_memory(V_mat, len_mat);
        return NULL;
    }

    values = calloc(len_mat, sizeof(double)); 
    if (values == NULL){
        free_memory(V_mat, len_mat);
        free_memory(the_p_mat, len_mat);
        return NULL;
    }
    count = 1, ret=1, first_A_mat = 1.00001, second_A_mat = 0;
    for(;((first_A_mat - second_A_mat > 0.00001) || (count == 0)) && (100 > count); count++){
        first_A_mat = calc_off_diag(len_mat, A);
        A_to_A_tag(len_mat, A, &i, &j);                   
        if (A[i][j] == 0)
            break;
        get_params(A, i, j, &d1, &d2);            
        calc_curr_P(len_mat, the_p_mat, i, j, d1, d2);      
        calc_A1(len_mat, A, d1, d2, i, j, &ret); 
        finish_V_mat=V_mat, V_mat = calc_mul(len_mat, V_mat, the_p_mat); 
        free_memory(finish_V_mat,len_mat);
        if(ret == 0){
            free(values);
            free_memory(the_p_mat, len_mat);
            return NULL;
        }
        if (V_mat == NULL){ 
            free(values);
            free_memory(the_p_mat, len_mat);
            return NULL;
        }
        second_A_mat = calc_off_diag(len_mat, A);
    }
    get_eigenvalues_from_A1(values, len_mat, A);           
    mat_of_jacobi = jacobi_eigen_merge(len_mat, values, V_mat); 
    free_memory(V_mat, len_mat);
    free_memory(the_p_mat, len_mat);
    free(values);
    return mat_of_jacobi;
}


void calc_A1(int len_mat, double **A, double d1, double d2, int i, int j, int *ret){
    int r;
    double **rows_cols = matrix_allocation(2, len_mat);
    if (NULL == rows_cols){
        *ret = 0;
        return;
    }
    r = 0;
    while(r < len_mat){
        if(r != i && r != j){
            rows_cols[0][r] = d1 * A[r][i] - d2 * A[r][j];
            rows_cols[1][r] = d1 * A[r][j] + d2 * A[r][i];
        }
        else if(r == i || r == j){
            if(r == i)
               rows_cols[0][r] = d2 * d2 * A[j][j] + d1 *d1 * A[i][i] - 2 * d2 * d1 * A[i][j];
            if(r == j)
               rows_cols[1][r] = 2 * d2 * d1 * A[i][j] + d2 * d2 * A[i][i] + d1 * d1 * A[j][j];
        }
        else {
            rows_cols[0][r] = 0;
            rows_cols[1][r] = 0;
        }
        r++;
    }
    r = 0;
    while(r < len_mat){
       if (r != i && r != j){
            A[r][i] = rows_cols[0][r]; 
            A[i][r] = A[r][i];           
            A[j][r] = rows_cols[1][r]; 
            A[r][j] = A[j][r];           
        }
        r++; 
    }
    A[i][j] = 0, A[j][i] = 0, A[i][i] = rows_cols[0][i], A[j][j] = rows_cols[1][j];
    free_memory(rows_cols,2);
    *ret = 1;
}


double calc_off_diag(int len_mat, double **A){
    int i, j;
    double squar = 0;
    i = 0;
    while(i < len_mat){
        j = 0;
        while(j < len_mat){
            if(i != j)
                squar = (A[i][j] * A[i][j]) + squar;
            j++;
        }
        i++;
    }
    return squar;
}

void get_mat_transe(double **mat, int N){
    int i, j;
    double t;
    i = 0;
    while(i < N){
        j = i +1;
        while(j < N){
            t = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = t;
            j++;
        }
        i++;
    }
}

/*calculate A tag according the file in the moudle*/
void A_to_A_tag(int N, double **A, int *i, int *j){
    int in1, in2;
    double max = -DBL_MAX;
    in1 = 0;
    while(in1 < N){
         for (in2 = in1 + 1; in2 < N; ++in2){
            if (fabs(A[in1][in2]) > max){
                max = fabs(A[in1][in2]);
                *i = in1;
                *j = in2;
            }
        }
        in1 = in1 +1;
    }
}

double calc_theta(double** mat, int i, int j){
    double ret = (mat[j][j] - mat[i][i]) / (mat[i][j] * 2);
    return ret;
}

double calc_t(int s, double t){
    double ret = s / (sqrt((t * t) + 1 + fabs(t)));
    return ret;
}

double divide(double t){
    double ret = 1 / (sqrt((t * t) + 1));
    return ret;
}

double calc_s(double t){
    return t / sqrt(1 + (t * t));
}

void get_params(double **A, int i, int j, double *p1, double *p2){
    double th, t;
    double sign = 1;
    if (A[i][j] == 0){
        *p1 = 1;
        *p2 = 0;
        return;
    }
    th = calc_theta(A, i, j);
    if (0 > th)
        sign = -1;
    t = calc_t(sign, th);
    *p1 = divide(t);
    *p2 = t / sqrt((t * t) + 1);
}
/*
void calc_curr_P(int max_iter, double **the_p_mat, int i, int j, double d1, double d2){
    int in1, in2;
    in1 = 0;
    while(in1 < max_iter){
        in2 = 0;
        while(in2 < max_iter){
            if(in1 == in2){
                if(in1 == i || in2 == j)
                   the_p_mat[in1][in2] = d1;
                else the_p_mat[in1][in2] = 1;   
            }
            else if (in1 == j && in2 == i)
                the_p_mat[in1][in2] = -d2;
            else{
                if(in1 == i && in2 == j)
                    the_p_mat[in1][in2] = d2;
                else the_p_mat[in1][in2] = 0;
            }
            in2++;
        }
        in1++;
    }
}
*/

void calc_curr_P(int max_iter, double **the_p_mat, int i, int j, double c, double s){
    int k, l;
    for (k = 0; k < max_iter; ++k){
        for (l = 0; l < max_iter; ++l){
            if (k == l)
                the_p_mat[k][l] = (k == i || l == j) ? c : 1; /* 1 on the diagonal*/
            else if (k == j && l == i)
                the_p_mat[k][l] = -s;
            else
                the_p_mat[k][l] = (k == i && l == j) ? s : 0;
        }
    }
}

/* Receives matrix A1, N- number of rows/columns and a pointer to eigenvalues's array
 * Updates eigenvalues according to the diagonal of A1 */
void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1){
    int i;
    for (i = 0; i < N; ++i)
        eigenvalues[i] = A1[i][i];
}

/* Receives matrix eigenVectors, N- number of rows/columns and an array of the eigenvalues
 * Returns matrix with first row- eigenValues, next rows- eigenVectors */
double **jacobi_eigen_merge(int len, double *eigenValues, double **eigenVectors){
    double **res = NULL;
    int i;
    int plus_one = len+1;
    res = matrix_allocation(plus_one, len);
    if (res == NULL)
        return NULL;
    for (i = 0; i < len; i++)
        res[0][i] = eigenValues[i];
    matrix_copy(len, len, &res[1], eigenVectors);
    return res;
}
