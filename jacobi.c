double **calc_jacob(int len_mat, double **A){
    int counter, i, j, ret;
    double d1, d2, first_A_mat, second_A_mat; 
    double **V_mat, **the_p_mat, **mat_of_jacobi, **finish_V_mat;
    double *eigenvalues; 
    
    /* Memory allocations- if matrix_allocation returns NULL-
     * an error occurred - free previous allocations and return NULL */
    V_mat = calc_id_mat(len_mat); /* In first iteration, V = the_p_mat1*/
    if (V_mat == NULL)
        return NULL;

    the_p_mat = matrix_allocation(len_mat, len_mat);
    if (the_p_mat == NULL){
        free_memory(V_mat, len_mat);
        return NULL;
    }

    eigenvalues = calloc(len_mat, sizeof(double)); /*len of diagonal of squared matrix (NxN) is always N*/
    if (eigenvalues == NULL){
        free_memory(V_mat, len_mat);
        free_memory(the_p_mat, len_mat);
        return NULL;
    }

    counter = 0, ret=1, first_A_mat = 1.00001, second_A_mat = 0;

    /* Will run up to 100 iterations if it doesn't reach convergence before + in the first iteration, convergence is not relevant*/
    while ((100 > counter) && ((first_A_mat - second_A_mat > 0.00001) || (counter == 0))){
        counter++;
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
            free(eigenvalues);
            free_memory(the_p_mat, len_mat);
            return NULL;
        }

        if (V_mat == NULL){ 
            free(eigenvalues);
            free_memory(the_p_mat, len_mat);
            return NULL;
        }

        second_A_mat = calc_off_diag(len_mat, A);
    }

    get_eigenvalues_from_A1(eigenvalues, len_mat, A);           
    mat_of_jacobi = jacobi_eigen_merge(len_mat, eigenvalues, V_mat); 

    free_memory(V_mat, len_mat);
    free_memory(the_p_mat, len_mat);
    free(eigenvalues);

    return mat_of_jacobi;
}


void calc_A1(int len_mat, double **A, double c, double s, int i, int j, int *ret){
    int r;
  
    double **row_col_i_j = matrix_allocation(2, len_mat);

    if (row_col_i_j == NULL){
        *ret = 0;
        return;
    }

    for (r = 0; r < len_mat; r++){
       
        row_col_i_j[0][r] = (r != i && r != j) ? (c * A[r][i] - s * A[r][j]) : ((r == i) ? (c * c * A[i][i] + s * s * A[j][j] - 2 * s * c * A[i][j]) : 0);
        row_col_i_j[1][r] = (r != i && r != j) ? (c * A[r][j] + s * A[r][i]) : ((r == j) ? (s * s * A[i][i] + c * c * A[j][j] + 2 * s * c * A[i][j]) : 0);
    }
    /* updates*/
    for (r = 0; r < len_mat; r++){
        if (r != j && r != i){
            A[r][i] = row_col_i_j[0][r]; 
            A[i][r] = A[r][i];           
            A[j][r] = row_col_i_j[1][r]; 
            A[r][j] = A[j][r];           
        }
    }

    A[i][j] = 0, A[j][i] = 0, A[i][i] = row_col_i_j[0][i], A[j][j] = row_col_i_j[1][j];

    free_memory(row_col_i_j,2);
    *ret = 1;
}


double calc_off_diag(int len_mat, double **A){
    int i, j;
    double off_A_squared = 0;

    for (i = 0; i < len_mat; i++){
        for (j = 0; j < len_mat; j++)
    
            if(i != j){
                off_A_squared = off_A_squared+(A[i][j] * A[i][j]);
            }
            
        
    }
    return off_A_squared;
}

 /* Updates mat to be get_mat_transed */
void get_mat_transe(double **mat, int N){
    int i, j;
    double temp;
    for (i = 0; i < N; i++){
        for (j = i + 1; j < N; j++){
            temp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = temp;
        }
    }
}

/*calculate A tag according the file in the moudle*/
void A_to_A_tag(int N, double **A, int *index_of_i, int *index_of_j){
    int q, l;
    double maximum = -DBL_MAX;
    for (q = 0; q < N; ++q){
        for (l = q + 1; l < N; ++l){
            
            if (fabs(A[q][l]) > maximum){
                maximum = fabs(A[q][l]);
                *index_of_i = q;
                *index_of_j = l;
            }
        }
    }
}

double calc_theta(double** mat, int i, int j){
    return (mat[j][j] - mat[i][i]) / (2 * mat[i][j]);
}

double calc_t(int sign, double theta){
    return (sign) / (fabs(theta) + sqrt((theta * theta) + 1));
}

double divide(double t){
    return (1) / (sqrt(1 + (t * t)));
}

double calc_s(double t){
    return t / sqrt(1 + (t * t));
}



/* Receives matrix A, i,j- the location of the pivot
 * Calculates c,s,t according to the given formulas */
void get_params(double **A, int i, int j, double *point_1, double *point_2){
    double theta, t;
    double signTheta = 1;

    if (A[i][j] == 0){
        *point_1 = 1;
        *point_2 = 0;
        return;
    }

    theta = calc_theta(A, i, j);
    if (theta < 0)
        signTheta = -1;
    t = calc_t(signTheta, theta);
    *point_1 = divide(t);
    *point_2 = t / sqrt(1 + (t * t));
}

/* Receives matrix the_p_mat, N- number of rows/columns, i,j- the location of the pivot, and c,s
 * Updates P according to the given instructions */
void calc_curr_P(int max_iter, double **the_p_mat, int i, int j, double c, double s){
    int k, l;
    /* P[i][i]=P[j][j]=c, P[i][j]=s, P[j][i]=-s 
     * on diagonal= 1, else=0*/
  

    for (k = 0; k < max_iter; ++k){
        for (l = 0; l < max_iter; ++l){
            if (k == l){
             if (k == i){
                 the_p_mat[k][l] = c;
            }
            if(l == j){
                the_p_mat[k][l] = c;
            }
            else{
                the_p_mat[k][l] = 1;
            }
        }

            else if (k == j && l == i)
                the_p_mat[k][l] = -s;
            else{
                if(k == i && l == j){
                    the_p_mat[k][l] = s;
                }
                else{
                    the_p_mat[k][l] = 0;
                }
            }
                
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
