double **calc_jacob(int len_mat, double **A){
    int count, i, j, ret;
    double d1, d2, first_A_mat, second_A_mat; 
    double **V_mat, **the_p_mat, **mat_of_jacobi, **finish_V_mat;
    double *values; 
    V_mat = calc_id_mat(len_mat);
    if (V_mat == NULL)
        return NULL;

    the_p_mat = alloc_for_mat(len_mat, len_mat);
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
        calc_first_mat(len_mat, A, d1, d2, i, j, &ret); 
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
    get_values_from_first_mat(values, len_mat, A);           
    mat_of_jacobi = miun_of_eig(len_mat, values, V_mat); 
    free_memory(V_mat, len_mat);
    free_memory(the_p_mat, len_mat);
    free(values);
    return mat_of_jacobi;
}


void calc_first_mat(int len_mat, double **the_mat, double d1, double d2, int i, int j, int *ret){
    int r;
    double **mat_2 = alloc_for_mat(2, len_mat);
    if (NULL == mat_2){
        *ret = 0;
        return;
    }
    r = 0;
    while(r < len_mat){
        if(r != i && r != j){
            mat_2[0][r] = d1 * the_mat[r][i] - d2 * the_mat[r][j];
            mat_2[1][r] = d1 * the_mat[r][j] + d2 * the_mat[r][i];
        }
        else if(r == i || r == j){
            if(r == i)
               mat_2[0][r] = d2 * d2 * the_mat[j][j] + d1 *d1 * the_mat[i][i] - 2 * d2 * d1 * the_mat[i][j];
            if(r == j)
               mat_2[1][r] = 2 * d2 * d1 * the_mat[i][j] + d2 * d2 * the_mat[i][i] + d1 * d1 * the_mat[j][j];
        }
        else {
            mat_2[0][r] = 0;
            mat_2[1][r] = 0;
        }
        r++;
    }
    r = 0;
    while(r < len_mat){
       if (r != i && r != j){
            the_mat[r][i] = mat_2[0][r]; 
            the_mat[i][r] = the_mat[r][i];           
            the_mat[j][r] = mat_2[1][r]; 
            the_mat[r][j] = the_mat[j][r];           
        }
        r++; 
    }
    the_mat[i][j] = 0, the_mat[j][i] = 0, the_mat[i][i] = mat_2[0][i], the_mat[j][j] = mat_2[1][j];
    free_memory(mat_2,2);
    *ret = 1;
}


double calc_off_diag(int len_mat, double **mat_of_the_calc_diag){
    int i, j;
    double squar = 0;
    i = 0;
    while(i < len_mat){
        j = 0;
        while(j < len_mat){
            if(i != j)
                squar = (mat_of_the_calc_diag[i][j] * mat_of_the_calc_diag[i][j]) + squar;
            j++;
        }
        i++;
    }
    return squar;
}

void get_mat_transe(double **mat, int max_iter){
    int i, j;
    double t;
    i = 0;
    while(i < max_iter){
        j = i +1;
        while(j < max_iter){
            t = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = t;
            j++;
        }
        i++;
    }

}

void A_to_A_tag(int len, double **mat, int *index_of_i, int *index_of_j){
    int in1, in2;
    double maximum = -DBL_MAX;
    in1 = 0;
    while(in1 < len){
        for (in2 = in1 + 1; in2 < len; ++in2){
            if (maximum < fabs(mat[in1][in2])){
                maximum = fabs(mat[in1][in2]);
                *index_of_i = in1;
                *index_of_j = in2;
            }
        }
        in1++;
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

void get_values_from_first_mat(double *values, int len, double **mat){
    int i = 0;
    while( i < len){
        values[i] = mat[i][i];
        i++;
    }
}

/*It's do merge of the eighens*/
double **miun_of_eig(int len, double *values, double **vectors){
    double **ret = NULL;
    int i;
    int plus_one = len+1;
    ret = alloc_for_mat(plus_one, len);
    if (ret == NULL)
        return NULL;
    i = 0;
    while(i <len){
        ret[0][i] = values[i];
        i++;
    }
    matrix_copy(len, len, &ret[1], vectors);
    return ret;
}
