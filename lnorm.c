double **calc_L_mat(double **d_mat, double **a_mat, int len){
   
    double **m1, **m2, **ret;
    calc_norm_mat(d_mat, len);
    m1 = calc_mul(len, d_mat, a_mat);
    if (m1 == NULL)
        return NULL;
    m2 = calc_mul(len, m1, d_mat);
    free_memory(m1, len);
    if (m2 == NULL)
        return NULL;
    ret = calc_id_mat(len);
    if (ret == NULL){
        free_memory(m2, len);
        return NULL;
    }
    sab_matrix(len, ret, m2);
    free_memory(m2, len);
    return ret;
}


double divide_lnorm(double** mat, int i, int j){
    return (1 / sqrt(mat[i][j]));
}

void calc_norm_mat(double **d_mat, int len){
    int i;
    for (i = 0; i < len; i++){
        d_mat[i][i] = divide_lnorm(d_mat, i, i);
        
    }
}

/*do mul between 2 mats and put the result in new ret mat*/
double **calc_mul(int dim_of_the_mats, double **mat1, double **mat2){
    int i, j, index2;
    double temp;
    double **ret = alloc_mat(dim_of_the_mats, dim_of_the_mats);
    if (ret == NULL)
        return NULL;

    i = 0;
    while(dim_of_the_mats>i){
        j = 0;
        while(dim_of_the_mats>j){
            ret[i][j] = 0;
            index2 = 0;
            while(dim_of_the_mats > index2){
                temp = mat1[i][index2] * mat2[index2][j];
                ret[i][j] += temp;
                index2++;
            }
            j++;
        }
        i++;
    }


    return ret;
}

/*sum mat1 and mat2 and put the result in mat1*/
void sab_matrix(int dim, double **mat1, double **mat2){
    int i, j;
    double temp;
    i = 0;
    while(i < dim){
        j = 0;
        while(j < dim){
            temp = mat2[i][j];
            mat1[i][j] = mat1[i][j]-temp;
            j++;
        }
        i++;
    }
  
}

/*calculate id-mat*/
double **calc_id_mat(int dim_of_mat){
    int i, j;
    double **id_mat = alloc_mat(dim_of_mat, dim_of_mat);
    if (id_mat == NULL)
        return NULL;
    i = 0, j = 0;
    while(i < dim_of_mat){
        j = 0;
        while(dim_of_mat > j){
            if(i == j){
                id_mat[i][j] = 1;
            }
            else{
                id_mat[i][j] = 0;
            }
            j = j + 1;;
        }
        i = i+1;
    }
    return id_mat;
}
