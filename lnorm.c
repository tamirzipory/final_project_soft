double **calc_L_mat(double **d_mat, double **a_mat, int len){
    double **m1;
    double **m2;
    double **ret;
    calcDMat(d_mat, len);
    m1 = calc_mul(len, d_mat, a_mat);
    if (NULL == m1)
        return NULL;
    m2 = calc_mul(len, m1, d_mat);
    if (NULL == m1)
        return NULL;
    free_memory(m1, len);
    if (NULL == m2)
        return NULL;
    ret = calc_id_mat(len);
    if (NULL == m2)
        return NULL;
    if (NULL == ret){
        free_memory(m1, len);
        return NULL;
    }
    sab_matrix(len, ret, m2);
    free_memory(m2, len);
    return ret;
}

double divide_lnorm(double** mat, int i, int j){
    return (1 / sqrt(mat[i][j]));
}


void calcDMat(double **diag_mat, int len){
    int i = 0;
    while (i < len)
    {
        diag_mat[i][i] = divide_lnorm(diag_mat, i, i);
        i++;
    }
}

/*create new mat with valus of the mult of mat1*mat2 when the mult is like algebra linearit*/
double **calc_mul(int dim_of_the_mats, double **mat1, double **mat2){
    int i, j, z;
    double temp;
    double **ret = alloc_for_mat(dim_of_the_mats, dim_of_the_mats);
    if (ret == NULL)
        return NULL;

    i = 0;
    while(dim_of_the_mats>i){
        j = 0;
        while(dim_of_the_mats>j){
            ret[i][j] = 0;
            z = 0;
            while(dim_of_the_mats > z){
                temp = mat1[i][z] * mat2[z][j];
                ret[i][j] += temp;
                z++;
            }
            j++;
        }
        i++;
    }
    return ret;
}

/*sub mat1-mat2 and put the result in mat1*/
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
    double **id_mat = alloc_for_mat(dim_of_mat, dim_of_mat);
    if (NULL == id_mat)
        return NULL;
    i = 0, j = 0;
    while(i < dim_of_mat){
        j = 0;
        while(dim_of_mat > j){
            if(i == j)
                id_mat[i][j] = 1;
            else id_mat[i][j] = 0;
            j = j + 1;;
        }
        i = i+1;
    }
    return id_mat;
}
