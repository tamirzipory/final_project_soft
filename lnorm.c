double **laplacian_matrix(double **diag_mat, double **adj_mat, int N){
    /* mul1 = diag_mat^(-0.5)*adj_mat^(-0.5)
     * mul2 = mul1*diag_mat^(-0.5) = diag_mat^(-0.5)*adj_mat*diag_mat^(-0.5)
     * lnorm=I - mul2 = I - diag_mat^(-0.5)*adj_mat*diag_mat^(-0.5) */
    double **mul1, **mul2, **lnorm;

    cal_D12(diag_mat, N);
    mul1 = calc_mul(N, diag_mat, adj_mat);
    if (mul1 == NULL)
        return NULL;

    mul2 = calc_mul(N, mul1, diag_mat);
    free_memory(mul1, N);
    if (mul2 == NULL)
        return NULL;

    lnorm = calc_id_mat(N);
    if (lnorm == NULL){
        free_memory(mul2, N);
        return NULL;
    }
    sab_matrix(N, lnorm, mul2);
    free_memory(mul2, N);
    return lnorm;
}

/* Receives D- diagonal degree matrix, N- number of rows/columns
 * Updates D to D^(-0.5)*/
void cal_D12(double **diag_mat, int N){
    int i;
    for (i = 0; i < N; i++){
        diag_mat[i][i] = 1 / sqrt(diag_mat[i][i]);
    }
}

/* Receives matrices: A,B and N- number of rows/columns
 * Returns matrix C= A*B
 * If an error occurred returns NULL*/
double **calc_mul(int N, double **A, double **B){
    int i, j, k;

    double **C = matrix_allocation(N, N);
    if (C == NULL)
        return NULL;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];   
        }
    }
    return C;
}

/* Receives matrices: A,B and N- number of rows/columns
 * Updates A= A-B */
void sab_matrix(int N, double **A, double **B){
    int i, j;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++)
            A[i][j] -= B[i][j];
    }
}

/*calculate id-mat*/
double **calc_id_mat(int N){
    int i, j;

    double **id_mat = matrix_allocation(N, N);
    if (id_mat == NULL)
        return NULL;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++)
            id_mat[i][j] = (i == j) ? 1 : 0;
    }
    return id_mat;
}