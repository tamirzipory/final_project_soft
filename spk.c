double **spk_algo(double **lnorm, int N, int *K)
{ /* Called after steps 1-2 have been made*/
    double **jacobi_output, **eigenvectors, **T, **sort_transpose;

    jacobi_output = (double **)jacobi_algo(N, lnorm);
    if (jacobi_output == NULL)
        return NULL;

    /* Transpose on eigenvectors- to make the sort easier*/
    eigenvectors = jacobi_output + 1; /* Jacobi without eigenvalues*/
    transpose(eigenvectors, N);

    /* in sort_transpose jacobi_output is being freed, there is no use in it again!*/
    sort_transpose = sort_matrix_values(jacobi_output, N);
    if (sort_transpose == NULL)
        return NULL;

    /* The Eigengap Heuristic- was told not to handle a case where k=1*/
    if (*K == 0)
        eigengap_heuristic(sort_transpose[0], N, K);

    /* Transpose on eigenvectors- get them to the right shape (vector in columns)*/
    eigenvectors = sort_transpose + 1; /* sort_transpose without eigenvalues*/
    transpose(eigenvectors, N);

    /* eigenvectors points to the start of eigenvectors, we will use only the first K vectors (first K columns) as U
     * and update eigenvectors (by renormalizing each of U’s rows) to be T */
    T = set_T(eigenvectors, N, *K);
    free_memory(sort_transpose, N + 1);

    return T;
}

/* Receives a jacobi's matrix
 * Sort first row and rows 1 to N according to the eigenvalues in first row */
double **sort_matrix_values(double **mat, int N)
{
    int i, j, max_index;
    double max_value;
    double **sort_mat = matrix_allocation(N + 1, N);
    if (sort_mat == NULL){
        free_memory(mat, N + 1);
        return NULL;
    }
    for (i = 0; i < N; i++){
        max_index = -1;
        max_value = -1;
        for (j = 0; j < N; j++){
            /* Found new max*/
            if (max_value < mat[0][j])
            {
                max_index = j;
                max_value = mat[0][j];
            }
        }
        /* Place the i'th eigenvalue in the i'th cell and it's correspoond eigenvectors in line number i+1 */
        sort_mat[0][i] = max_value;
        sort_mat[i + 1] = mat[max_index + 1];
        mat[0][max_index] = -1;
    }
    /* free (mat=jacobi_output) */
    free(mat[0]);
    free(mat);
    return sort_mat;
}

/* Receives U (created by largest eigenvectors of jacobi), N- number of rows, K- number of columns
 * Returns T- by renormalizing each of U’s rows to have unit length */
double **set_T(double **U, int N, int K){
    int i, j, q;
    double sum = 0;

    double **T = matrix_allocation(N, K);
    if (T == NULL)
        return NULL;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            /* Calculate sum once for each new row!*/
            if (j == 0)
            {
                sum = 0;
                for (q = 0; q < K; q++)
                    sum += pow(U[i][q], 2);
            }
            T[i][j] = (sum != 0) ? (U[i][j] / sqrt(sum)) : 0;
        }
    }
    return T;
}

/* Receives eigenvalues (in decreasing order), N- number of values, K- number of clusters
 * Calculate K as described (Eigengap Heuristic algorithm) */
void eigengap_heuristic(double *eigenvalues, int N, int *K)
{ /* (Assumption) lnorm formed as a decreasing ordered eigenvalues*/
    int i;
    double curr_max_gap = DBL_MIN;

    /* lmda(1)= E[0]>=lmda(2)=E[1]>=...>=lmda(n/2)=E[(N/2)-1]>=0*/
    for (i = 1; i <= (int)(N / 2); i++)
    {
        if (curr_max_gap < fabs(eigenvalues[i - 1] - eigenvalues[i])){
            curr_max_gap = fabs(eigenvalues[i - 1] - eigenvalues[i]);
            *K = i;
        }
    }
}