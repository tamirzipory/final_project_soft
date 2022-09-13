#include "spkmeans.h"


void invalid_input(){
    printf("Invalid Input!\n");
    exit(1);
}

void err_print(){
    printf("An Error Has Occurred\n");
    exit(1);
}


#include "spk.c"

/* ================================== WAM (Weighted Adjacency Matrix) ================================== */
/* Receives N datapoints and their dimension
 * Returns the corresponding weighted adjacency matrix.
 * If an error occurred returns NULL*/
#include "wam.c"
/* ================================== Done WAM ==================================*/

/* ================================== DDG (Diagonal Degree Matrix) ================================== */
/* Receives weighted adjacency matrix and N- number of rows/columns
 * Returns the corresponding diagonal degree matrix.
 * If an error occurred returns NULL*/
#include "ddg.c"
/* ================================== Done DDG ==================================*/

/* ================================== LNORM (Normalized Graph Laplacian) ================================== */
/* Receives diagonal degree matrix,weighted adjacency matrix and N- number of rows/columns
 * Returns the corresponding normalized graph Laplacian matrix.
 * If an error occurred returns NULL*/
#include "lnorm.c"
/* ================================== Done LNORM ==================================*/

/* ================================== JACOBI ================================== */
/* Receives symetric matrix A- size N*N
 * Returns matrix with first row= eigenvalues, next rows are the
 * corresponding eigenvectors (each column is a vector)
 * If an error occurred returns NULL*/
#include "jacobi.c"
/* ================================== Done JACOBI ==================================*/

/* ================================== General/ Main's Functions ==================================*/
double **matrix_allocation(int num_rows, int num_cols){
    /*allocation of memory of size (nxn) and return 1 if there was a failure!*/
    int i;
    double **mat = malloc((sizeof(double *)) * num_rows);
    if (mat == NULL)
        return NULL;
    for (i = 0; i < num_rows; i++){
        mat[i] = malloc((sizeof(double)) * (num_cols));
        if (mat[i] == NULL)
            return NULL;
    }
    return mat;
}

/* Receives file's pointer and an int- find_who (2 options: N (FIND_N) or D (FIND_D))
 * Returns N or D according to find_who */
int get_n_d_parameters(FILE *ifp, int find_who){
    int count;
    char c;
    count = 0;
    while ((c = fgetc(ifp)) != EOF){
        /* N- calculated by number of rows ("\n")*/
        if (find_who == FIND_N){
            if (c == '\n')
                count++;
        }
        else{
            /* D- calculated by number of comas (+ 1) in first row.
             * After done reading first row- returns calculated D */
            if (c == '\n'){
                rewind(ifp);
                count++;
                return count;
            }
            else{
                if (c == ',')
                    count++;
            }
        }
    }
    rewind(ifp);
    return count;
}

/* Receives matrices dest_mat,src_mat and number of rows and columns
 * Updates dest_mat to be a copy of src_mat */
void matrix_copy(int num_rows, int num_cols, double **dest_mat, double **src_mat){
    int i, j;
    for (i = 0; i < num_rows; i++){
        for (j = 0; j < num_cols; j++)
            dest_mat[i][j] = src_mat[i][j];
        
    }
}

/* Receives file's pointer, pointer to an empty matrix of datapoints and num of rows, num of cols
 * Updated data_input according to the file's data */
void set_input(FILE *ifp, double **data_input, int num_rows, int num_cols){
    int i, j;
    double curr_value;
    i = 0, j = 0;

    for (i = 0; i < num_rows; i++){
        for (j = 0; j < num_cols; j++){
            if (fscanf(ifp, "%lf", &curr_value) == 1)
                data_input[i][j] = curr_value;
            else
                j--;
            fgetc(ifp);
        }
    }
    rewind(ifp);
}

/* Receives a matrix: mat_to_free and num_rows- matrix's number of rows
 * Frees mat_to_free*/
void free_memory(double **mat_to_free, int num_rows){
    int i;
    for (i = 0; i < num_rows; i++)
        free(mat_to_free[i]);
    free(mat_to_free);
}

/* Receives an int: error_type (2 options: invalid (INVALID_TYPE) or error (ERROR_TYPE)) and is_error- boolean number
 * If is_error is 1 (true- an error occurred), print the correspond message (according to error_type) and exit
 * Else- continue (do nothing) */
void msg_and_exit(int error_type, int is_error){
    if (is_error == 1){
        if (error_type == 0)
           invalid_input();
        else
            err_print();
        
    }
}

/* Receives a matrix, number of rows and columns, and enum Goal
 * Prints the matrix (if the goal is jacobi then updates num_rows +1) */
void print_result(double **mat, int num_rows, int num_cols, enum Goal goal)
{
    int i, j;
    if (goal == JACOBI)
        num_rows = num_rows+1;
    for (i = 0; i < num_rows; i++){
        for (j = 0; j < num_cols; j++){
            if (j == num_cols - 1)
                printf("%.4f", mat[i][j]);
            else
                printf("%.4f,", mat[i][j]);
        }
        if (i != num_rows - 1)
            printf("\n");
    }
}

/* Receives an enum Goal, matrix with given data (from file), N- number of rows (or number of points),
 * D- number of columns (or point's dimension) and a pointer to K
 * Returns corresponed matrix according to the given goal
 * If an error occured returns NULL*/
double **run_goal(enum Goal goal, double **data_input, int N, int D, int *K){
    double **data_output, **wam_matrix, **ddg_matrix, **lnorm_matrix;

    if (goal == 4){
        data_output = calc_jacob(N, data_input);
        return data_output;
    }

    /* Run WAM*/
    data_output = adjacency_matrix(data_input, D, N);
    if (goal == 1 || data_output == NULL)
        return data_output;

    wam_matrix = data_output;

    /* Run DDG*/
    data_output = diagonal_matrix(wam_matrix, N);
    if (data_output == NULL || goal == 2){
        free_memory(wam_matrix, N);
        return data_output;
    }

    ddg_matrix = data_output;

    /* Run LNORM*/
    data_output = laplacian_matrix(ddg_matrix, wam_matrix, N);
    free_memory(wam_matrix, N);
    free_memory(ddg_matrix, N);
    if (data_output == NULL || goal == 3)
        return data_output;

    lnorm_matrix = data_output;

    /* run SPK*/
    data_output = spk_algo(lnorm_matrix, N, K);
    free_memory(lnorm_matrix, N);
    return data_output;
}

/* Receives k, goal and file_name from user
 * Calculates needed information from file and call run_goal, prints result at end*/
int main(int argc, char *argv[]){
    char *file_name;
    int N, D, K;
    double **data_input, **data_output;
    FILE *ifp;
    enum Goal goal = 0;
    K = 0;

    /* invalid number of arguments*/
    msg_and_exit(0, argc != 3);

    /* set goal correct enum*/
    if (!strcmp("wam", argv[1]))
        goal = wam_g;
    if (!strcmp("ddg", argv[1]))
        goal = ddg_g;
    if (!strcmp("lnorm", argv[1]))
        goal = lnorm_g;
    if (!strcmp("jacobi", argv[1]))
        goal = jacobi_g;
    msg_and_exit(0, goal == 0);

    file_name = argv[2];
    ifp = fopen(file_name, "r");
    msg_and_exit(1, ifp == NULL);
    N = get_n_d_parameters(ifp, 1);
    D = get_n_d_parameters(ifp, 2);

    /* Creates matrix for input*/
    data_input = matrix_allocation(N, D);
    msg_and_exit(1, data_input == NULL);

    /* Sets the N points/symmetric matrix in data_input*/
    set_input(ifp, data_input, N, D);

    /* Sets the goal's result in data_output*/
    data_output = run_goal(goal, data_input, N, D, &K);
    if (data_output == NULL){ /* An error has occurred*/
        free_memory(data_input, N);
        msg_and_exit(1, 1);
    }

    print_result(data_output, N, N, goal);
    printf("\n");
    free_memory(data_input, N);
    if (goal == jacobi_g)
        free_memory(data_output, N + 1);
    else
        free_memory(data_output, N);
    fclose(ifp);
    exit(0);
}
