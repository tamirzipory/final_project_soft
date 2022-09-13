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

#include "wam.c"
/* ================================== Done WAM ==================================*/

/* ================================== DDG (Diagonal Degree Matrix) ================================== */

#include "ddg.c"
/* ================================== Done DDG ==================================*/

/* ================================== LNORM (Normalized Graph Laplacian) ================================== */

#include "lnorm.c"
/* ================================== Done LNORM ==================================*/

/* ================================== JACOBI ================================== */

#include "jacobi.c"
/* ================================== Done JACOBI ==================================*/

/* ================================== General/ Main's Functions ==================================*/
double **matrix_allocation(int rows, int cols){
    /*alloc of mat in size of rows * cols*/
    int i;
    double **mat = calloc(rows, (sizeof(double *)));
    if (NULL == mat)
        err_print();
    i = 0;
    while (i < rows){
        mat[i] = calloc(cols, (sizeof(double)) );
        if (NULL == mat[i])
            err_print();
        i++;
    }
    return mat;
}

/*the program recive file and purpose (purpose is calculate dim or num of vectors) and calculate the purpose*/
int get_n_d_parameters(FILE *ifp, int situattion){
    char ch;
    int count = 0;
    ch = 0;
    while ((ch = fgetc(ifp)) != EOF){
        if(situattion != 1){
            if ('\n' == ch){
                rewind(ifp);
                return (1 + count);
            }
            else{
                if (',' == ch)
                    count = count + 1;
                }
        }
        else{
            if ('\n' == ch)
                count = count + 1;
        }
    }
    rewind(ifp);
    return count;
}

/* Receives matrices copy_mat,original_mat and number of rows and columns
 * Updates copy_mat to be a copy of original_mat */
void matrix_copy(int rows, int cols, double **copy_mat, double **original_mat){
    int i, j;
    for (i = 0; i < rows; i++){
        j = 0;
        while(j < cols){
            copy_mat[i][j] = original_mat[i][j];
            j++;
        }        
    }
}

/* Receives file's pointer, pointer to an empty matrix of datapoints and num of rows, num of cols
 * Updated data according to the file's data */
void set_input(FILE *ifp, double **data, int rows, int cols){
    int i, j;
    double value;
    i = 0, j = 0;
    for (i = 0; i < rows; i++){
        j = 0;
        while(j < cols){
            if(1 != fscanf(ifp, "%lf", &value))
                  j--;
            else data[i][j] = value;;
            fgetc(ifp);
            j++;
        }
    }
    rewind(ifp);
}

/* Receives a matrix: mat_to_free and rows- matrix's number of rows
 * Frees mat_to_free*/
void free_memory(double **mat, int rows){
    int i = 0;
    while(i < rows){
        free(mat[i]);
    }
    free(mat);
}

/* Receives an int: error_type (2 options: invalid (INVALID_TYPE) or error (ERROR_TYPE)) and is_error- boolean number
 * If is_error is 1 (true- an error occurred), print the correspond message (according to error_type) and exit
 * Else- continue (do nothing) */
void msg_and_exit(int type_of_err, int err){
    if (err == 1){
        if (type_of_err != 0)
            err_print();
        else invalid_input();
    }
}

/* Receives a matrix, number of rows and columns, and enum Goal
 * Prints the matrix (if the goal is jacobi then updates rows +1) */
void print_result(double **mat, int rows, int cols, enum Goal target)
{
    int i, j;
    if (target == JACOBI)
        rows++;
    i = 0;
    while(i < rows){
        j = 0;
        while(j < cols){
            if (cols - 1 == j)
                printf("%.4f", mat[i][j]);
            else printf("%.4f,", mat[i][j]);
            j++;
        }
        if ((rows - 1) != i)
            printf("\n");
        i++;
    }
}

/* Receives an enum Goal, matrix with given data (from file), N- number of rows (or number of points),
 * D- number of columns (or point's dimension) and a pointer to K
 * Returns corresponed matrix according to the given goal
 * If an error occured returns NULL*/
double **run_goal(enum Goal goal, double **data, int N, int D, int *K){
    double **data_output, **wam_matrix, **ddg_matrix, **lnorm_matrix;

    if (goal == 4){
        data_output = calc_jacob(N, data);
        return data_output;
    }

    /* Run WAM*/
    data_output = adjacency_matrix(data, D, N);
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
    double **data, **data_output;
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
    data = matrix_allocation(N, D);
    msg_and_exit(1, data == NULL);

    /* Sets the N points/symmetric matrix in data*/
    set_input(ifp, data, N, D);

    /* Sets the goal's result in data_output*/
    data_output = run_goal(goal, data, N, D, &K);
    if (data_output == NULL){ /* An error has occurred*/
        free_memory(data, N);
        msg_and_exit(1, 1);
    }

    print_result(data_output, N, N, goal);
    printf("\n");
    free_memory(data, N);
    if (goal == jacobi_g)
        free_memory(data_output, N + 1);
    else
        free_memory(data_output, N);
    fclose(ifp);
    exit(0);
}
