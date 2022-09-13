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
        i++;
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
double **run_goal(enum Goal target, double **data, int n1, int n2, int *n3){
    double **ret;
    double **mat_dd,**mat_wam, **mat_lnorm;
    if (target == 4)
        return calc_jacob(n1, data);
    
    ret = adjacency_matrix(data, n2, n1);
    if (ret == NULL)
        return NULL;
    if(target == 1)
        return ret;
    mat_wam = ret;
    ret = diagonal_matrix(mat_wam, n1);
    if (target == 2 ||ret == NULL){
        free_memory(mat_wam, n1);
        return ret;
    }
    mat_dd = ret;
    ret = laplacian_matrix(mat_dd, mat_wam, n1);
    free_memory(mat_wam, n1);
    free_memory(mat_dd, n1);

    if (target == 3 || ret == NULL)
        return ret;
    mat_lnorm = ret;
    ret = spk_algo(mat_lnorm, n1, n3);
    free_memory(mat_lnorm, n1);
    return ret;
}

/* Receives k, goal and file_name from user
 * Calculates needed information from file and call run_goal, prints result at end*/
int main(int argc, char *argv[]){
    double **data, **ret;
    int n1, n2, n3;
    FILE *ifp;
    enum Goal target = 0;
    n3 = 0;
    msg_and_exit(0, argc != 3);

    if (strcmp("wam", argv[1]) == 0)
        target = wam_g;
    else if (strcmp("ddg", argv[1]) == 0)
        target = ddg_g;
    else if (strcmp("lnorm", argv[1]) == 0)
        target = lnorm_g;
    else if (strcmp("jacobi", argv[1]) == 0)
        target = jacobi_g;
    msg_and_exit(0, 0 == target);

    ifp = fopen(argv[2], "r");
    msg_and_exit(1, ifp == NULL);
    n1 = get_n_d_parameters(ifp, 1);
    n2 = get_n_d_parameters(ifp, 2);
    data = matrix_allocation(n1, n2);
    msg_and_exit(1, data == NULL);
    set_input(ifp, data, n1, n2);
    ret = run_goal(target, data, n1, n2, &n3);
    if (NULL == ret){ 
        free_memory(data, n1);
        msg_and_exit(1, 1);
    }

    print_result(ret, n1, n1, target);
    printf("\n");
    free_memory(data, n1);
    if (jacobi_g != target)
        free_memory(ret, n1);
    else free_memory(ret, n1 + 1);
    fclose(ifp);
    exit(0);
}
