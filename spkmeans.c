#include "spkmeans.h"
#include "spk.c"
#include "wam.c"
#include "ddg.c"
#include "lnorm.c"
#include "jacobi.c"

void invalid_input(){
    printf("Invalid Input!\n");
    exit(1);
}

void err_print(){
    printf("An Error Has Occurred\n");
    exit(1);
}

/*alloc mat according the dim*/
double **alloc_mat(int rows, int cols){
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
/*return the dimentins according the purpose (dim or num of vectors*/
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

/* receive the mats copy_mat,original_mat and dimentions of them and copy the value of copy_mat of the values of original_mat */
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

/* set the inputs of the file to the mat of data according the dimentions*/
void set_input(FILE *ifp, double **data, int rows, int cols){
    int i, j;
    double value;
    i = 0, j = 0;
    for (i = 0; i < rows; i++){
        j = 0;
        while(j < cols){
            if(1 != fscanf(ifp, "%lf", &value))
                  j--;
            else data[i][j] = value;
            fgetc(ifp);
            j++;
        }
    }
    rewind(ifp);
}

/* recive mat and free it*/
void free_memory(double **mat, int rows){
    int i = 0;
    while(i < rows){
        free(mat[i]);
        i++;
    }
    free(mat);
}

/*Replace the */
void msg_and_exit(int type_of_err, int err){
    if (err == 1){
        if (type_of_err != 0)
            err_print();
        else invalid_input();
    }
}

/* print matrix according the pourpose that we got from the user*/
void print_result(double **mat, int rows, int cols, enum Goal target)
{
    int i, j;
    if (target == 4)
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

/*the running function */
double **run_goal(enum Goal target, double **data, int n1, int n2, int *n3){
    double **ret;
    double **mat_dd,**mat_wam, **mat_lnorm;
    if (target == 4)
        return calc_jacob(n1, data);
    
    ret = mat_adj(data, n2, n1);
    if (ret == NULL)
        return NULL;
    if(target == 1)
        return ret;
    mat_wam = ret;
    ret = diag_mat(mat_wam, n1);
    if (target == 2 ||ret == NULL){
        free_memory(mat_wam, n1);
        return ret;
    }
    mat_dd = ret;
    ret = calc_L_mat(mat_dd, mat_wam, n1);
    free_memory(mat_wam, n1);
    free_memory(mat_dd, n1);

    if (target == 3 || ret == NULL)
        return ret;
    mat_lnorm = ret;
    ret = calc_spk_method(mat_lnorm, n1, n3);
    free_memory(mat_lnorm, n1);
    return ret;
}

/*The main function */
int main(int argc, char *argv[]){
    double **data, **ret;
    int n1, n2, n3;
    FILE *ifp;
    enum Goal target = 0;
    n3 = 0;
    msg_and_exit(0, argc != 3);

    if (strcmp("wam", argv[1]) == 0)
        target = code_of_wam;
    else if (strcmp("ddg", argv[1]) == 0)
        target = code_of_ddg;
    else if (strcmp("lnorm", argv[1]) == 0)
        target = code_of_lnorm;
    else if (strcmp("jacobi", argv[1]) == 0)
        target = code_of_jacobi;
    msg_and_exit(0, 0 == target);

    ifp = fopen(argv[2], "r");
    msg_and_exit(1, ifp == NULL);
    n1 = get_n_d_parameters(ifp, 1);
    n2 = get_n_d_parameters(ifp, 2);
    data = alloc_mat(n1, n2);
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
    if (code_of_jacobi != target)
        free_memory(ret, n1);
    else free_memory(ret, n1 + 1);
    fclose(ifp);
    exit(0);
}
