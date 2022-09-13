#include "spkmeans.h"
#include <string.h>

int rows, cols;

static void print_matrix(double** matrix){
    int i, j;
    if (eigen_values != NULL){
        for (i = 0; i < rows-1; i++){
            if (eigen_values[i].val <= 0 && eigen_values[i].val > -0.0001)
                printf("%.4f%c", 0.0000, ',');
            else
                printf("%.4f%c",eigen_values[i].val, ',');
            }
        if (eigen_values[i].val <= 0 && eigen_values[i].val > -0.0001)
            printf("%.4f%s", 0.0000, "\n");
        else
            printf("%.4f%s", eigen_values[i].val, "\n");
        free(eigen_values);
    }

    for (i = 0; i < rows; i++) {
        for (j = 0; j < rows-1; j++)
            printf("%.4f%c", matrix[i][j], ',');
       printf("%.4f%s", matrix[i][j], "\n");
    }
}

static double** get_matrix(char* file_name){
    int i, j;
    double d;
    double** matrix;
    char c = 0;
    FILE *file_input;
    cols = 0, rows = 1;

    file_input = fopen(file_name, "r");
    if (file_input == NULL)
        return NULL;

    while (c != '\n') {
        fscanf(file_input, "%*f%c", &c);
        cols++;
    }
    while (fscanf(file_input, "%*s\n") != EOF)
        rows++;
    rewind(file_input);
    matrix = Array_2D(rows, cols, 1);
    if (matrix == NULL)
        return NULL;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            fscanf(file_input, "%lf%*c", &d);
            matrix[i][j] = d;
        }
    }
    fclose(file_input);
    return matrix;
}

int main(int argc, char* argv[]) {
    char goal[6];
    double **matrix, **W, ** temp;

    if (argc != 3){
        printf("Invalid Input!");
        return 1;
    }
    sscanf(argv[1], "%s", goal);
    matrix = get_matrix(argv[2]);
    if (matrix == NULL) {
        printf("Invalid Input!");
        return 1;
    }

    temp = matrix;
    if (strcmp(goal, "wam") == 0)
        matrix = wam(matrix, rows, cols);
    else if (strcmp(goal, "ddg") == 0) {
        W = wam(matrix, rows, cols);
        matrix = ddg(W, rows);
        free_Mem(W, rows);
    }
    else if (strcmp(goal, "lnorm") == 0)
        matrix = L_norm(matrix, rows, cols);
    else if (strcmp(goal, "jacobi") == 0){
        matrix = jacobi(matrix, rows);
        }
    else{
        printf("Invalid Input!");
        return 1;
    }

    if (matrix == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    print_matrix(matrix);
    free_Mem(temp, rows);
    free_Mem(matrix, rows);
    return 0;
}
