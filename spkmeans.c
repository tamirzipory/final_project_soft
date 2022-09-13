#include "spkmeans.h"

/*main*/
int main(int argc, char *argv[])
{
    char *input_file_name;
    GOAL goal = invalid;
    matrix *data_points = NULL;
    matrix *out = NULL;
    eigen_data *j_out = NULL;

    /*validate the input*/
    argc = argc - 1;
    if (argc != 2)
    {
        return invalid_input();
    }

    goal = string_to_goal(argv[1]);
    input_file_name = argv[2];
    data_points = read_input_file(input_file_name);
    if (data_points == NULL)
    {
        return error();
    }

    switch (goal)
    {
    case wam:
        out = calculate_w(data_points);
        break;
    case ddg:
        out = calculate_d(data_points);
        break;
    case lnorm:
        out = calculate_l(data_points);
        break;
    case jacobi:
        /*check that data_points is a square matrix*/
        if (data_points->rows != data_points->columns)
        {
            free_matrix(data_points);
            return invalid_input();
        }
        j_out = calculate_jacobi(data_points);
        break;
    default: /*goal wasn't valid*/
        free_matrix(data_points);
        return invalid_input();
    }

    if (goal != jacobi)
    {
        if (out == NULL) /*an error occurred in matrix calculation*/
        {
            free_matrix(data_points);
            return error();
        }
        else
        {
            print_matrix(out);
            free_matrix(out);
        }
    }
    else
    {
        if (j_out == NULL) /*an error occurred in jacobi algorithm calculation*/
        {
            free_matrix(data_points);
            return error();
        }
        else
        {
            /*prints the eigenvalues without the "-" before zeros*/
            print_matrix_formatted(j_out->eigenvalues, 1);
            print_matrix(j_out->eigenvectors_mat);
            free_eigen_data(j_out);
        }
    }

    free_matrix(data_points);
    return 0;
}

/*implementation - general*/
/*converts from string to enum GOAL.
from https://stackoverflow.com/questions/16844728/converting-from-string-to-enum-in-c*/
GOAL string_to_goal(const char *str)
{
    int j;
    for (j = 0; j < NUMBER_OF_VALID_ENUM_OPTIONS; j++)
    {
        if (!strcmp(str, conversion[j].str))
        {
            return conversion[j].val;
        }
    }
    return invalid;
}

/*implementation - functions of kmeans++*/
/*  gets a csv/txt file name (assumes that the file is in the right format and that it exists),
 *  returns a matrix that represents the parsed data points from the file.
 *  each row in the matrix is a point and each cell is a coordinate.
 *  @param input_file_name the name of the input file.
 *  @return struct matrix containing the data points. returns NULL on error.
 */
matrix *read_input_file(char *input_file_name)
{
    FILE *ifp = NULL;
    matrix *data_points = NULL;
    char c;
    double coordinate = 0;
    int i = 0, j = 0, d = 1, n = 0;

    /*read input file*/
    ifp = fopen(input_file_name, "r");

    /*cannot open file*/
    if (ifp == NULL)
    {
        return NULL;
    }

    /*find d*/
    while ((c = fgetc(ifp)) != '\n') /*there must be \n before EOF because file is in the right format*/
    {
        if (c == ',')
            d++;
    }
    rewind(ifp);

    /*find n*/
    while ((c = fgetc(ifp)) != EOF)
    {
        if (c == '\n')
            n++;
    }
    rewind(ifp);

    /*create matrix*/
    data_points = create_matrix(n, d);
    if (data_points == NULL)
    {
        fclose(ifp);
        return NULL;
    }

    /*read points from file*/
    while (fscanf(ifp, "%lf,", &coordinate) > 0)
    {
        if (j == d) /*we finished a point, move to a the next point*/
        {
            j = 0;
            i++;
        }
        data_points->values[i][j] = coordinate;
        j++;
    }

    /*close the file*/
    fclose(ifp);

    return data_points;
}

/*  the main part of the kmeans algorithm.
 *  the function gets the parsed data points and the different parameters,
 *  and returns the final centroids.
 *  @param data_points the parsed data points.
 *  @param centroids the chosen centroids (from python). IMPORTANT - the function changes the given centroids!
 *  @param max_iter the maximal number of iterations.
 *  @param epsilon the precision of the result.
 *  @return a pointer to the given centroids (after the function changed them). returns NULL on error.
 */
matrix *fit_kmeanspp_c(matrix *data_points, matrix *centroids, int max_iter, double epsilon)
{
    int mu_flag = 1;
    int i;
    int iter_counter = 1;
    int closest_centroid_index;
    int d = data_points->columns, n = data_points->rows, k = centroids->rows;
    double *point = NULL;
    double *new_centroid = NULL;
    cluster *clusters = NULL;

    /*allocate current centroid*/
    new_centroid = (double *)calloc(d, sizeof(double));
    if (new_centroid == NULL) /*calloc of new_centroid failed*/
    {
        return NULL;
    }

    /*create clusters*/
    clusters = create_clusters(k, n);
    if (clusters == NULL)
    {
        free(new_centroid);
        return NULL;
    }

    while (mu_flag && iter_counter <= max_iter)
    {
        reset_clusters(clusters, k);
        for (i = 0; i < n; i++)
        {
            /*put points into clusters*/
            point = data_points->values[i];
            closest_centroid_index = find_closest_centroid(point, centroids);
            add_point(&(clusters[closest_centroid_index]), i);
        }

        /*update centroids*/
        mu_flag = 0;
        for (i = 0; i < k; i++)
        {
            if (update_centroid(&(clusters[i]), centroids->values[i], new_centroid, data_points, epsilon) == 1)
            {
                mu_flag = 1;
            }
            memset(new_centroid, 0, d * sizeof(double));
        }
        iter_counter++;
    }

    /*free memory*/
    free_clusters(clusters, k);
    free(new_centroid);

    return centroids;
}

/*  gets the clusters array and resets each cluster.
 *  @param clusters array of clusters.
 *  @param k the number of clusters.
 */
void reset_clusters(cluster *clusters, int k)
{
    int i;
    for (i = 0; i < k; i++)
    {
        clusters[i].num_points = 0;
    }
}

/*  allocates the clusters array and retuns pointer to the array.
 *  @param k the number of clusters.
 *  @param n the number of data points.
 */
cluster *create_clusters(int k, int n)
{
    cluster *clusters = NULL;
    int *points_indexes = NULL;
    int i;

    clusters = (cluster *)calloc(k, sizeof(cluster));
    if (clusters == NULL) /*calloc of clusters failed*/
    {
        return NULL;
    }

    for (i = 0; i < k; i++)
    {
        points_indexes = (int *)calloc(n, sizeof(int));
        if (points_indexes == NULL) /*calloc of points_indexes failed*/
        {
            free_clusters(clusters, k);
            return NULL;
        }

        clusters[i].num_points = 0;
        clusters[i].points_indexes = points_indexes;
    }
    return clusters;
}

/*  gets a data point and the centroids array
 *  and returns the index of the closest centroid (in euclidian norm).
 *  @param point a data point.
 *  @param centroids the current centroids.
 */
int find_closest_centroid(double *point, matrix *centroids)
{
    int i;
    int k = centroids->rows, d = centroids->columns;
    int min_i = 0;
    double min_norm, norm;

    min_norm = euclidean_norm_squared(centroids->values[0], point, d); /*k must be > 1*/
    for (i = 1; i < k; i++)
    {
        norm = euclidean_norm_squared(centroids->values[i], point, d);
        if (norm < min_norm)
        {
            min_norm = norm;
            min_i = i;
        }
    }
    return min_i;
}

/*  gets a data point index and adds it as a point in the given cluster.
 *  @param cluster poinetr to the cluster.
 *  @param point_index a point index in the data points matrix.
 */
void add_point(cluster *cluster, int point_index)
{
    int i = cluster->num_points;
    cluster->points_indexes[i] = point_index;
    cluster->num_points++;
}

/*  gets a cluster and the current centroids and updates it according to the points in the cluster.
 *  if the difference between the old centroid and the new centroid is smaller then epsilon, 0 will be returened.
 *  else 1 will be returned indicating that we need to keep updating this centroid.
 *  @param current_cluster pointer to the cluster.
 *  @param current_centroid.
 *  @param new_centroid must be initiallized with 0s!!!
 *  @param data_points the data points.
 *  @param epsilon the algorithm precision.
 */
int update_centroid(cluster *current_cluster, double *current_centroid, double *new_centroid, matrix *data_points, double epsilon)
{
    double *current_point = NULL;
    int i, j, current_point_index;
    int d = data_points->columns;
    double dist;

    for (i = 0; i < current_cluster->num_points; i++)
    {
        current_point_index = current_cluster->points_indexes[i];
        current_point = data_points->values[current_point_index];

        for (j = 0; j < d; j++)
        {
            new_centroid[j] += current_point[j] / current_cluster->num_points;
        }
    }

    dist = sqrt(euclidean_norm_squared(new_centroid, current_centroid, d));
    /*update current centroid*/
    memcpy(current_centroid, new_centroid, sizeof(double) * d);
    if (dist < epsilon)
    {
        return 0;
    }
    return 1;
}

/*  frees the clusters array.
 *  @param clusters the clusters array.
 *  @param k the number of clusters.
 */
void free_clusters(cluster *clusters, int k)
{
    int i;
    for (i = 0; i < k; i++)
    {
        free(clusters[i].points_indexes);
    }
    free(clusters);
}

/*print error message and return (for main in c).*/
int invalid_input()
{
    printf("Invalid Input!");
    return 1;
}
int error()
{
    printf("An Error Has Occurred");
    return 1;
}

/*implementation - functions of The Normalized Spectral Clustering Algorithm*/
/*calculates euclidean norm between two points (that are represented by one dimensional array).
the returned result is squared.*/
double euclidean_norm_squared(double *point1, double *point2, int d)
{
    double sum = 0;
    int i;

    for (i = 0; i < d; i++)
    {
        sum += pow(point1[i] - point2[i], 2);
    }
    return sum;
}

/*prints a matrix in the desired format.
if format_zero is 1 then "-0.0000" will be printed as "0.0000".
else, the output will be formatted only with "%.4f".*/
void print_matrix_formatted(matrix *mat, int format_zero)
{
    int i, j;
    double val;

    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->columns - 1; j++)
        {
            val = mat->values[i][j];
            printf((format_zero && fabs(val) < 0.00005) ? "0.0000," : "%.4f,", val);
        }

        val = mat->values[i][mat->columns - 1];
        /*print the last coordinate in the row seperatly with \n*/
        printf((format_zero && fabs(val) < 0.00005) ? "0.0000\n" : "%.4f\n", val);
    }
}

/*prints a matrix in the default format.*/
void print_matrix(matrix *mat)
{
    print_matrix_formatted(mat, 0);
}

/*creates a new matrix with the result of dot(mat1, mat2).
if there is an allocation error - returns NULL.
the dimensions of the matrixes are assumed to be legal.*/
matrix *dot_product(matrix *mat1, matrix *mat2)
{
    /*the number of rows in mat1 is the number of rows in the result matrix.
    the number of columns in mat2 is the number of columns in the result matrix.*/
    int i, j, k;
    matrix *mat = NULL;
    mat = create_matrix(mat1->rows, mat2->columns);
    if (mat == NULL)
    {
        return NULL;
    }
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->columns; j++)
        {
            for (k = 0; k < mat1->columns; k++)
            {
                mat->values[i][j] += mat1->values[i][k] * mat2->values[k][j];
            }
        }
    }
    return mat;
}

/*creates the Weighted Adjancency Matrix from the given data points.
the diagonal is all 0s, and the rest of the cells are calculated according to the formula.
the dimensions of the created matrix are n*n, as n = data_points->rows (the number of points).
in case of an allocation error, returns NULL.*/
matrix *calculate_w(matrix *data_points)
{
    int i, j;
    double norm;
    matrix *w_mat = NULL;
    w_mat = create_matrix(data_points->rows, data_points->rows);
    if (w_mat == NULL)
    {
        return NULL;
    }
    for (i = 0; i < w_mat->rows; i++)
    {
        for (j = i + 1; j < w_mat->columns; j++)
        {
            {
                norm = sqrt(euclidean_norm_squared(data_points->values[i], data_points->values[j], data_points->columns));
                w_mat->values[i][j] = exp(norm * -0.5);
                w_mat->values[j][i] = w_mat->values[i][j];
            }
        }
    }
    return w_mat;
}

/*creates the Diagonal Degree Matrix from a calculted Weighted Adjancency Matrix.
the diagonal is calculated according to the formula, the rest of the cells are 0s.
the dimensions of the created matrix are the same as w matrix.
in case of an allocation error, returns NULL.*/
matrix *calculate_d_from_w(matrix *w_mat)
{
    int i, j;
    matrix *d_mat = NULL;
    /*create d matrix*/
    d_mat = create_matrix(w_mat->rows, w_mat->columns);
    if (d_mat == NULL)
    {
        return NULL;
    }
    for (i = 0; i < d_mat->rows; i++)
    {
        for (j = 0; j < d_mat->columns; j++)
        {
            d_mat->values[i][i] += w_mat->values[i][j];
        }
    }

    return d_mat;
}

/*creates the Diagonal Degree Matrix from the data points.
the diagonal is calculated according to the formula, the rest of the cells are 0s.
the dimensions of the created matrix are the same as w matrix.
in case of an allocation error, returns NULL.*/
matrix *calculate_d(matrix *data_points)
{
    matrix *d_mat = NULL, *w_mat = NULL;
    /*calculate w in order to calculate d*/
    w_mat = calculate_w(data_points);
    if (w_mat == NULL)
    {
        return NULL;
    }
    /*create d matrix*/
    d_mat = calculate_d_from_w(w_mat);
    if (d_mat == NULL)
    {
        free_matrix(w_mat);
        return NULL;
    }

    /*free w_mat*/
    free_matrix(w_mat);

    return d_mat;
}

/*calculates for each i between 0 and mat->rows: 1/sqrt(mat[i][i]),
and sets mat[i][i] to be the calculated value.
the changes are in-place.*/
void sqrt_diagonal(matrix *mat)
{
    int i;
    for (i = 0; i < mat->rows; i++)
    {
        mat->values[i][i] = 1 / sqrt(mat->values[i][i]);
    }
}

/*creates the Normalized Graph Laplacian from a calculted Weighted Adjancency Matrix and Diagonal Degree Matrix.
the matrix is calculated according to the formula.
in case of an allocation error, returns NULL.*/
matrix *calculate_l(matrix *data_points)
{
    int i, j;
    matrix *l_mat = NULL, *d_mat = NULL, *w_mat = NULL;
    matrix *temp = NULL;

    /*calculate w matrix*/
    w_mat = calculate_w(data_points);
    if (w_mat == NULL)
    {
        return NULL;
    }
    /*calculate d matrix*/
    d_mat = calculate_d_from_w(w_mat);
    if (d_mat == NULL)
    {
        free_matrix(w_mat);
        return NULL;
    }
    /*calculate D^(-0.5), changes D in-place!!!*/
    sqrt_diagonal(d_mat);

    /*calculate D^(-0.5) * W * D^(-0.5).*/
    temp = dot_product(d_mat, w_mat);
    if (temp == NULL)
    {
        free_matrix(w_mat);
        free_matrix(d_mat);
        return NULL;
    }
    l_mat = dot_product(temp, d_mat);
    if (l_mat == NULL)
    {
        free_matrix(w_mat);
        free_matrix(d_mat);
        free_matrix(temp);
        return NULL;
    }

    /*free all*/
    free_matrix(w_mat);
    free_matrix(d_mat);
    free_matrix(temp);

    /*calculate I - D^(-0.5) * W * D^(-0.5) (== (-l_mat) + I)*/
    for (i = 0; i < l_mat->rows; i++)
    {
        for (j = 0; j < l_mat->columns; j++)
        {
            l_mat->values[i][j] *= -1;
            if (i == j)
            {
                l_mat->values[i][j] += 1;
            }
        }
    }
    return l_mat;
}

/*normalizes the rows of the given matrix according to the formula.
the change is in-place.*/
void normalize_matrix(matrix *mat)
{
    int i, j;
    double sum = 0;
    for (i = 0; i < mat->rows; i++)
    {
        /*for each row, calculate the sum of squares*/
        for (j = 0; j < mat->columns; j++)
        {
            sum += pow(mat->values[i][j], 2);
        }
        sum = sqrt(sum);
        /*normalize*/
        for (j = 0; j < mat->columns; j++)
        {
            mat->values[i][j] = sum > 0 ? mat->values[i][j] / sum : 0;
        }
        /*init sum for next row*/
        sum = 0;
    }
}

/*creates identity matrix with dimensions of n*n*/
matrix *create_identity_matrix(int n)
{
    int i;
    matrix *i_mat = NULL;
    i_mat = create_matrix(n, n);
    if (i_mat == NULL)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        i_mat->values[i][i] = 1;
    }
    return i_mat;
}

/*creates the eigenvalues array and the eigenvectors matrix V.
the input a, is a valid symmetric matrix.
in case of an allocation error, returns NULL.*/
eigen_data *calculate_jacobi(matrix *a)
{
    /*i and j are indexes for the pivot in the jacobi algorithm.
    we preform the algorithm only on matrixes bigger then 1X1 so it is safe to assume that j=1 exists.
    the initialization is 0,1 because the defult pivot is 0,1 entry of the matrix*/
    int i = 0, j = 1;

    int l = 0, count_iter = 0;
    double off_a, off_a_tag;
    matrix *v = NULL, *p = NULL, *temp = NULL;
    eigen_data *result = NULL;

    /*create identity matrix to initialize V matrix*/
    v = create_identity_matrix(a->rows);
    if (v == NULL)
    {
        return NULL;
    }

    /*check if the matrix is allready diagonal.
    if the matrix a is diagonal - don't preform the jacobi algorithm.
    in that case v will be the identity matrix and the eigenvalues will be on the diagonal of a.*/
    if (off(a) != 0)
    {
        /*preforming the algorithm*/
        do
        {
            count_iter++;
            /*find maximal off-diagonal a[i][j] in a matrix*/
            find_pivot(a, &i, &j);
            /*calculate off before rotation*/
            off_a = off(a);
            /*calculate rotation matrix p*/
            p = calculate_p(a, i, j);
            if (p == NULL)
            {
                free_matrix(v);
                return NULL;
            }
            /*rotate a, this changes a (in-place) to be a' = p^t*a*p */
            rotate_matrix(a, p, i, j);
            /*calculate off after rotation*/
            off_a_tag = off(a);
            /*calculte v = v*p */
            temp = dot_product(v, p);
            if (temp == NULL)
            {
                free_matrix(v);
                free_matrix(p);
                return NULL;
            }
            free_matrix(v); /*free the old v*/
            v = temp;       /*update the current v to be the dot product*/
            free_matrix(p); /*free the current rotation matrix*/

        } while (off_a - off_a_tag > CONVERGENCE_EPSILON && count_iter < MAX_ITER_JACOBI);
    }

    /*prepering the result*/
    /*allocate struct eigen_data for the returned value*/
    result = (eigen_data *)calloc(1, sizeof(eigen_data));
    if (result == NULL) /*calloc of result failed*/
    {
        free_matrix(v);
        return NULL;
    }
    /*allocate eigenvalues one dimensional matrix*/
    result->eigenvalues = create_matrix(1, a->rows);
    if (result->eigenvalues == NULL)
    {
        free(result); /*only struct exists at this point*/
        free_matrix(v);
        return NULL;
    }
    /*make the eigenvalues array from the diagonal of final a after all rotations*/
    for (l = 0; l < a->rows; l++)
    {
        result->eigenvalues->values[0][l] = a->values[l][l];
    }

    result->eigenvectors_mat = v; /*v is the eigenvectors of a*/

    return result;
}

/*creates the rotation matrix P according to the formula.
the dimensions of the created matrix are the same as a.
i and j are row and column of the pivot element.
in case of an allocation error, returns NULL.*/
matrix *calculate_p(matrix *a, int i, int j)
{
    double theta, t, c, s;
    int sign_theta;

    matrix *p = NULL;
    p = create_identity_matrix(a->rows);
    if (p == NULL)
    {
        return NULL;
    }

    /*calculate according to the formula*/
    theta = (a->values[j][j] - a->values[i][i]) / (2 * a->values[i][j]);
    sign_theta = (theta < 0) ? -1 : 1;
    t = sign_theta / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    c = 1 / sqrt(pow(t, 2) + 1);
    s = t * c;
    /*put the values in their place in p*/
    p->values[i][i] = c;
    p->values[j][j] = c;
    p->values[i][j] = s;
    p->values[j][i] = -1 * s;
    return p;
}

/*finds the off-diagonal element with the largset absolute value in a (a must be symmetric matrix).
i will be the index of it's row and j of it's column*/
void find_pivot(matrix *a, int *i, int *j)
{
    int row, col;
    double max = 0, curr;

    /*initialize i,j to be defult values*/
    *i = 0;
    *j = 1;

    for (row = 0; row < a->rows; row++)
    {
        for (col = row + 1; col < a->columns; col++)
        {
            curr = fabs(a->values[row][col]);
            if (curr > max)
            {
                max = curr;
                *i = row;
                *j = col;
            }
        }
    }
}

/*creates the rotated matrix A' according to the formula.
the change is in-place.
i and j are row and column of the pivot element.*/
void rotate_matrix(matrix *a, matrix *p, int i, int j)
{
    double c = p->values[i][i];
    double s = p->values[i][j];
    double temp_ri, temp_ii, temp_jj;
    int r = 0;

    for (r = 0; r < a->rows; r++)
    {
        if (r != i && r != j)
        {
            temp_ri = a->values[r][i];
            a->values[r][i] = c * temp_ri - s * a->values[r][j];
            a->values[i][r] = a->values[r][i];
            a->values[r][j] = c * a->values[r][j] + s * temp_ri;
            a->values[j][r] = a->values[r][j];
        }
    }

    temp_ii = a->values[i][i];
    temp_jj = a->values[j][j];
    a->values[i][i] = pow(c, 2) * temp_ii + pow(s, 2) * temp_jj - 2 * s * c * a->values[i][j];
    a->values[j][j] = pow(s, 2) * temp_ii + pow(c, 2) * temp_jj + 2 * s * c * a->values[i][j];
    a->values[i][j] = 0;
    a->values[j][i] = 0;
}

/*calculats off(a)^2 according to the formula.*/
double off(matrix *a)
{
    double sum = 0;
    int i, j;
    for (i = 0; i < a->rows; i++)
    {
        for (j = 0; j < a->columns; j++)
        {
            if (i != j)
            {
                sum += pow(a->values[i][j], 2);
            }
        }
    }
    return sum;
}

/*compare function for qsort of two eigan_pair pointers. sort is according to eigenvalue*/
int compare_eigen_pair(const void *p1, const void *p2)
{
    const eigen_pair *ep1 = p1, *ep2 = p2;
    double comp;

    comp = (ep1->eigenvalue) - (ep2->eigenvalue);
    /*negative if ev1 < ev2. 0 if equal. positive if ev1 > ev2*/
    if (comp < 0)
    {
        return -1;
    }
    if (comp > 0)
    {
        return 1;
    }
    /*equal*/
    return 0;
}

/*creates the pairs of eigenvalue-eigenvector from the eigen data.
the result is sorted according to the eigenvalues, ascending order*/
eigen_pair *extract_sorted_eigen_pairs(eigen_data *jacobi_data)
{
    eigen_pair *eigen_pairs = NULL;
    int i;
    int row, col;
    int length = jacobi_data->eigenvalues->columns;
    double *eigenvalues_jacobi_arr = jacobi_data->eigenvalues->values[0];
    matrix *eigenvectors_jacobi_mat = jacobi_data->eigenvectors_mat;

    /*allocate the result array and check allocation*/
    eigen_pairs = (eigen_pair *)calloc(length, sizeof(eigen_pair));
    if (eigen_pairs == NULL) /*calloc of eigen_pairs failed*/
    {
        return NULL;
    }

    /*allocate eigenvectors in all pairs and set eigenvalues*/
    for (i = 0; i < length; i++)
    {
        /*allocation:*/
        /*the lentgh of an eigenvector is the same as the number of eigenvalues,
        becuase the number of eigenvalues is fit to the dimensions of the symmetric matrix.*/
        eigen_pairs[i].eigenvector = create_matrix(1, length);
        if (eigen_pairs[i].eigenvector == NULL)
        {
            free_eigen_pairs(eigen_pairs, length);
            return NULL;
        }

        /*set eigenvalues:
        the eigenvalue in index i of eigen_pairs represents
        the eigenvalue connected to eigenvector in column i*/
        eigen_pairs[i].eigenvalue = eigenvalues_jacobi_arr[i];
    }

    /*extract the eigenvectors one cooridinate at a time*/
    for (row = 0; row < length; row++)
    {
        for (col = 0; col < length; col++)
        {
            /*assign to each vector the relevant coordinate*/
            eigen_pairs[col].eigenvector->values[0][row] = eigenvectors_jacobi_mat->values[row][col];
        }
    }

    /*now eigen_pairs is the array of the connected eigenvalue-eigenvector.
    sort it in-place*/
    qsort(eigen_pairs, length, sizeof(eigen_pair), compare_eigen_pair);

    return eigen_pairs;
}

/*recieves an array of eigen_pair structs sorted(!!!) according to the eigenvalues, ascending order.
finds k according to the formula.
n is the number of eigenvalues, which is the length of eigen_pairs*/
int find_k_heuristic(eigen_pair *eigen_pairs, int n)
{
    double max_delta = 0, curr_delta;
    int k = 1;
    int i;
    for (i = 0; i < n / 2; i++)
    {
        curr_delta = fabs(eigen_pairs[i].eigenvalue - eigen_pairs[i + 1].eigenvalue);
        if (curr_delta > max_delta)
        {
            max_delta = curr_delta;
            k = i + 1;
        }
    }
    return k;
}

/*creates the T Matrix from sorted eigen pairs.
T is the normalized matrix of first k eigenvectors, according to the algorithm.
the dimensions of the created matrix are (n = number of eigen pairs, k = chosen number of clusters)
in case of an allocation error, returns NULL.*/
matrix *calculate_t(eigen_pair *eigen_pairs, int n, int k)
{
    matrix *u = NULL;
    int row, col;

    /*creates u - the matrix of k first eigenvectors (not normalized)*/
    u = create_matrix(n, k);
    if (u == NULL)
    {
        return NULL;
    }
    for (row = 0; row < n; row++)
    {
        for (col = 0; col < k; col++)
        {
            /*assign the relevant coordinate of the eigenvector to the matrix*/
            u->values[row][col] = eigen_pairs[col].eigenvector->values[0][row];
        }
    }

    /*normalize u*/
    normalize_matrix(u);

    return u;
}

/*make the full process of the algorithm except K-means algorithn itself.
if k=0 use hueristic to find k,
else use the given k (this k was given as input from the user and is valid).
returns T - the normalized matrix of first k eigenvectors.
in case of an allocation error, returns NULL.*/
matrix *fit_nsc_c(matrix *data_points, int k)
{
    matrix *t = NULL, *lnorm = NULL;
    eigen_data *ed_jacobi = NULL;
    eigen_pair *eigen_pairs = NULL;
    int n = data_points->rows; /*number of data points*/

    /*jacobi algorithm - calculate eigenvalues and eigenvectors from data points*/
    lnorm = calculate_l(data_points);
    if (lnorm == NULL)
    {
        return NULL;
    }

    ed_jacobi = calculate_jacobi(lnorm);
    if (ed_jacobi == NULL)
    {
        free_matrix(lnorm);
        return NULL;
    }

    /*extract the sorted eigen pairs*/
    eigen_pairs = extract_sorted_eigen_pairs(ed_jacobi);
    if (eigen_pairs == NULL)
    {
        free_matrix(lnorm);
        free_eigen_data(ed_jacobi);
        return NULL;
    }

    /*find k, if k=0*/
    if (k == 0)
    {
        k = find_k_heuristic(eigen_pairs, n);
        /*k=1 is an error in the algorithm according to the instructions*/
        if (k == 1)
        {
            free_matrix(lnorm);
            free_eigen_data(ed_jacobi);
            free_eigen_pairs(eigen_pairs, n);
            return NULL;
        }
    }

    /*calculate t*/
    t = calculate_t(eigen_pairs, n, k);
    if (t == NULL)
    {
        free_matrix(lnorm);
        free_eigen_data(ed_jacobi);
        free_eigen_pairs(eigen_pairs, n);
        return NULL;
    }

    /*free memory*/
    free_matrix(lnorm);
    free_eigen_data(ed_jacobi);
    free_eigen_pairs(eigen_pairs, n);

    return t;
}

/*allocates an empty matrix filled with 0s.
returns a pointer to the matrix.
if the allocation fails, returns NULL*/
matrix *create_matrix(int rows, int columns)
{
    int i = 0;
    matrix *mat = NULL;
    /*allocate struct matrix and check allocation*/
    mat = (matrix *)calloc(1, sizeof(matrix));
    if (mat == NULL) /*calloc of mat failed*/
    {
        return NULL;
    }
    /*set values of matrix*/
    mat->rows = rows;
    mat->columns = columns;
    /*allocate values of matrix and check allocation*/
    mat->values = NULL;
    mat->values = (double **)calloc(rows, sizeof(double *));
    if (mat->values == NULL) /*calloc of mat->values failed*/
    {
        free(mat);
        return NULL;
    }
    for (i = 0; i < rows; i++)
    {
        mat->values[i] = (double *)calloc(columns, sizeof(double));
        if (mat->values[i] == NULL) /*calloc of mat->values[i] failed*/
        {
            free_matrix(mat);
            return NULL;
        }
    }
    return mat;
}

/*free an allocated matrix.*/
void free_matrix(matrix *mat)
{
    int i;
    if (mat != NULL)
    {
        if (mat->values != NULL)
        {
            for (i = 0; i < mat->rows; i++)
            {
                free(mat->values[i]);
            }
            free(mat->values);
        }
        free(mat);
    }
}

/*free an allocated struct of eigen_data*/
void free_eigen_data(eigen_data *jacobi_data)
{
    free_matrix(jacobi_data->eigenvalues);
    free_matrix(jacobi_data->eigenvectors_mat);
    free(jacobi_data);
}

/*free an allocated array of struct eigen_pair*/
void free_eigen_pairs(eigen_pair *eigen_pairs, int length)
{
    int i;
    for (i = 0; i < length; i++)
    {
        free_matrix(eigen_pairs[i].eigenvector);
    }
    free(eigen_pairs);
}
