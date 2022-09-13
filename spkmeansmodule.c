#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

/*declarations - spkmeansmodule*/
matrix *parse_python_2Dlist_to_matrix(PyObject *unparsed_matrix, int rows, int cols);
PyObject *pack_matrix(matrix *mat);
PyObject *calculate_matrix(PyObject *unparsed_data_points, int n, int d, matrix *(*f)(matrix *));
static PyObject *calculate_matrix_capi(PyObject *self, PyObject *args);
static PyObject *calculate_jacobi_capi(PyObject *self, PyObject *args);
static PyObject *fit_kmeanspp_capi(PyObject *self, PyObject *args);
static PyObject *fit_nsc_capi(PyObject *self, PyObject *args);

/*implementation*/

/*parses a PyObject into struct matrix. returns NULL if fails.*/
matrix *parse_python_2Dlist_to_matrix(PyObject *unparsed_matrix, int rows, int cols)
{
    matrix *mat = NULL;
    PyObject *row = NULL;
    int i, j;

    /*create an empty c struct matrix*/
    mat = create_matrix(rows, cols);
    if (mat == NULL)
    {
        return NULL;
    }
    /*fills the matrix according to the unparsed_matrix from python*/
    for (i = 0; i < rows; i++)
    {
        row = PyList_GetItem(unparsed_matrix, i);
        for (j = 0; j < cols; j++)
        {
            mat->values[i][j] = PyFloat_AsDouble(PyList_GetItem(row, j));
        }
    }

    return mat;
}

/*parses a struct matrix into PyObject. returns NULL if fails*/
PyObject *pack_matrix(matrix *mat)
{
    int i, j;
    PyObject *row = NULL;
    PyObject *packed_matrix = PyList_New((Py_ssize_t)mat->rows);
    if (packed_matrix == NULL)
    {
        return NULL;
    }
    for (i = 0; i < mat->rows; i++)
    {
        row = PyList_New((Py_ssize_t)mat->columns);
        if (row == NULL)
        {
            Py_XDECREF(packed_matrix);
            return NULL;
        }
        for (j = 0; j < mat->columns; j++)
        {
            /*add coordinate to row*/
            PyList_SetItem(row, j, PyFloat_FromDouble(mat->values[i][j]));
        }
        /*add packed row to packed_matrix*/
        PyList_SetItem(packed_matrix, i, row);
    }

    return packed_matrix;
}

/*  gets the unparsed data points from the relevant capi function
 *  and returns the packed result according to the passed function.
 *  @param unparsed_data_points the data points from the capi function.
 *  @param n the number of data points.
 *  @param d the number of coordinates of each point.
 *  @param f the desired function.
 *  @return struct matrix containing the data points. returns NULL on error.
 */
PyObject *calculate_matrix(PyObject *unparsed_data_points, int n, int d, matrix *(*f)(matrix *))
{
    PyObject *final_out_mat = NULL;
    matrix *data_points = NULL, *out_mat = NULL;

    /*parse data points*/
    data_points = parse_python_2Dlist_to_matrix(unparsed_data_points, n, d);
    if (data_points == NULL)
    {
        return NULL;
    }
    /*calculate out matrix*/
    out_mat = (*f)(data_points);
    if (out_mat == NULL)
    {
        free_matrix(data_points);
        return NULL;
    }
    /*pack matrix*/
    final_out_mat = pack_matrix(out_mat);
    if (final_out_mat == NULL)
    {
        free_matrix(data_points);
        free_matrix(out_mat);
        return NULL;
    }
    /*free memory*/
    free_matrix(data_points);
    free_matrix(out_mat);

    return final_out_mat;
}

/*python-c api*/

/* This part actually defines the functions for the python module using a wrapper C API functions.
 * each function has input PyObject *args from Python.
 * The methods parse the PyObject *args and then call the right c function with the parsed parameters.
 * all functions return PyObject that contains the final result.
 */

/*The method receives data points (two dimensional python list) and calculates the desired matrix.
The method returns the matrix represented by two dimensional python list.
recieves the parameters from the python code in the following order:
1. the data points (as a two dimensional list).
2. the number of rows (which is N)
3. the number of columns (which is d)
4. the desired matrix parameter name (given as a string)*/
static PyObject *calculate_matrix_capi(PyObject *self, PyObject *args)
{
    PyObject *unparsed_data_points = NULL;
    char *enum_str = NULL;
    int n, d;
    GOAL goal = invalid;

    if (!PyArg_ParseTuple(args, "Oiis", &unparsed_data_points, &n, &d, &enum_str))
    {
        return NULL;
    }
    /*converts given string to goal, supposed to be valid*/
    goal = string_to_goal(enum_str);

    /*returns the desired matrix*/
    switch (goal)
    {
    case wam:
        return calculate_matrix(unparsed_data_points, n, d, calculate_w);
        break;
    case ddg:
        return calculate_matrix(unparsed_data_points, n, d, calculate_d);
        break;
    case lnorm:
        return calculate_matrix(unparsed_data_points, n, d, calculate_l);
        break;
    default:
        return NULL;
        break;
    }
}

/*The method receives data points (two dimensional python list) and calculates the eigenvalues and the eigenvectors.
The method returns a tuple which contains (eigenvalues matrix, eigenvectors matrix).
recieves the parameters from the python code in the following order:
1. the data points (as a two dimensional list).
2. the number of rows (which is N)
3. the number of columns (which is d)*/
static PyObject *calculate_jacobi_capi(PyObject *self, PyObject *args)
{
    PyObject *unparsed_data_points = NULL, *packed_eigenvalues = NULL, *packed_eigenvectors_mat = NULL, *out_data = NULL;
    matrix *data_points = NULL;
    eigen_data *jacobi_data = NULL;
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &unparsed_data_points, &n, &d))
    {
        return NULL;
    }
    /*parse data points*/
    data_points = parse_python_2Dlist_to_matrix(unparsed_data_points, n, d);
    if (data_points == NULL)
    {
        return NULL;
    }
    /*calculate jacobi data*/
    jacobi_data = calculate_jacobi(data_points);
    if (jacobi_data == NULL)
    {
        free_matrix(data_points);
        return NULL;
    }
    /*pack jacobi data*/
    /*pack eigenvalues and eigenvectors*/
    packed_eigenvalues = pack_matrix(jacobi_data->eigenvalues);
    if (packed_eigenvalues == NULL)
    {
        free_matrix(data_points);
        free_eigen_data(jacobi_data);
        return NULL;
    }
    packed_eigenvectors_mat = pack_matrix(jacobi_data->eigenvectors_mat);
    if (packed_eigenvectors_mat == NULL)
    {
        Py_XDECREF(packed_eigenvalues);
        free_matrix(data_points);
        free_eigen_data(jacobi_data);
        return NULL;
    }
    /*create a tuple containing the above*/
    out_data = PyTuple_New(2);
    if (out_data == NULL)
    {
        Py_XDECREF(packed_eigenvalues);
        Py_XDECREF(packed_eigenvectors_mat);
        free_matrix(data_points);
        free_eigen_data(jacobi_data);
        return NULL;
    }
    PyTuple_SetItem(out_data, 0, packed_eigenvalues);
    PyTuple_SetItem(out_data, 1, packed_eigenvectors_mat);

    /*free memory and return*/
    free_matrix(data_points);
    free_eigen_data(jacobi_data);

    return out_data;
}

/*The method receives data points (two dimensional python list) and initial centroids and
returns the final centroids after applying the kmeans++ algorithm.
recieves the parameters from the python code in the following order:
1. the data points (as a two dimensional list).
2. the initial centroids (as a two dimensional list).
3. n, the number of data points.
4. d, the number of coordinates of each point.
5. k, the number of clusters.
6. max_iter.
7. epsilon.*/
static PyObject *fit_kmeanspp_capi(PyObject *self, PyObject *args)
{
    PyObject *unparsed_data_points = NULL, *unparsed_centroids = NULL, *packed_centroids = NULL;
    matrix *data_points = NULL, *centroids = NULL, *final_centroids = NULL;
    double eps;
    int k, n, d, max_iter;

    if (!PyArg_ParseTuple(args, "OOiiiid", &unparsed_data_points, &unparsed_centroids, &n, &d, &k, &max_iter, &eps))
    {
        return NULL;
    }
    /*parse data*/
    data_points = parse_python_2Dlist_to_matrix(unparsed_data_points, n, d);
    if (data_points == NULL)
    {
        return NULL;
    }
    centroids = parse_python_2Dlist_to_matrix(unparsed_centroids, k, d);
    if (centroids == NULL)
    {
        free_matrix(data_points);
        return NULL;
    }
    /*get the final centroids - final centroids and centroids are pointers to the same memory
    if fit_kmeanspp_c does not return NULL!*/
    final_centroids = fit_kmeanspp_c(data_points, centroids, max_iter, eps);
    if (final_centroids == NULL)
    {
        free_matrix(centroids);
        free_matrix(data_points);
        return NULL;
    }
    /*pack centroids*/
    packed_centroids = pack_matrix(final_centroids);
    if (packed_centroids == NULL)
    {
        free_matrix(centroids);
        free_matrix(data_points);
        return NULL;
    }
    /*free memory - final centroids and centroids are pointers to the same memory!*/
    free_matrix(data_points);
    free_matrix(centroids);

    return packed_centroids;
}

/*The method receives data points (two dimensional python list) and initial centroids and
returns the normilized final data points matrix (t matrix) according to the Normalized Spectral Clustering Algorithm.
recieves the parameters from the python code in the following order:
1. the data points (as a two dimensional list).
2. n, the number of data points.
3. d, the number of coordinates of each point.
4. k, the desired number of clusters. if k=0, use heuristic to find k*/
static PyObject *fit_nsc_capi(PyObject *self, PyObject *args)
{
    PyObject *unparsed_data_points = NULL, *packed_t_matrix = NULL;
    matrix *data_points = NULL, *t_mat = NULL;
    int n, d, k;

    if (!PyArg_ParseTuple(args, "Oiii", &unparsed_data_points, &n, &d, &k))
    {
        return NULL;
    }
    /*parse data*/
    data_points = parse_python_2Dlist_to_matrix(unparsed_data_points, n, d);
    if (data_points == NULL)
    {
        return NULL;
    }
    /*get the final data points (t matrix)*/
    t_mat = fit_nsc_c(data_points, k);
    if (t_mat == NULL)
    {
        free_matrix(data_points);
        return NULL;
    }
    /*pack final data points (t matrix)*/
    packed_t_matrix = pack_matrix(t_mat);
    if (packed_t_matrix == NULL)
    {
        free_matrix(t_mat);
        free_matrix(data_points);
        return NULL;
    }
    /*free memory*/
    free_matrix(data_points);
    free_matrix(t_mat);

    return packed_t_matrix;
}

/*This array tells Python what methods this module has.*/
static PyMethodDef ckmeansMethods[] = {
    {"calculate_matrix",
     (PyCFunction)calculate_matrix_capi,
     METH_VARARGS,
     PyDoc_STR("The method receives data points (two dimensional python list) and calculates\
     the desired matrix according to the given parameter. The method returns the matrix represented by two dimensional python list.")},
    {"calculate_jacobi",
     (PyCFunction)calculate_jacobi_capi,
     METH_VARARGS,
     PyDoc_STR("The method receives data points (two dimensional python list) and calculates\
     the eigenvalues and the eigenvectors. The method returns a tuple which contains (eigenvalues matrix, eigenvectors matrix).")},
    {"fit_kmeanspp",
     (PyCFunction)fit_kmeanspp_capi,
     METH_VARARGS,
     PyDoc_STR("The method recieves data points, centroids and other configuration parameters and applies the kmeans++ algorithm.\
     Returns the final centroids.")},
    {"fit_nsc",
     (PyCFunction)fit_nsc_capi,
     METH_VARARGS,
     PyDoc_STR("The method recieves data points, and returns the normilized final data points matrix \
     according to the Normalized Spectral Clustering Algorithm.")},
    {NULL, NULL, 0, NULL}};

/*This initiates the module using the above definitions.*/
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL,
    -1,
    ckmeansMethods};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}
