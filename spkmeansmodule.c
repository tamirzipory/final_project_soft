#include <Python.h>
#include "spkmeans.h"

static PyObject *fit(PyObject *self, PyObject *args){
    PyObject *points_op_py, *cent_of_python, *d_points_curr, *curr_of_number, *arr_of_output, *vector_of_now, *arr_cent_of_now;

    int len, K, dis, i, j
    int num_of_clusters, num_of_vectors, num_vec_alloc, what_client_get;
    double **mat_of_points, **mat_of_cent, **mat_of_goals;
    enum Goal goal;
    arr_cent_of_now=NULL, mat_of_cent=NULL;
    if (!PyArg_ParseTuple(args, "iiiOiO", &len, &K, &dis, &points_op_py, &goal, &cent_of_python)){
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    num_vec_alloc=dis;
    if(goal==6)
        num_vec_alloc=dis+1;
    mat_of_points = alloc_mat(len, num_vec_alloc);
    if (mat_of_points == NULL){
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    if(goal == 6){
        mat_of_cent = alloc_mat(K, num_vec_alloc);
        if (mat_of_cent == NULL){
            free_memory(mat_of_points, len);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
    }
     i=0;
     while(i < len){
        d_points_curr = PyList_GetItem(points_op_py, i);
        if (goal == 6){
            if(i < K)
            arr_cent_of_now = PyList_GetItem(cent_of_python, i);
        }
        j = 0;
        while (j < dis)
        {
            curr_of_number = PyList_GetItem(d_points_curr, j);
            mat_of_points[i][j] = PyFloat_AsDouble(curr_of_number);
            if (goal == 6){
                if(i < k){
                curr_of_number = PyList_GetItem(arr_cent_of_now, j);
                mat_of_cent[i][j] = PyFloat_AsDouble(curr_of_number);
                }
            }
            j++;
        }

        if(goal == 6){
            mat_of_points[i][j] = 0;
            if (i < K) mat_of_cent[i][j] = 0;
        }
        i++;
     }
    if (goal == 6){
        what_client_get = kMeans(len, K, mat_of_points, mat_of_cent, dis);
        if (what_client_get == -1){
            free_memory(mat_of_points,len);
            free_memory(mat_of_cent,K);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        mat_of_goals = mat_of_cent;
        num_of_clusters = K; 
        num_of_vectors=dis; 
    }
    else
    {
        mat_of_goals = run_goal(goal, mat_of_points, len, D, &K);
        if (mat_of_goals == NULL){
            free_memory(mat_of_points, len);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }

        if(goal == 4)
            num_of_clusters = len+1;
        else
            num_of_clusters = len;

        num_of_vectors=len; 

        if(goal==5)
            num_of_vectors=K; 
    }

    arr_of_output = PyList_New(num_of_clusters);
    for (i = 0; i < num_of_clusters; ++i){
        vector_of_now = PyList_New(num_of_vectors);
        for (j = 0; j < num_of_vectors; j++)
            PyList_SetItem(vector_of_now, j, Py_BuildValue("d", mat_of_goals[i][j]));
        PyList_SetItem(arr_of_output, i, Py_BuildValue("O", vector_of_now));
    }

    free_memory(mat_of_points, len);
    free_memory(mat_of_goals, num_of_clusters);

    return arr_of_output;
}


static PyMethodDef Methods[] = {
        {"fit",
                (PyCFunction)fit,
                     METH_VARARGS,
                        NULL},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moudledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeans_module",
        NULL,
        -1,
        Methods
};

PyMODINIT_FUNC PyInit_spkmeans_module(void){
    PyObject *m;
    m = PyModule_Create(&moudledef);
    if (!m)
        return NULL;
    
    return m;
}
