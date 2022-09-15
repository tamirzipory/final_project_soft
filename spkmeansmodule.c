#include <Python.h>
#include "spkmeans.h"


static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *points_of_py, *cent_py;
    PyObject *curr_points, *d_curr, *res, *vector_curr, *curr_cent;
    int cluster_count, vectors_count, vectors_count_alloc, ret;
    int len, K, D, i, j;
    double **Datapoints, **Centroids, **goal_result;
    enum Goal goal;
    curr_cent=NULL, Centroids=NULL;
    
    if (!PyArg_ParseTuple(args, "iiiOiO", &len, &K, &D, &points_of_py, &goal, &cent_py)){
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    vectors_count_alloc=D;
    if(goal==6){
        vectors_count_alloc = D;
        vectors_count_alloc = vectors_count_alloc +1;
    }
    Datapoints = alloc_mat(len, vectors_count_alloc);
    if (Datapoints == NULL){
        /*I dont know if we only need to return NULL, but my friends that we help them say that it's reccomend*/
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    if(goal == 6){
        Centroids = alloc_mat(K, vectors_count_alloc);
        if (Centroids == NULL){
            free_memory(Datapoints, len);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
    }
    i = 0;
    while (i < len)
    {
        curr_points = PyList_GetItem(points_of_py, i);
        if (goal == 6){
            if(i < K)
                curr_cent = PyList_GetItem(cent_py, i);
        }
        j = 0;
        while (j < D)
        {
            d_curr = PyList_GetItem(curr_points, j);
            Datapoints[i][j] = PyFloat_AsDouble(d_curr);
            if (oal == 6){
                if(i < k){
                   d_curr = PyList_GetItem(curr_cent, j);
                   Centroids[i][j] = PyFloat_AsDouble(d_curr);
                }
            }
            j++;
        }
        if(goal == 6){
            Datapoints[i][j] = 0;
            if (i < K)
                Centroids[i][j] = 0;
        }
        i++;
    }
    if(goal != 6){
        goal_result = run_goal(goal, Datapoints, len, D, &K);
        if (goal_result == NULL){
            free_memory(Datapoints, len);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        if(goal == 4){
            cluster_count = len;
            cluster_count = cluster_count +1;
        }
        else cluster_count = len;
        num_of_vectors=len; 
        if(goal==5)
            num_of_vectors=K; 
    } else {
        ret = kMeans(len, K, Datapoints, Centroids, D);
        if (ret == -1){
            free_memory(Datapoints,len);
            free_memory(Centroids,K);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        goal_result = Centroids;
        cluster_count = K; 
        num_of_vectors=D; 
    }
    res = PyList_New(cluster_count);
    i = 0;
    while (i < cluster_count)
    {
        vector_curr = PyList_New(num_of_vectors);
        j = 0;
        while (j < num_of_vectors)
        {
            PyList_SetItem(vector_curr, j, Py_BuildValue("d", goal_result[i][j]));
            j++;
        }
        PyList_SetItem(res, i, Py_BuildValue("O", vector_curr));
        i++;
    }
    free_memory(Datapoints, len);
    free_memory(goal_result, cluster_count);
    return res;
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
