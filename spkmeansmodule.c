#include <Python.h>
#include "spkmeans.h"

/* Receives N, K, D, Datapoints/matrix, goal, Centroids from python
 * Returns correspond matrix according to given goal*/
static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *Datapoints_PyObject, *Centroids_PyObject, *current_datapoint, *current_double, *returned_result, *current_vector, *current_centroid;

    int len, K, D, i, j, rows, cols, cols_allocation, return_value;
    double **Datapoints, **Centroids, **goal_result;
    enum Goal goal;
    current_centroid=NULL, Centroids=NULL;
    

    if (!PyArg_ParseTuple(args, "iiiOiO", &len, &K, &D, &Datapoints_PyObject, &goal, &Centroids_PyObject)){
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    cols_allocation=D;
    /* Only in spk-ex2 (when we return to fit in the second time), Datapoints need to have D+1 cols*/
    if(goal==6)
        cols_allocation=D+1;

    /* Set up Datapoints and Centroids's matrix*/
    Datapoints = matrix_allocation(len, cols_allocation);
    if (Datapoints == NULL){
        /*I dont know if we only need to return NULL, but my friends that we help them say that it's reccomend*/
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    
   
    if(goal == 6){
        Centroids = matrix_allocation(K, cols_allocation);
        if (Centroids == NULL){
            free_memory(Datapoints, len);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
    }
    
    /* Fill matrix values as given list from python*/
    for (i = 0; i < len; i++)
    {
        current_datapoint = PyList_GetItem(Datapoints_PyObject, i);
        if (i < K && goal == SPK_EX2)
            current_centroid = PyList_GetItem(Centroids_PyObject, i);

        /*Set up each of vector*/
        for (j = 0; j < D; j++)
        {
            current_double = PyList_GetItem(current_datapoint, j);
            Datapoints[i][j] = PyFloat_AsDouble(current_double);
            if (i < K && goal == SPK_EX2)
            {
                current_double = PyList_GetItem(current_centroid, j);
                Centroids[i][j] = PyFloat_AsDouble(current_double);
            }
        }

        /* Only in spk-ex2: Zero in last cell [dimension]*/
        if(goal == 6)
        {
            Datapoints[i][j] = 0;
            if (i < K)
            
                Centroids[i][j] = 0;
            
        }
    }

    /* If goal is spk_ex2 (spk second run) run kMeans from ex2 else- goal is wam, ddg, lnorm or spk (first run)- use run_goal function*/
    if (goal == 6){
        return_value = kMeans(len, K, Datapoints, Centroids, D);
        if (return_value == -1){
            free_memory(Datapoints,len);
            free_memory(Centroids,K);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        goal_result = Centroids;
        rows = K; /*rows=number of clusters/centroids=K*/
        cols=D; /*cols= dimension*/
    }
    else
    {
        goal_result = run_goal(goal, Datapoints, len, D, &K);
        if (goal_result == NULL){
            free_memory(Datapoints, len);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        rows = (goal == 4) ? (len + 1) : len; /*Jacobi needs N+1 rows*/
        cols=len; /* In wam,ddg,lnorm,jacobi*/
        if(goal==5)
            cols=K; /* In spk (first run)- T's dimesnions are N*K (updated/original k) */
    }

    /* Converts result_matrix to an array list (python)*/
    returned_result = PyList_New(rows);
    for (i = 0; i < rows; ++i){
        current_vector = PyList_New(cols);
        for (j = 0; j < cols; j++)
            PyList_SetItem(current_vector, j, Py_BuildValue("d", goal_result[i][j]));
        PyList_SetItem(returned_result, i, Py_BuildValue("O", current_vector));
    }

    free_memory(Datapoints, len);
    free_memory(goal_result, rows);

    return returned_result;
}

/* ============== C Python API ============== */
static PyMethodDef Methods[] = {
        {"fit",
                (PyCFunction)fit,
                     METH_VARARGS,
                        NULL},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moudledef = {
        PyModuleDef_HEAD_INIT,
        "my_spkmeans",
        NULL,
        -1,
        Methods
};

PyMODINIT_FUNC PyInit_my_spkmeans(void){
    PyObject *m;
    m = PyModule_Create(&moudledef);
    if (!m)
        return NULL;
    
    return m;
}
