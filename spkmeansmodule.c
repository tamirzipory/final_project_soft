#include <Python.h>
#include "spkmeans.h"


static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *Datapoints_PyObject, *Centroids_PyObject, *current_datapoint, *current_double, *returned_result, *current_vector, *current_centroid;

    int len, K, D, i, j, num_of_clusters, num_of_vectors, num_of_vectors_allocation, return_value;
    double **Datapoints, **Centroids, **goal_result;
    enum Goal goal;
    current_centroid=NULL, Centroids=NULL;
    

    if (!PyArg_ParseTuple(args, "iiiOiO", &len, &K, &D, &Datapoints_PyObject, &goal, &Centroids_PyObject)){
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    num_of_vectors_allocation=D;
    if(goal==6)
        num_of_vectors_allocation=D+1;

    Datapoints = alloc_for_mat(len, num_of_vectors_allocation);
    if (Datapoints == NULL){
        /*I dont know if we only need to return NULL, but my friends that we help them say that it's reccomend*/
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    
   
    if(goal == 6){
        Centroids = alloc_for_mat(K, num_of_vectors_allocation);
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
        if (i < K && goal == 6)
            current_centroid = PyList_GetItem(Centroids_PyObject, i);

        /*Set up each of vector*/
        for (j = 0; j < D; j++){
            current_double = PyList_GetItem(current_datapoint, j);
            Datapoints[i][j] = PyFloat_AsDouble(current_double);
            if (i < K && goal == 6){
                current_double = PyList_GetItem(current_centroid, j);
                Centroids[i][j] = PyFloat_AsDouble(current_double);
            }
        }

      
        if(goal == 6){
            Datapoints[i][j] = 0;
            if (i < K)
                Centroids[i][j] = 0;
            
        }
    }

  
    if (goal == 6){
        return_value = kMeans(len, K, Datapoints, Centroids, D);
        if (return_value == -1){
            free_memory(Datapoints,len);
            free_memory(Centroids,K);
            PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
            return NULL;
        }
        goal_result = Centroids;
        num_of_clusters = K; 
        num_of_vectors=D; 
    }
    else
    {
        goal_result = target_runner(goal, Datapoints, len, D, &K);
        if (goal_result == NULL){
            free_memory(Datapoints, len);
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

    /* Converts result_matrix to an array list (python)*/
    returned_result = PyList_New(num_of_clusters);
    for (i = 0; i < num_of_clusters; ++i){
        current_vector = PyList_New(num_of_vectors);
        for (j = 0; j < num_of_vectors; j++)
            PyList_SetItem(current_vector, j, Py_BuildValue("d", goal_result[i][j]));
        PyList_SetItem(returned_result, i, Py_BuildValue("O", current_vector));
    }

    free_memory(Datapoints, len);
    free_memory(goal_result, num_of_clusters);

    return returned_result;
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
