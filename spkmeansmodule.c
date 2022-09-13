#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"


static PyObject *spkmeans(PyObject *self, PyObject *args)
{
    /* Terminates invalid parameters from spkmeans.py */
    if (PyArg_ParseTuple(args, "isii:spkmeans", &goal_int, &file_name, &n, &d)) {
        /*
        if (goal_int==0){
             This builds the answer ("d" = Convert a C double to a Python floating point number) back into a python object 
            return Py_BuildValue("i", Spkmeans(file_name, n, d));   Py_BuildValue(...) returns a PyObject*  
        } else {
        */
            Matrix M =  MatOutput(goal_int, file_name, n, d);
            k = M.M;
            PyObject* nodes = PyList_New(M.N);
            for (i = 0; i < M.N; ++i){
                PyObject* node = Py_BuildValue("d", M.nodes[i]);
                PyList_SetItem(nodes, i, node);
            }

            PyObject* edges = PyList_New(M.N);
            for (i = 0; i < M.N; ++i){
                PyObject* row = PyList_New(M.M);
                for (j = 0; j < M.M; ++j){
                    PyObject* edge = Py_BuildValue("d", M.edges[i][j]);
                    PyList_SetItem(row, j, edge);
                }
                PyList_SetItem(edges, i, row);
            }

            free_matrix(M);
            return Py_BuildValue("{s:O,s:O,s:i}", "nodes", nodes, "edges", edges, "k", k);
    }
    
   
    
    return NULL;
}


static PyMethodDef capiMethod[] = {
    {"fit",                                                    /* the Python method name that will be used */
     (PyCFunction)spkmeans,                                      /* the C-function that implements the Python function and returns static PyObject*  */
     METH_VARARGS,                                             /* flags indicating parameters accepted for this function */
     PyDoc_STR("spkmeans algorithm for given data of vectors")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}                                      /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans_c", /* name of module */
    NULL,       /* module documentation, may be NULL */
    -1,         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethod  /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC
PyInit_spkmeans_c(void)
{
    return PyModule_Create(&moduledef);
}
