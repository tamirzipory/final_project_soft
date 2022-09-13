#define PY_SSIZE_T_CLEAN /* For all # variants of unit formats (s#, y#, etc.) use Py_ssize_t rather than int. */
#include <Python.h>      /* MUST include <Python.h>, this implies inclusion of the following standard headers: <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <structmember.h>
#include <numpy/arrayobject.h>


#include "spkmeans.h"

/* fit function for kmeansPP */
static PyObject * fit(PyObject *self, PyObject *args)
{
    int err;
	int dim, k, max_iter;
    PyObject *x = NULL;
    PyObject *results = NULL;
    PyObject *centroids = NULL;
    PyArrayObject *x_array = NULL;
    PyArrayObject *results_array = NULL;
    PyArrayObject *centroids_array = NULL;
    npy_intp *shape_x;

    if (!PyArg_ParseTuple(args, "OOOiii", &x, &results, &centroids, &dim, &k, &max_iter)) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }
    x_array = (PyArrayObject*) PyArray_FROM_OTF(x, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (x_array == NULL) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }
    results_array = (PyArrayObject*) PyArray_FROM_OTF(results, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (results_array == NULL) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }
    centroids_array = (PyArrayObject*) PyArray_FROM_OTF(centroids, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (centroids_array == NULL) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }

    shape_x = PyArray_DIMS(x_array);

    /* running the real calculation,
    using pointers to the memory that numpy has allocated */
	
    err = fitToKmeans(
        (double*) PyArray_DATA(x_array),
        (int) shape_x[0],
		dim,
		k,
		max_iter,
		(double*) PyArray_DATA(centroids_array),
        (double*) PyArray_DATA(results_array));


    if (err != 0) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }
    
out:
    Py_XDECREF(x_array);
    Py_XDECREF(results_array);
    Py_XDECREF(centroids_array);
	
	return Py_BuildValue("");
}

/* fit function for SPKmeans */
static PyObject * fitToSPKmeans(PyObject *self, PyObject *args)
{
    int err;
	char *c1;
	char *c2;
	char *c3;
	PyObject *results = NULL;
	PyObject *k_ans = NULL;
    PyArrayObject *results_array = NULL;
    PyArrayObject *k_array = NULL;

    if (!PyArg_ParseTuple(args, "sssOO", &c1, &c2, &c3, &results, &k_ans)) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }
	
	results_array = (PyArrayObject*) PyArray_FROM_OTF(results, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (results_array == NULL) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }
	
	k_array = (PyArrayObject*) PyArray_FROM_OTF(k_ans, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (results_array == NULL) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }
	
    err = do_main(c1, 
	c2, 
	c3, 
	(double*) PyArray_DATA(results_array),
	(double*) PyArray_DATA(k_array));

    if (err != 0) {
        PyErr_Format(PyExc_ValueError, "An Error Has Occured");
        goto out;
    }

out:
	return Py_BuildValue("");
}


static PyMethodDef demo_methods[] = {
    {"fit", (PyCFunction) fit, METH_VARARGS,
         "Runs an algorithm defined in a local C file."},
    {"fitToSPKmeans", (PyCFunction) fitToSPKmeans, METH_VARARGS,
         "Runs an algorithm defined in a local C file."},
    {NULL, NULL, 0}   /* sentinel */
};

static struct PyModuleDef demomodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    demo_methods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *mod = NULL;
    import_array();
    mod = PyModule_Create(&demomodule);
    return mod;
}
