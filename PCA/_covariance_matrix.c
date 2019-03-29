#include <Python.h>
#include <numpy/arrayobject.h>
#include "covariance_matrix.h"

static char module_docstring[] = 
    "This module provides an interface for calculating chi-squared using C."
static char coefficients_docstring[] = 
    "Compute the coefficients of the H0 expansion"

static PyObject *covariance_matrix_coefficients(PyObject *self, PyObject *args)

static PyMethodDef module_methods[] = {
  {"coefficients", covariane_matrix_coefficients, METH_VARARGS, coefficients_docstring},
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC init_chi2(void)
{
    PyObject *m = Py_InitModule3("_covariance_matrix", module_methods, module_docstring);
    if (m==NULL)
      return;

    /* Load `numpy` functionaluty. */
    import_array()
}



static PyObject *covariance_matrix_coefficients(PyObject *self, PyObject *args)
{
    int n_points, nmax, lmax;
    PyObject *r, *theta, *phi, *M;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dd00", &n_points, &r, &theta, &phi, &M, &nmax, &lmax))
      return NULL;
    /* interpret the input objects as numpy arrays */
    PyObject *r = PyArray_FROM_OTF(r, NPY_double, NPY_IN_ARRAY);
    PyObject *theta = PyArray_FROM_OTF(theta, NPY_double, NPY_IN_ARRAY);
    PyObject *phi = PyArray_FROM_OTF(phi, NPY_double, NPY_IN_ARRAY);
    PyObject *M = PyArray_FROM_OTF(M, NPY_double, NPY_IN_ARRAY);


}






