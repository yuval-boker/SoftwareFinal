#define PY_SSIZE_T_CLEAN
#include </System/Library/Frameworks/Python.framework/versions/2.7/include/python2.7/Python.h>
#include "spkmeans.h"



void pyList_to_array(PyObject *list, double* vector, int dim){
    PyObject *item;
    Py_ssize_t coor;
    for(coor = 0; coor < dim; coor++){
        item = PyList_GetItem(list, coor);
        vector[coor] = PyFloat_AsDouble(item);
    }
}


void create_matrix(PyObject* points_list, Point* points, int dim, int n) {
    int i;
    PyObject *item, *coordinate;
    for (i = 0; i < n; i++) {
        item = PyList_GetItem(points_list, i);
        pyList_to_array(item, points[i].vector, dim);
    }
}

/*
 * Transfers 2D PyList of floats to 2D array of doubles.
 */
void pyList_to_matrix(PyObject* vector_list, double** a, int dim, int n){
    Py_ssize_t i, j;
    PyObject *row, *item;
    for(i = 0; i < n; i++){
        row = PyList_GetItem(vector_list, i);
        for(j = 0; j < dim; j++){
            item = PyList_GetItem(item, j);
            a[i][j] = PyFloat_AsDouble(item);
        }
    }
}

PyObject* clusters_to_PyList(Cluster *clusters, int k, int dim){
    Py_ssize_t i, j;
    PyObject *PyList = PyList_New((Py_ssize_t)(k));
    PyObject *item;
    for(i = 0; i < k; i++){
        item = PyList_New((Py_ssize_t)(dim));
        for(j = 0; j < dim; j++){
            PyList_SetItem(item, j, PyFloat_FromDouble(clusters[i].centroid[j]));
        }
        PyList_SetItem(PyList, i, item);
    }
    return PyList;
}


PyObject* matrix_to_PyList(double **arr, int n){
    Py_ssize_t i, j;
    PyObject *PyList = PyList_New((Py_ssize_t)(n));
    PyObject *item;
    for(i = 0; i < n; i++){
        item = PyList_New((Py_ssize_t)(n));
        for(j = 0; j < n; j++){
            PyList_SetItem(item, j, PyFloat_FromDouble(arr[i][j]));
        }
        PyList_SetItem(PyList, i, item);
    }
    return PyList;
}

/*
 * Return WAM as Python list 
 */
static PyObject *get_WAM(PyObject *self, PyObject *args) {
    double **WAM;
    Point *points;
    PyObject *data_points, *py_WAM;
    Py_ssize_t n, dim;
    if (!PyArg_ParseTuple(args, "(Oii):wam", &data_points, &n, &dim))
        return NULL;
    if (!PyList_Check(data_points))
        return PyErr_Format(PyExc_TypeError, "Invalid input!");
    points = allocate_mem(dim, n);
    create_matrix(data_points, points, dim, n);
    WAM = matrix_init(n, n);
    assert(WAM != NULL);
    set_WAM(points, WAM, dim, n);
    py_WAM = matrix_to_PyList(WAM, n);
    free_2D(WAM, n);
    free_data_points(n, points);
    return py_WAM;
}

/*
 * Return DDG as Python list 
 */
static PyObject *get_DDG(PyObject *self, PyObject *args) {
    double **WAM, **DDG;
    Point *points;
    PyObject *data_points, *py_DDG;
    Py_ssize_t n, dim;
    if (!PyArg_ParseTuple(args, "(Oii):wam", &data_points, &n, &dim))
        return NULL;
    if (!PyList_Check(data_points))
        return PyErr_Format(PyExc_TypeError, "Invalid input!");
    points = allocate_mem(dim, n);
    create_matrix(data_points, points, dim, n);
    WAM = matrix_init(n, n);
    DDG = matrix_init(n, n);
    assert(WAM != NULL || DDG != NULL);
    set_WAM(points, WAM, dim, n);
    set_DDG(WAM, DDG, n);
    py_DDG = matrix_to_PyList(DDG, n);
    free_2D(WAM, n);
    free_2D(DDG, n);
    free_data_points(n, points);
    return py_DDG;
}

/*
 * Return L_norm as Python list 
 */
static PyObject *get_L_norm(PyObject *self, PyObject *args) {
    double **WAM, **DDG, **L_norm;
    Point *points;
    PyObject *data_points, *py_L_norm;
    Py_ssize_t n, dim;
    if (!PyArg_ParseTuple(args, "(Oii):wam", &data_points, &n, &dim))
        return NULL;
    if (!PyList_Check(data_points))
        return PyErr_Format(PyExc_TypeError, "Invalid input!");
    points = allocate_mem(dim, n);
    create_matrix(data_points, points, dim, n);
    WAM = matrix_init(n, n);
    DDG = matrix_init(n, n);
    L_norm = matrix_init(n, n);
    assert(WAM != NULL || DDG != NULL || L_norm != NULL);
    set_WAM(points, WAM, dim, n);
    set_DDG(WAM, DDG, n);
    set_L_norm(WAM, DDG, L_norm, n);
    py_L_norm = matrix_to_PyList(L_norm, n);
    free_2D(WAM, n);
    free_2D(DDG, n);
    free_2D(L_norm, n);
    free_data_points(n, points);
    return py_L_norm;
}

static PyObject *run_jacobi(PyObject *self, PyObject *args){
    double  **eigen_vectors, **vectors;
    double *eigen_values;
    int n, dim;
    PyObject *py_vectors, *res_values, *res_vectors, *res_tuple;
    if(!PyArg_ParseTuple(args, "(Oii):jacobi", &py_vectors, &n, &dim)){
        return NULL;
    }
    if(!PyList_Check(py_vectors)){
        return PyErr_Format(PyExc_TypeError, "Invalid input!");
    }
    vectors = matrix_init(n, dim);
    eigen_vectors = matrix_init(n,n);
    eigen_values = calloc(n, sizeof(double));
    if(!vectors || !eigen_vectors || !eigen_values) { //not sure what to do
        return PyErr_NoMemory();
    }
    PyList_to_matrix(py_vectors, vectors, n, dim);
    eigen_vectors = jacobi(vectors, n);
    if (!eigen_vectors){
        return PyErr_NoMemory();
    }
    eigen_values = get_diag(vectors,eigen_values, n);
    if(!eigen_values){
       return PyErr_NoMemory(); 
    }
    print_Jacobi(eigen_vectors, eigen_values, n);
    free(eigen_values);
    free_2D(eigen_vectors, n);
    free_2D(vectors, n);
    //return? 
}
#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
//        FUNC(METH_VARARGS, fit_kmeans, "run kmeans"),
//        FUNC(METH_VARARGS, get_T, "get T"),
        FUNC(METH_VARARGS, get_WAM, "get WAM"),
        FUNC(METH_VARARGS, get_DDG, "get DDG"),
        FUNC(METH_VARARGS, get_L_norm, "get L_norm"),
        FUNC(METH_VARARGS, run_jacobi, "run jacobi"),
        {NULL, NULL, 0, NULL}   /* sentinel */
};


static struct PyModuleDef _moduledef = {
       PyModuleDef_HEAD_INIT,
       "spkmeans_capi",
       NULL,
       -1,
       _methods
};


PyMODINIT_FUNC
PyInit_spkmeans_capi(void) {
   PyObject *m;
   m = PyModule_Create(&_moduledef);
   if (!m) {
       return NULL;
   }
   return m;
}