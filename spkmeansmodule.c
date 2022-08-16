#define PY_SSIZE_T_CLEAN
#include <Python.h>
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
    PyObject *item;
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
            item = PyList_GetItem(row, j);
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


PyObject* matrix_to_PyList(double **arr, int n, int k){
    Py_ssize_t i, j;
    PyObject *PyList = PyList_New(n);
    PyObject *item;
    for(i = 0; i < n; i++){
        item = PyList_New(k);
        for(j = 0; j < k; j++){
            printf("value: %f\n", arr[i][j]);
            PyList_SetItem(item, j, PyFloat_FromDouble(arr[i][j]));
        }
        PyList_SetItem(PyList, i, item);
    }
    return PyList;
}

/*
 * Return WAM as Python list 
 */

void free_objects(double** WAM, Point* points, int n){
    free_2D(WAM);
    free_data_points(n, points);
    printf("All objects are free\n");
    WAM = NULL;
    points = NULL;
    printf("points pointer is: %p\n", points);
    printf("WAM pointer is: %p\n", WAM);
}

static PyObject *get_WAM(PyObject *self, PyObject *args) {
    double **WAM;
    Point *points;
    PyObject *data_points, *py_WAM;
    int n, dim;
    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &dim)) {
        return NULL;
    }
    points = allocate_mem(dim, n);
    printf("points pointer is: %p\n", points);
    create_matrix(data_points, points, dim, n);
    WAM = matrix_init(n, n);
    printf("WAM pointer is: %p\n", WAM);
    set_WAM(points, WAM, dim, n);
    py_WAM = NULL;
    py_WAM = matrix_to_PyList(WAM, n, n);
    printf("py_WAM: %p\n", py_WAM);
    free_objects(WAM, points, n);
    return py_WAM;
}

/*
 * Return DDG as Python list 
 */
static PyObject *get_DDG(PyObject *self, PyObject *args) {
    double **WAM, **DDG;
    Point *points;
    PyObject *data_points, *py_DDG;
    int n, dim;
    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &dim)) {
        return NULL;
    }
    points = allocate_mem(dim, n);
    create_matrix(data_points, points, dim, n);
    WAM = matrix_init(n, n);
    DDG = matrix_init(n, n);
    set_WAM(points, WAM, dim, n);
    set_DDG(WAM, DDG, n);
    py_DDG = matrix_to_PyList(DDG, n, n);
    free_2D(WAM);
    free_2D(DDG);
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
    int n, dim;
    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &dim)) {
        return NULL;
    }
    points = allocate_mem(dim, n);
    create_matrix(data_points, points, dim, n);
    WAM = matrix_init(n, n);
    DDG = matrix_init(n, n);
    L_norm = matrix_init(n, n);
    set_WAM(points, WAM, dim, n);
    set_DDG(WAM, DDG, n);
    set_L_norm(WAM, DDG, L_norm, n);
    py_L_norm = matrix_to_PyList(L_norm, n, n);
    free_2D(WAM);
    free_2D(DDG);
    free_2D(L_norm);
    free_data_points(n, points);
    return py_L_norm;
}

static PyObject *run_jacobi(PyObject *self, PyObject *args){
    double  **eigen_vectors, **vectors;
    double *eigen_values;
    int n, dim;
    PyObject *py_vectors;
//    *res_values, *res_vectors, *res_tuple;
    if(!PyArg_ParseTuple(args, "Oii", &py_vectors, &n, &dim)){
        return NULL;
    }
    if(!PyList_Check(py_vectors)){
        return PyErr_Format(PyExc_TypeError, "Invalid input!");
    }
    vectors = matrix_init(n, dim);
    eigen_values = calloc(n, sizeof(double));
    if(!eigen_values) {
        return PyErr_NoMemory();
    }
    pyList_to_matrix(py_vectors, vectors, dim, n);
    eigen_vectors = jacobi(vectors, n);
    if (!eigen_vectors){
        return PyErr_NoMemory();
    }
    get_diag(vectors,eigen_values, n);
    print_Jacobi(eigen_vectors, eigen_values, n);
    free(eigen_values);
    free_2D(eigen_vectors);
    free_2D(vectors);
    return 0;
}
#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
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