#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"


/* 
 * Assigns the current point, given as PyList, into the given vector
 * input: PyList, empty vector in which the data will be saved and dimension
 */
void pyList_to_array(PyObject *list, double* vector, int dim){
    PyObject *item;
    Py_ssize_t coor;
    for(coor = 0; coor < dim; coor++){
        item = PyList_GetItem(list, coor);
        vector[coor] = PyFloat_AsDouble(item);
    }
}

/* 
 * Assigns the n data points, stored in points_list, to the Points in points
 */
void create_matrix(PyObject* points_list, Point* points, int dim, int n) {
    Py_ssize_t i;
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

/*
 * matrix_to_Pylist transformes 2D matrix of doubles to a pylist of pylists containing Python floats
 */
PyObject* matrix_to_PyList(double **arr, int n, int m){
    Py_ssize_t i, j;
    PyObject *PyList = PyList_New((Py_ssize_t)(n));
    PyObject *item;
    for(i = 0; i < n; i++){
        item = PyList_New((Py_ssize_t)(m));
        for(j = 0; j < m; j++){
            PyList_SetItem(item, j, PyFloat_FromDouble(arr[i][j]));
        }
        PyList_SetItem(PyList, i, item);
    }
    return PyList;
}

/*
 * array_to_Pylist transformes kmeans final centroids array of arrays to a pylist of pylists containing Python floats
 * input: clusters = list of clusters, k, dim = dimension
 * output: PyList = a python lists of lists containing final centroids points as PyFloat
 */
PyObject* array_to_PyList(Cluster *clusters, int k, int dim){
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

/* 
 * Allocates memory for k Clusters, and assigns each Cluster a unique centroid, zero vector and count = 0
 * input: dimension and PyList of centroids
 * output: list of k Clusters  
*/
Cluster* createClusters(int dim, PyObject *centroids, int k) {
    Py_ssize_t i;
    Cluster *clusters;
    Cluster clus;
    PyObject *item;
    clusters = (Cluster*) calloc(k, sizeof(Cluster));
    Mem_Assertion(clusters != NULL);
    for (i = 0; i < k; i++) {
        int count = 0;
        double* newPoints = (double*) calloc(dim, sizeof(double));
        double* centroid = (double*) calloc(dim, sizeof(double));
        Mem_Assertion(newPoints != NULL && centroid != NULL);
        item = PyList_GetItem(centroids, i); /* extract the point in index i */
        pyList_to_array(item, centroid, dim); /* assign the point to the new centroid */
        clus.count = count;
        clus.newPoints = newPoints;
        clus.centroid = centroid;
        *(clusters + i) = clus;
    }
    return clusters;
}

/*
 * Transforms array of doubles to 1D Pylist.
 */
PyObject *arr_to_PyList(double *vectors, int n) {
    Py_ssize_t i;
    PyObject *PyList = PyList_New((Py_ssize_t)(n));
    if(!PyList){
        printf("An Error Has Occurred\n");
        return NULL;
    }
    for (i = 0; i < n; i++) {
        PyList_SetItem(PyList, i, PyFloat_FromDouble(vectors[i]));
    }
    return PyList;
}

/*
 * Returns T, given k = 0 returns k aswell
 */
static PyObject *get_T(PyObject *self, PyObject *args){
    Point *points;
    PyObject *data_points, *py_T, *py_tuple;
    double **T;
    int n, dim, k;
    if (!PyArg_ParseTuple(args, "Oiii", &data_points, &n, &dim, &k)) {
        return NULL;
    }
    points = allocate_mem(dim, n);
    create_matrix(data_points, points, dim, n);
    T = create_T(points, dim, n, &k);
    py_T = matrix_to_PyList(T, n, k);
    free_2D(T);
    free_data_points(n, points);
    py_tuple = Py_BuildValue("Oi", py_T, k);
    Py_DECREF(py_T);
    return py_tuple;
}

/* 
 * fit function:
 * input: k, n - the number of points in data, dim - the dimension of each data point 
 * centroids_list - python list of intial centroids, points_list - a python list of lists where each list is a data point in the given data
 * output: raises an error if one or more of the input arguments is not in correct format,
 * otherwise return final- a python list of lists with final centroids caculated by k-means algorithem implemented in kmeans.
 */
static PyObject* fit(PyObject *self, PyObject *args) {
    PyObject *centroids_list, *points_list, *final;
    int k, n, dim;
    Point *points;
    Cluster *clusters;
    if (!PyArg_ParseTuple(args, "iiiOO", &k, &n, &dim, &centroids_list, &points_list)) {
        return NULL;
    }
    points = allocate_mem(dim, n);
    create_matrix(points_list, points, dim, n);
    clusters = createClusters(dim, centroids_list, k);
    kmeans(n, clusters, points, dim, k);
    final = array_to_PyList(clusters, k, dim);
    free_memory(k, n, points, clusters);
    return final;
}

/*
 * Return WAM as Python list 
 */

void free_objects(double** WAM, Point* points, int n){
    free_2D(WAM);
    free_data_points(n, points);
}

static PyObject *get_WAM(PyObject *self, PyObject *args) {
    double **WAM;
    Point *points;
    PyObject *data_points, *py_WAM, *py_list_res;
    int n, dim;
    if (!PyArg_ParseTuple(args, "Oii", &data_points, &n, &dim)) {
        return NULL;
    }
    points = allocate_mem(dim, n);
    create_matrix(data_points, points, dim, n);
    WAM = matrix_init(n, n);
    set_WAM(points, WAM, dim, n);
    py_WAM = matrix_to_PyList(WAM, n, n);
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
/*
 * Jacobi called from Python
 * Returns a tuple containing eigen_values matrix and corresponding eigen_vectors
 */
static PyObject *run_jacobi(PyObject *self, PyObject *args){
    double  **eigen_vectors, **vectors;
    double *eigen_values;
    int n, dim;
    PyObject *py_vectors, *py_values, *py_tuple;
    if(!PyArg_ParseTuple(args, "Oii", &py_vectors, &n, &dim)){
        return NULL;
    }
    // if(!PyList_Check(py_vectors)){
    //     return PyErr_Format(PyExc_TypeError, "Invalid input!");
    // }
    vectors = matrix_init(n, dim);
    eigen_values = calloc(n, sizeof(double));
    if(!eigen_values) {
        printf("An Error Has Occurred\n");
        return NULL;
    }
    pyList_to_matrix(py_vectors, vectors, dim, n);
    eigen_vectors = jacobi(vectors, n);
    py_vectors = matrix_to_PyList(eigen_vectors, n, n);
    get_diag(vectors,eigen_values, n);
    py_values = arr_to_PyList(eigen_values, n);
    free(eigen_values);
    free_2D(eigen_vectors);
    free_2D(vectors);
    py_tuple = Py_BuildValue("OO", py_values, py_vectors);
    Py_DECREF(py_values);
    Py_DECREF(py_vectors);
    return py_tuple;
}

#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
        FUNC(METH_VARARGS, fit, "fit"),
        FUNC(METH_VARARGS, get_T, "get T"),
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