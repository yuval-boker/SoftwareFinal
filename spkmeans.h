#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#define Mem_Assertion(x) if (!(x)){printf("An Error Has Occurred\n"); exit(1);}
#define MAX_ITER 300
#define EPS 0.0
#define EPSILON_J 0.00001
#define MAX_ITER_J 100

/** Structs **/
/* Cluster definition:
 * count = number of vectors assigned to the cluster 
 * newPoints = sum of all vectors assigned to the cluster, instead of matrix of all vectors
 * centroid = current centroid vector
 */
typedef struct{
    int count;
    double* newPoints;
    double* centroid;
} Cluster;

/* 
 * Point definition:
 * vector = the point's data vector
 */
typedef struct{
    double* vector;
} Point;

typedef struct{
    int n;
    int dim;
} Info;

typedef struct{
    double* vector;
    double val;
} EigenData;

/*
 * Memory Allocation/Deallocation 
 */
double **matrix_init(int n, int m);
void free_data_points(int n, Point* points);
void free_2D(double **arr);
Point* allocate_mem(int dim, int n);
void data_to_matrix(Point *points, double **matrix, int n, int dim);
void free_clusters(int k, Cluster* clusters);
void free_memory(int k, int n, Point* points, Cluster* clusters);

/*
 * Print
 */
void print_double(double num);
void print_row(double *row, int n);
void print_column(double *col, int dim, int n);
void print_matrix(double **array, int n, int dim);
void print_Jacobi(double **eigen_vectors, double *eigen_values, int n);

/*
 * WAM
 */
double **create_WAM(Point *points, int dim, int n);
void set_WAM(Point *points, double **matrix, int dim, int n);
double norm_computation(double *v1, double *v2, int dim);

/*
 * DDG 
 */
double **create_DDG(Point *points, int dim, int n);
void set_DDG(double **WAM, double **matrix, int n);
double sum_row(double **matrix, int row, int n);

/*
 * L_norm 
 */
double **create_L_norm(Point *points, int dim, int n);
void set_L_norm(double **WAM, double **DDG, double **L_norm, int n);
double *process_DDG(double **DDG, int n);

/*
 * Jacobi 
 */
void create_Jacobi(Point *points,int dim, int n);
double **I_matrix(int n);
int find_max_indices_off_diag(double **mat, int *i_val, int *j_val, int n);
double sign(double num);
int transform_matrix(double **mat, double **v, int n, int i, int j);
double off_square(double **mat, int n);
double **jacobi(double **A, int n);
void get_diag(double **mat, double *diag, int n);

/*
 * Eigengap Heuristic
 */
void sort(EigenData eigen_arr[], int n);
EigenData *create_sorted_eigen_arr(double **eigen_vectors, double *eigen_values, int n);
int eigengap(EigenData *eigen_arr, int n);

/*
 * Spk
 */
double **create_T(Point *points, int dim, int n, int *k);
double **create_U(EigenData *eigen_arr, int n, int k);
void normalize_U(double **U, int n, int k);
/** Kmeans **/
double distance(Point* point1, const double* centroid, int dim);
int min_distance(Point* point, Cluster* clusters, int dim, int k);
double euclidean_norm(const double* vector, int dim);
void add_point(Point* point, Cluster* cluster, int dim);
int centroid_update(Cluster* cluster, int dim, double *tmp_vector);
int clusters_update(Cluster* clusters, int k, int dim);
void kmeans(int n, Cluster* clusters, Point* points, int dim, int k);

/*
 * C 
 */
Info extractInfo(FILE* file);
int processFile(Info info, Point* points, FILE* file);
void get_goal(char *goal, Point *points, int dim, int n);