#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#define Mem_Assertion(x) if (!(x)){printf("An Error Has Occurred\n"); exit(1);}//changed to exit() from abort inorder to close open files
#define MAX_ITER 300
#define EPSILON_J 0.00001
#define MAX_ITER_J 100

/*
 * Structs 
 */
typedef struct{
    int count;
    double* newPoints;
    double* centroid;
} Cluster;

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
} EignData;

/*
 * Memory Allocation/Deallocation 
 */
double **matrix_init(int n, int m);
void free_data_points(int n, Point* points);
void free_2D(double **arr);
Point* allocate_mem(int dim, int n);
void data_to_matrix(Point *points, double **matrix, int n, int dim);

/*
 * Print
 */
void print_double(double num);
void print_row(double *row, int n);
void print_column(double *col, int dim, int n);
void print_matrix(double **array, int n, int dim);
void print_Jacobi(double **eign_vectors, double const *eign_values, int n);

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
int get_diag(double const **mat, double *diag,int n);

/*
 * C 
 */
int is_num(char* arg);
Info extractInfo(FILE* file); //changed to file pointer and not path
int processFile(Info info, Point* points, FILE* file);
void get_goal(char *goal, Point *points, int dim, int n);