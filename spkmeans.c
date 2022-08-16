#include "spkmeans.h"


/** print funcs **/
void print_double(double num) {
    if (num < 0 && num > -0.00005)
        num = 0.0;
    printf("%.4f", num);
}

void print_row(double *row, int n) {
    int i;
    for (i = 0; i < n; i++) {
        print_double(row[i]);
        if (i < n - 1)
            printf("%c", ',');
    }
}

void print_column(double *col, int dim, int n) {
    int i;
    for (i = 0; i < n; i++) {
        print_double(col[i * dim]);
        if (i < n - 1)
            printf("%c", ',');
    }
}


void print_matrix(double **array, int n, int dim) {
    int i;
    for (i = 0; i < n; i++) {
        print_row(array[i], dim);
        if (i < n - 1) //change to i < n?
            printf("\n");
    }
}

/*
 * Prints the eigenvalues as the first line, second line onward is the corresponding eigenvectors
 */
void print_Jacobi(double **eign_vectors, double *eign_values, int n){ //may change const
    print_row(eign_values, n);
    printf("\n");
    print_matrix(eign_vectors, n, n);
}

/*
 * Allocates memory for n data points stored as matrix of Point 
 */
Point* allocate_mem(int dim, int n){
    int i;
    double* v;
    Point* points = (Point*)calloc(n, sizeof(Point));
    Mem_Assertion(points != NULL);
    for(i = 0; i < n; i++){
        v = (double*) calloc(dim,sizeof(double));
        Mem_Assertion(v != NULL);
        points[i].vector = v;
    }
    return points;
}

/*
 * Deallocates the memory that was previously allocated 
 */
void free_data_points(int n, Point* points){
    int i;
    if (points != NULL){
        for (i = 0; i < n; i++){
            free(points[i].vector);
        }
        free(points);
    }
}

/*
 * Deallocates the memory that was previously dynamically allocated 
 */
void free_2D(double **matrix) {
    free(matrix[0]);
    free(matrix);
}

/*
 * Euclidean norm computation 
 */
double norm_computation(double *v1, double *v2, int dim) {
    double res = 0.0;
    double diff;
    int i;
    for (i = 0; i < dim; i++) {
        diff = v1[i] - v2[i];
        res += (diff * diff);
    }
    return res;
}

/*
 * Sets WAM using Euclidean norm computation 
 */
void set_WAM(Point *points, double **matrix, int dim, int n) {
    int i, j;
    double w_ij, tmp;
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            tmp = sqrt(norm_computation(points[i].vector, points[j].vector, dim));
            w_ij = exp((-0.5) * tmp);
            matrix[i][j] = w_ij;
            matrix[j][i] = w_ij;
        }
    }
}


/*
 * Allocates a matrix dynamically (as an array of arrays)
 * The matrix will be stored in memory continuously, as shown in class 
 */
double **matrix_init(int row, int col){
    double *content, **matrix;
    int i;
    content = (double*)calloc(row*col, sizeof(double));
    Mem_Assertion(content != NULL);
    matrix = (double**) calloc(row, sizeof(double*));
    Mem_Assertion(matrix != NULL);
    for (i = 0; i < row; i++){
        matrix[i] = content + (i * col);
    }
    return matrix;
}

/*
 * Sums the values in a given matrix row 
 */
double sum_row(double **matrix, int row, int n){
    int i;
    double res = 0.0;
    for (i = 0; i < n; i++){
        res += matrix[row][i];
    }
    return res;
}

/** JACOBI **/
/*
 * Tranforms data structure from array of points to 2D matrix for design reasons 
 */
void data_to_matrix(Point *points, double **matrix, int n, int dim){
    int i, j;
    for (i = 0; i < n; i++){
        for (j = 0; j < dim; j++){
            matrix[i][j] = points[i].vector[j];
        }
    }
}

/*
 * creates nxn Identity matrix
 */ 
double **I_matrix(int n){
    double **mat;
    int i;
    mat = matrix_init(n,n);
    for(i = 0; i < n; i++){
        mat[i][i]= 1.0;
    }
    return mat;
}

/*
 * Finds max off diagonal element's indices and updates i_val and j_val accordingly
 */
int find_max_indices_off_diag(double **mat, int *i_val, int *j_val, int n){
    int i,j;
    double max_val = 0.0;
    for(i = 0; i < n ; i++){
        for(j = i + 1; j < n; j++){
            if(fabs(mat[i][j]) > max_val){
                max_val = fabs(mat[i][j]);
                *i_val = i;
                *j_val = j;
            }
        }
    }
    return 0;
}

/*
 * Return the sign of given num
 */
double sign(double num){
    if(num < 0){
        return -1.0;
    }
    return 1.0;
}

/* Transforms A as described in 1.2.1 and updates V
 * we used these refrences:
 * 1. http://phys.uri.edu/nigh/NumRec/bookfpdf/f11-1.pdf
 * 2. https://github.com/mateuv/MetodosNumericos/blob/master/python/NumericalMethodsInEngineeringWithPython/jacobi.py
 */
int transform_matrix(double **mat, double **v, int n, int i, int j){
    int r;
    double tmp1, tmp2, theta, t, c, s;
    theta = (mat[j][j] - mat[i][i])/(2 * mat[i][j]);
    t = sign(theta)/(fabs(theta) + sqrt((theta * theta)+1.0));
    c = 1.0/(sqrt((t * t) + 1.0));
    s = t * c;

    for(r = 0; r < n; r++){
        if (r != i && r != j){
            tmp1 = mat[r][i]; /* Ari */
            tmp2 = mat[r][j]; /* Arj */
            mat[r][i] = (c * tmp1) - (s * tmp2);
            mat[i][r] = (c * tmp1) - (s * tmp2); /* symmetry */
            mat[r][j] = (c * tmp2) + (s * tmp1);
            mat[j][r] = (c * tmp2) + (s * tmp1); /* symmetry */
        }
    }

    tmp1 = mat[i][i]; /* Aii */
    tmp2 = mat[j][j]; /* Ajj */
    mat[i][i] = ((c * c) * tmp1 ) + ((s * s) * tmp2) - (2 * s * c * mat[i][j]);
    mat[j][j] = ((s * s) * tmp1 ) + ((c * c) * tmp2) + (2 * s * c * mat[i][j]);
    mat[i][j] = 0.0;
    mat[j][i] = 0.0;

    for(r = 0; r < n; r++){ /* updates V */
        tmp1 = v[r][i];
        tmp2 = v[r][j];
        v[r][i] = (c * tmp1) - (s * tmp2);
        v[r][j] = (s * tmp1) + (c * tmp2);
    }
    return 0;
}

/* 
 * Calculate sum of squares of all off-diagonal elements of given matrix
 * going through half of matrix elments due to symmetry
 */
double off_square(double **mat, int n) {
    int i, j;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            sum += 2 * (mat[i][j] * mat[i][j]);
        }
    }
    return sum;
}
/*
 * creates jacobi
 */
double **jacobi(double **A, int n) {
    double **V;
    double off_A = 0.0, off_A_prime;
    int i, j, l = 0;
    V = I_matrix(n); /* initializing V to Identity matrix */
    while (l < MAX_ITER_J) {
        if (l == 0) {
            off_A = off_square(A, n);
        }
        find_max_indices_off_diag(A, &i, &j, n);
        transform_matrix(A, V, n, i, j);
        off_A_prime = off_square(A, n);
        if((l != 0) && ((off_A - off_A_prime) <= EPSILON_J)){
            return V;
        }
        off_A = off_A_prime;
        l++;
    }
    return V;
}

/*
* Stores given matrix diagonal values in a 1D array
*/
void get_diag(double **mat, double *diag, int n){
    int i;
    for(i = 0; i < n; i++){
        diag[i] = mat[i][i];
    }
}

/** EIGENGAP HEURISTIC **/
/*
* Given a EignData array sorts inplace by eignvalues in decreasing order.
* we used this reference:
* https://www.javatpoint.com/c-program-to-sort-the-elements-of-an-array-in-descending-order
*/
void sort(EignData eign_arr[], int n){
    int i, j;
    EignData temp;
    for(i = 0; i < n; i++){
        for(j = i + 1; j < n; j++){
            if(eign_arr[i].val < eign_arr[j].val){
                temp = eign_arr[i];
                eign_arr[i] = eign_arr[j];
                eign_arr[j] = temp;
            }
        }
    }
}

/*
 * Given an n eignvectors 2D array and a 1D array of n corresponding eignvalues,
 * creates a new 1D EignData (contains eignvalue it's eignvector) array
 * which is decreasingly ordered by eignvalues and returns the new array.
 */
EignData *create_sorted_eign_arr(double **eign_vectors, double const *eign_values, int n){
    int i;
    EignData* eign_arr =(EignData*)calloc(n, sizeof(EignData));
    Mem_Assertion(eign_arr != NULL);
    for(i = 0; i < n; i++){
        eign_arr[i].val = eign_values[i];
        eign_arr[i].vector = eign_vectors[i]; //check with goni that its not &
    }
    sort(eign_arr, n);
    return eign_arr;
}

/*
 * Given an array of n EignData elements (each element contains eignvalue it's eignvector)
 * implements eigengap heuristic as described in 1.3 and returns k.
 */
int eigngap(EignData *eign_arr, int n){
    int i, max_i = 0;
    double curr_delta, max_delta = 0.0;
    for(i = 0; i < (n / 2); i++){
        curr_delta = fabs(eign_arr[i].val - eign_arr[i + 1].val);
        if(curr_delta > max_delta){
            max_delta = curr_delta;
            max_i = i + 1; /* i starts from 1 in eigngap */
        }
    }
    return max_i;
}

/*
 * Creates WAM using points 
 */
double **create_WAM(Point *points, int dim, int n){
    double **matrix;
    matrix = matrix_init(n, n);             /* WAM ∈ R^(n×n) */
    set_WAM(points, matrix, dim, n);
    return matrix;
}

/*
 * Sets DDG using WAM
 */
void set_DDG(double **WAM, double **matrix, int n){
    int i;
    for (i = 0; i < n; i++){
        matrix[i][i] = sum_row(WAM, i, n);
    }
}

/*
 * Creates DDG using WAM creation function 
 */
double **create_DDG(Point *points, int dim, int n){
    double **matrix_ddg, **matrix_wam;
    matrix_wam = create_WAM(points, dim, n);
    matrix_ddg = matrix_init(n, n);            /* DDG ∈ R^(n×n) */
    set_DDG(matrix_wam, matrix_ddg, n);
    free_2D(matrix_wam);
    return matrix_ddg;
}

/*
 * Based on DDG's diagonal, creates an array 
 */
double *process_DDG(double **DDG, int n){
    double *res;
    int i;
    res = (double*) calloc(n,sizeof(double));
    Mem_Assertion(res != NULL);
    for (i = 0; i < n; i++){
        res[i] = pow(DDG[i][i], -0.5);
    }
    return res;
}

/*
 * Sets l_norm using WAM and DDG
 */
void set_L_norm(double **WAM, double **DDG, double **L_norm, int n){
    double *D;
    int i, j;
    D = process_DDG(DDG, n);
    for (i = 0; i < n; i++){
        for (j = 0; j <= i; j++){
            if (j < i){
                L_norm[i][j] = (-1) * D[i] * WAM[i][j] * D[j];
                L_norm[j][i] = L_norm[i][j];
            } else {
                L_norm[i][j] = 1 - (D[i] * WAM[i][j] * D[j]);
            }
        }
    }
    free(D);
}

/*
 * Creates l_norm using WAM's and DDG's creation functions 
 */
double **create_L_norm(Point *points, int dim, int n) {
    double **matrix_L_norm, **matrix_ddg, **matrix_wam;
    matrix_wam = create_WAM(points, dim, n);
    matrix_ddg = create_DDG(points, dim, n);
    matrix_L_norm = matrix_init(n, n);            /* L_norm ∈ R^(n×n) */
    set_L_norm(matrix_wam, matrix_ddg, matrix_L_norm, n);
    free_2D(matrix_wam);
    free_2D(matrix_ddg);
    return matrix_L_norm;
}

/*
 * Jacobi called from c,
 * Prints a matrix where The first line contains the eigenvalues,
 * and the second line onward contains the corresponding eigenvectors 
 */
void create_Jacobi(Point *points,int dim, int n){
    double **A, **eigen_vectors;
    double *eigen_values;
    A = matrix_init(n,dim);
    data_to_matrix(points, A, n, dim);
    eigen_vectors = jacobi(A, n);
    eigen_values = calloc(n, sizeof(double));
    Mem_Assertion(eigen_values != NULL);
    get_diag(A, eigen_values, n);
    print_Jacobi(eigen_vectors, eigen_values, n);
    free(eigen_values);
    free_2D(eigen_vectors);
    free_2D(A);
}

/*
 * Returns n and dim 
 */
Info extractInfo(FILE* file){
    char c;
    int dim, n;
    Info inf;
    dim = 1;
    n = 1;
    for (c = getc(file); c != '\n'; c = getc(file)){
        if (c == ','){
            dim++;
        }
    }
    for (c = getc(file); c != EOF; c = getc(file)){
        if (c == '\n'){
            n++;
        }
    }
    inf.dim = dim;
    inf.n = n;
    return inf;
}

/*
 * Stores data from file in points 
 */
int processFile(Info info, Point* points, FILE* file){
    int i,j;
    double coordinate; //changed defenition at start, was inside loop.
    for(i = 0; i < info.n; i++){
        for(j = 0; j < info.dim; j++){
            fscanf(file, "%lf", &coordinate);
            points[i].vector[j] = coordinate;
            getc(file);
        }
    }
    fclose(file);
    return 0;
}

/*
 * Runs code according to user input 
 */
void get_goal(char *goal, Point *points, int dim, int n){
    double **tmp;
    if(strcmp(goal, "wam") == 0) {
        tmp = create_WAM(points, dim, n);
        print_matrix(tmp, n, n);
        free_2D(tmp);
    }
    else if(strcmp(goal, "ddg") == 0) {
        tmp = create_DDG(points, dim, n);
        print_matrix(tmp, n, n);
        free_2D(tmp);
    }
    else if(strcmp(goal, "lnorm") == 0) {
        tmp = create_L_norm(points, dim, n);
        print_matrix(tmp, n, n);
        free_2D(tmp);
    }
    else if(strcmp(goal, "jacobi") == 0){
        create_Jacobi(points, dim, n);
    }
    else printf("Invalid Input!");
}

int main(int argc, char *argv[]) {
    char* goal;
    int dim, n;
    FILE* fp;
    Info info;
    Point* points;
    if (argc != 2) { //changed from 4 to 2: we get goal and file!
        printf("Invalid Input!\n");
        exit(1);
    }
    goal = argv[1];
    fp = fopen(argv[2], "r");
    if (fp == NULL){
        printf("Invalid Input!\n");
        fclose(fp);
        exit(1);
    }
    info = extractInfo(fp);
    dim = info.dim;
    n = info.n;
    points = allocate_mem(dim, n);
    rewind(fp); //added by Dana
    processFile(info, points, fp);
    get_goal(goal, points, dim, n);
    return 0;
}