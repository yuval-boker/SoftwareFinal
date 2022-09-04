#include "spkmeans.h"


/** print funcs **/
void print_double(double num) {
    if (num < 0 && num > -0.00005)
        num = 0.0;
    printf("%.4f", num);
}

void print_row(double *row, int n){
    int i;
    for (i = 0;i < n - 1; i++){
        printf("%.4f,",row[i]);
    }
    printf("%.4f\n", row[n-1]);
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
    int i,j;
    for (i = 0; i < n; i++){
        for (j = 0;j < dim - 1; j++){
            printf("%.4f,", array[i][j]);
        }
        printf("%.4f\n", array[i][dim-1]);
    }
}
/* for test use
void print_transpose(double **mat, int n, int dim){
    int i,j;
    for (j = 0; j<dim;j++){
        for(i = 0; i<n; i++){
            print_double(mat[i][j]);
            if(i<n-1){
                printf("%c",',');
            }
        }
        printf("\n");
    }
}
 */
/*
 * Prints the eigenvalues as the first line, second line onward is the corresponding eigenvectors
 */
void print_Jacobi(double **eigen_vectors, double *eigen_values, int n) { 
    print_row(eigen_values, n);
    printf("\n");
    print_matrix(eigen_vectors, n, n);
    /*print_transpose(eigen_vectors,n,n);*/
}

/*
 * Allocates memory for n data points stored as matrix of Point 
 * input: dimension and size
 * output: an empty matrix of size n*m
 */
Point* allocate_mem(int dim, int n) {
    int i,j;
    double* v;
    Point* points = (Point*)calloc(n, sizeof(Point));
    Mem_Assertion(points != NULL);
    for (i = 0; i < n; i++){
        v = (double*) calloc(dim,sizeof(double));
        /*Mem_Assertion(v != NULL);*/
        if (!v){
            for (j = 0; j < i; j++){
                free(points[j].vector);
            }
            free(points);
            Mem_Assertion(v);
        }
        points[i].vector = v;
    }
    return points;
}

/*
 * Frees points memory, used before returning the output to Python
 */
void free_data_points(int n, Point* points) {
    int i;
    if (points != NULL){
        for (i = 0; i < n; i++){
            free(points[i].vector);
        }
        free(points);
    }
}

/*
 * Frees clusters memory, used before returning the output to Python
 */
void free_clusters(int k, Cluster* clusters) {
    int j;
    for (j = 0; j < k; j++){
        free(clusters[j].newPoints);
        free(clusters[j].centroid);
    }
    free(clusters);
}

void free_memory(int k, int n, Point* points, Cluster* clusters) {
    free_data_points(n, points);
    free_clusters(k, clusters);
}

/*
 * Deallocates the memory that was previously dynamically allocated 
 */
void free_2D(double **matrix) {
    if(matrix){
        free(matrix[0]);
        free(matrix);
    }
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
double sum_row(double **matrix, int row, int n) {
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
void data_to_matrix(Point *points, double **matrix, int n, int dim) {
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
int find_max_indices_off_diag(double **mat, int *i_val, int *j_val, int n) {
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
int transform_matrix(double **mat, double **v, int n, int i, int j) {
    int r;
    double tmp1, tmp2, theta, t, c, s;
    theta = (mat[j][j] - mat[i][i])/(2 * mat[i][j]);
    t = sign(theta)/(fabs(theta) + sqrt((theta * theta)+1.0));
    c = 1.0/(sqrt((t * t) + 1.0));
    s = t * c;
    /*printf("theta= %lf t= %lf c= %lf s= %lf\n",theta,t,c,s);*/
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
            sum += pow(mat[i][j],2);
        }
    }
    return 2 * sum;
}
/*
 * creates jacobi
 */
double **jacobi(double **A, int n) {
    /*printf("n = %d",n);*/
    double **V;
    double off_A = 0.0, off_A_prime;
    int i, j, l = 0;
    V = I_matrix(n); /* initializing V to Identity matrix */
    /*printf("V start:\n");
    print_matrix(V,n,n);*/
    while (l < MAX_ITER_J) {
        if(l==0){
            off_A = off_square(A, n);
        }
        find_max_indices_off_diag(A, &i, &j, n);
        /*printf("i is %d, j is %d \n",i,j);*/
        transform_matrix(A, V, n, i, j);
        /*printf(" A after %d iteration is:\n",l);
        print_matrix(A,n,5);
        printf("V after %d iteration is:\n",l);
        print_matrix(V,n,n);*/
        off_A_prime = off_square(A, n);
        if((l != 0) && (off_A - off_A_prime <= EPSILON_J)) {
            /*printf("off is:%lf\n",off_A-off_A_prime);*/
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
void get_diag(double **mat, double *diag, int n) {
    int i;
    for(i = 0; i < n; i++){
        diag[i] = mat[i][i];
    }
}

/** EIGENGAP HEURISTIC **/
/*
* Given a EigenData array sorts inplace by eigenvalues in decreasing order.
* we used this reference:
* https://www.javatpoint.com/c-program-to-sort-the-elements-of-an-array-in-descending-order
*/
void sort(EigenData eigen_arr[], int n) {
    int i, j;
    EigenData temp;
    for(i = 0; i < n; i++){
        for(j = i + 1; j < n; j++){
            if(eigen_arr[i].val < eigen_arr[j].val) {
                temp = eigen_arr[i];
                eigen_arr[i] = eigen_arr[j];
                eigen_arr[j] = temp;
            }
        }
    }
}

/*
 * Given an n eigenvectors 2D array and a 1D array of n corresponding eigenvalues,
 * creates a new 1D EigenData (contains eigenvalue it's eigenvector) array
 * which is decreasingly ordered by eigenvalues and returns the new array.
 */
EigenData *create_sorted_eigen_arr(double **eigen_vectors, double *eigen_values, int n) {
    int i;
    EigenData* eigen_arr =(EigenData*)calloc(n, sizeof(EigenData));
    Mem_Assertion(eigen_arr != NULL);
    for(i = 0; i < n; i++){
        eigen_arr[i].val = eigen_values[i];
        eigen_arr[i].vector = eigen_vectors[i];
    }
    sort(eigen_arr, n);
    return eigen_arr;
}

/*
 * Given an array of n EigenData elements (each element contains eigenvalue it's eigenvector)
 * implements eigengap heuristic as described in 1.3 and returns k.
 */
int eigengap(EigenData *eigen_arr, int n) {
    int i, max_i = 0;
    double curr_delta, max_delta = 0.0;
    for(i = 0; i < (n / 2); i++){
        curr_delta = fabs(eigen_arr[i].val - eigen_arr[i + 1].val);
        if(curr_delta > max_delta){
            max_delta = curr_delta;
            max_i = i + 1; /* i starts from 1 in eigengap */
        }
    }
    return max_i;
}

/** SPK **/
/*
 * Creates U, an nxk matrix containing the first k columns of sorted Lnorm,
 * Returns U
 */
double **create_U(EigenData *eigen_arr, int n, int k) {
    int i, j;
    double **U;
    U = matrix_init(n, k);
    for(j = 0; j < k; j++){
        for(i = 0; i < n; i++){
            U[i][j] = eigen_arr[j].vector[i];
        }
    }
    return U;
}

/*
 * Given U, this function renormalizes each of U’s rows to have unit length
 * and thus creates T as described in step 5 of The Normalized Spectral Clustering Algorithm
 */
void normalize_U(double **U, int n, int k) {
    int i, j;
    double sum;
    for(i = 0; i < n; i++){
        sum = 0.0;
        for(j = 0; j < k; j++){
            sum += pow(U[i][j], 2);
        }
        sum = sqrt(sum);
        if (sum != 0){
            for(j = 0; j < k; j++){
                U[i][j] /= sum;
            }
        }
    }
}
/** KMEANS **/
/*
 * Computes the distance between the given point and the given centroid
 * input: point, centroid, dimension
 * output: distance, double  
 */
double distance(Point* point1, const double* centroid, int dim){
    int i;
    double sum, p;
    sum = 0;
    for (i = 0; i < dim; i++){
        p = point1->vector[i] - *(centroid + i);
        sum += p*p;
    }
    return sqrt(sum);
}

/*
 * Finds the closest cluster to the given point
 * input: point, list of clusters, dimension, k
 * output: index of the required cluster, int  
 */
int min_distance(Point* point, Cluster* clusters, int dim, int k){
    int min_index, i;
    double min_val, curr;
    min_index = 0;
    min_val = distance(point, clusters[0].centroid, dim); /* define the distance from the first cluster as min_val */
    for (i = 1; i < k; i++){
        curr = distance(point, clusters[i].centroid, dim);
        if (curr < min_val){
            min_val = curr;
            min_index = i;
        }
    }
    return min_index;
}

/*
 * Compute the euclidean norm for KMEANS
 */
double euclidean_norm(const double* vector, int dim){
    int i;
    double result,p;
    result = 0;
    for (i = 0; i < dim; i++){
        p = vector[i];
        result += p*p;
    }
    return sqrt(result);
}

/*
 * Add the given point to the given cluster 
 * The insert operation is implemented by adding each coordinate to its adequate index in new_points vector
 * input: point, cluster, dimension
 */
void add_point(Point* point, Cluster* cluster, int dim){
    int i;
    for (i = 0; i < dim; i++) {
        cluster->newPoints[i] += point->vector[i];
    }
    cluster->count += 1; /* increase the number of points in the given cluster by 1 */ 
}
/*
 * Update the centroid of the given cluster by computing the average value of each coordinate in new_points
 * Check convergence
 * input: cluster, dimension, an empty vector used for holding the result of the computation
 * output: indicator for convergence, int
 */
int centroid_update(Cluster* cluster, int dim, double *tmp_vector){
    int has_changed, i, l;
    double norm_check;
    has_changed = 1;
    if (cluster->count == 0){
        return 1;
    }
    for (i = 0; i < dim; i++) {
        tmp_vector[i] = cluster->newPoints[i]/cluster->count;
    }
    norm_check = euclidean_norm(cluster->centroid, dim) - euclidean_norm(tmp_vector, dim);
    if (norm_check >= EPS || norm_check <= -EPS){ /*?*/
        has_changed = 0;
    }
    for (l = 0; l < dim; l++) {
        cluster->centroid[l] = tmp_vector[l];
        cluster->newPoints[l] = 0;
    }
    cluster-> count = 0;
    return has_changed;
}

/*
 * Update the centroids of the clusters, and initialize count and new_points for each cluster
 * Called at the end of each iteration in Kmeans
 * input: list of clusters, k, dimension
 * output: indicator for convergence, int
 */
int clusters_update(Cluster* clusters, int k, int dim) {
    int changed, i, epsilon_indicator;
    double *tmp_vector;
    changed = 1;
    tmp_vector = (double *) calloc(dim, sizeof(double)); /* Create a temporary vector to hold the new values before assigning them to the centroid */
    Mem_Assertion(tmp_vector != NULL);
    for (i = 0; i < k; i++) {
        epsilon_indicator = centroid_update(&clusters[i], dim, tmp_vector);
        changed = ((changed) && (epsilon_indicator));
    }
    free(tmp_vector); /* free 'tmp_vector' */
    return changed;
}

/*
 * kmeans algorithem implementaion
 * input: n, clusters = list of clusters, points = list of data points, dim = dimension, k
 */
void kmeans(int n, Cluster* clusters, Point* points, int dim, int k) {
    int epsilon_check, iter, i, index;
    epsilon_check = 0; 
    /* epsilon_check is an identicator that has a value of 1 iff the euclidean norm of each centroids doesn't change by more then epsilon. */
    iter = 0; /* iterations counter */
    while ((iter < MAX_ITER) && (1 - epsilon_check)) {
        for (i = 0; i < n; i++) {
            index = min_distance(&points[i], clusters, dim, k);
            add_point(&points[i], &clusters[index], dim);
        }
        epsilon_check = clusters_update(clusters, k, dim);
        iter++;
    }
}
/*
 * Creates WAM using points 
 */
double **create_WAM(Point *points, int dim, int n) {
    double **matrix;
    matrix = matrix_init(n, n);             /* WAM ∈ R^(n×n) */
    set_WAM(points, matrix, dim, n);
    return matrix;
}

/*
 * Sets DDG using WAM
 */
void set_DDG(double **WAM, double **matrix, int n) {
    int i;
    for (i = 0; i < n; i++){
        matrix[i][i] = sum_row(WAM, i, n);
    }
}

/*
 * Creates DDG using WAM creation function 
 */
double **create_DDG(Point *points, int dim, int n) {
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
void set_L_norm(double **WAM, double **DDG, double **L_norm, int n) {
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
void create_Jacobi(Point *points, int dim, int n) {
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
 * This function is used only for "spk" goal from Python.
 * Creates T matrix according to the instructions presented in the project.
 * Calculates and sets k using eigengap heuristic if given k value is zero.
 */
double **create_T(Point *points, int dim, int n, int *k) {
    double **matrix_L_norm, **U, **eigen_vectors;
    double *eigen_values;
    EigenData *eigen_arr;
    matrix_L_norm = create_L_norm(points, dim, n);
    eigen_vectors = jacobi(matrix_L_norm, n);
    eigen_values = calloc(n, sizeof(double));
    Mem_Assertion(eigen_values != NULL);
    get_diag(matrix_L_norm, eigen_values, n);
    eigen_arr = create_sorted_eigen_arr(eigen_vectors, eigen_values, n);
    if(*k == 0){
        *k = eigengap(eigen_arr, n);
    }
    U = create_U(eigen_arr, n, *k);
    normalize_U(U, n, *k);
    free(eigen_values);
    return U;
}

/*
 * Returns n and dim 
 */
Info extractInfo(FILE* file) {
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
int processFile(Info info, Point* points, FILE* file) {
    int i,j;
    double coordinate; 
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
void get_goal(char *goal, Point *points, int dim, int n) {
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
    if (argc != 3) { 
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
    rewind(fp); 
    processFile(info, points, fp);
    get_goal(goal, points, dim, n);
    return 0;
}
