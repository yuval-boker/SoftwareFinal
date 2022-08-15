import sys
import numpy as np
import pandas as pd
import spkmeans_capi

MAX_ITER = 300
JACOBI_MAX_ITER = 100
np.random.seed(0)

# def print_matrix(mat):
#     np.set_printoptions(precision=4)
#     print(mat)
def double_to_str(num):
    """
    Print a double with 4 point persicion,
    includes a fix for situations where -0.0000 is printed
    """
    if -0.00005 < num < 0:
        num = 0.0
    return f'{num:.4f}'


def print_list(lst):
    for i, vec in enumerate(lst):
        print(','.join(map(double_to_str, vec)), end='' if i == len(lst)-1 else '\n')

def prepare():
    assert len(sys.argv) == 4, 'Invalid Input!'
    k = validity_check_k()
    goal = sys.argv[2]
    file = sys.argv[3]
    return k, goal, file


def validity_check_k():
    assert sys.argv[1].isnumeric(), 'Invalid Input!'
    k = int(sys.argv[1])
    return k


def process_file(file, k):
    data = pd.read_csv(file, header=None)
    n = data.shape[0]
    dim = data.shape[1]
    data_array = data.to_numpy().flatten().tolist()
    assert k < n, 'Invalid Input!'
    return data_array, n, dim


def get_goal(goal, data_array, n, dim):
    if goal == "wam":
        res = spkmeans_capi.get_WAM((data_array, n, dim))
        # res = np.array(res).reshape(n, n)
    elif goal == "ddg":
        res = spkmeans_capi.get_DDG((data_array, n, dim))
        # res = np.array(res).reshape(n, n)
    elif goal == "lnorm":
        res = spkmeans_capi.get_L_norm((data_array, n, dim))
        # res = np.array(res).reshape(n, n)
    elif goal == "jacobi":
        spkmeans_capi.run_jacobi((data_array, n, dim))
    print_list(res)

# add goal = spk (needs to get data matrix from c, and call kmeans_pp), goal = jacobi

# implementation of k-means++, n = number of points, dim = dimension of each point,
# init_centroids = array with the index of k points selected to be the initial centroids
# min_dis = at first is a constant of infinity, after first iteration will be an np array of size n containing
# the min distances of the points given from all of the centroids that have been chosen.
# p = the probability to be selected of each point in points, at start the probability of each point is equal.
# def k_means_pp(points, k):
#         n, dim = points.shape
#         init_centroids = []
#         min_dis = np.inf
#         p = None
#         for j in range(k):  # initializing k centroids as in k-means++ initialization
#             curr = np.random.choice(n, p=p)  # picking a random index of points provided
#             init_centroids.append(curr)
#             distances = np.power((points-points[curr]), 2).sum(axis=1)
#             min_dis = np.minimum(distances, min_dis)
#             p = np.divide(min_dis, min_dis.sum())
#         res = spkmeans_capi.fit(k, n, dim, max_iter, eps,  points[init_centroids].tolist(), points.tolist())
#         print(','.join([str(i) for i in pIndicies[init_centroids]]))  # prints the indices of observations chosen by
#         # the K-means++ algorithm as the initial centroids.
#         for centroid in res:  # prints the final centroids from the K-means algorithm executed in c
#             print(",".join('{:.4f}'.format(np.round(coord, 4)) for coord in centroid))

if __name__ == '__main__':
    k, goal, file = prepare()
    data_points, n, dim = process_file(file, k)
    get_goal(goal, data_points, n, dim)
