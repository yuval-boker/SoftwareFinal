import sys
import numpy as np
import pandas as pd
import spkmeans_capi

MAX_ITER = 300
JACOBI_MAX_ITER = 100
np.random.seed(0)

def prepare():
    assert len(sys.argv) == 3, 'Invalid Input!'
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
        res = spkmeans_capi.get_WAM(data_array, n, dim)
        res = np.array(res).reshape(n, n)
    elif goal == "ddg":
        res = spkmeans_capi.get_DDG(data_array, n, dim)
        res = np.array(res).reshape(n, n)
    elif goal == "lnorm":
        res = spkmeans_capi.get_L_norm(data_array, n, dim)
        res = np.array(res).reshape(n, n)
    elif goal == "jacobi":
        spkmeans_capi.run_jacobi(data_array, n, dim)

# add goal = spk (needs to get data matrix from c, and call kmeans_pp), goal = jacobi


def distance_calc(column, index, auxiliary, points):
    x = points - points.loc[index]
    auxiliary[column] = np.sum(x * x, axis=1)


def probability_calc(auxiliary):
    total = auxiliary["min_D"].sum()
    auxiliary["probability"] = auxiliary["min_D"]/total


def kmeans_pp():
    np.random.seed(0)
    c = np.random.choice(n - 1, 1)[0]
    distance_calc("min_D", c, auxiliary, points_matrix)
    centroids_by_index.append(c)
    centroids.append(points_list[c])
    z = 1
    while z < k:
        curr_miu = centroids_by_index[z-1]
        distance_calc("curr_D", curr_miu, auxiliary, points_matrix)
        auxiliary["min_D"] = np.minimum(auxiliary["min_D"], auxiliary["curr_D"])
        probability_calc(auxiliary)
        tmp = np.random.choice(n, 1, p=auxiliary["probability"])[0]
        centroids_by_index.append(tmp)
        centroids.append(points_list[tmp])
        z += 1


# def main():
#     k, goal, file = prepare()
#     data_points, n, dim = process_file(file, k)
#     data_array = data_points.to_numpy().flatten().tolist()
#     get_goal(goal, data_array, n, dim)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    k, goal, file = prepare()
    data_points, n, dim = process_file(file, k)
    data_array = data_points.to_numpy().flatten().tolist() #why again???
    get_goal(goal, data_array, n, dim)

