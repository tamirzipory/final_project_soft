import math
import numpy as np
import pandas as pd
import sys
from enum import Enum
import spkmeans_module

class Goal(Enum):
    wam = 1
    ddg = 2
    lnorm = 3
    jacobi = 4
    spk = 5
    second_kmeans = 6


def print_output(mat_from_fit, len_of, goal):
    str_ret = ""
    num_rows = len_of

    if(goal == Goal.jacobi):
        num_rows += 1

    for i in range(num_rows):
        for j in range(len_of):
            str_ret += str('%.4f' % (mat_from_fit[i][j]))
            if (j != len_of-1):
                str_ret += ","
        if (i != num_rows-1):
            str_ret += "\n"
    print(str_ret)

def handle_errors():
    print("An Error Has Occurred")
    sys.exit()
    
def handle_errors_input():
    print("Invalid Input!")
    sys.exit()

def read_file(filename):
    try:
     data = pd.read_csv(filename, header=None)
     return data
    except:
        handle_errors()


def get_goal_input(filename):
    try:
        data = pd.read_csv(filename, header=None)
        return (data.to_numpy(), data.shape[0], data.shape[1])
    except:
        handle_errors()


def print_output_spk(mat_from_fit, K, arr2, list_cent_index):
    str_ret = ""
    for i in range(K):
        str_ret += str(int(list_cent_index[i]))
        if (i != K-1):
            str_ret += ","
    str_ret += "\n"

    
    for i in range(K):
        for j in range(arr2):
            str_ret += str('%.4f' % (mat_from_fit[i][j]))
            if (j != arr2-1):
                str_ret += ","
        if ( i != K-1):
            str_ret += "\n"
    print(str_ret)


def check_input(given_input, argc):
    if (argc != 4):
        handle_errors_input()
   
    is_valid = given_input[1].isnumeric()
    if is_valid:
        is_valid = int(given_input[1]) >= 0 
    
    if is_valid:
        is_valid = given_input[2] in [curr_goal.name for curr_goal in Goal]

    if not is_valid:
        handle_errors_input()

    return ((int)(given_input[1]), Goal[given_input[2]])


def start_kmeans(K, data_points_array):

    Centroids_array = []  
    Centroids_index_array = []

    len_of = len(data_points_array)
    arr_of_distance = np.array([0.0 for i in range(len_of+1)])
    arr_of_vectors = np.array([0.0 for i in range(len_of)])
    arr_of_index = np.array([i for i in range(len_of)])

    np.random.seed(0)

    index=np.random.choice(arr_of_index)
    Centroids_index_array.append(index)
    Centroids_array.append(data_points_array[index]) 

    for i in range(1, K):  
        find_D(arr_of_distance, data_points_array, len_of, Centroids_array)  
        
        arr_of_vectors = np.array([(arr_of_distance[l] / arr_of_distance[len_of]) for l in range(len_of)])
        index = np.random.choice(arr_of_index, p=arr_of_vectors)
        Centroids_index_array.append(index)
        Centroids_array.append(data_points_array[index])
    return Centroids_index_array,Centroids_array


def find_D(arr_of_distance, datapoints_array, len_of, Centroids_array):
    arr_of_distance[len_of] = 0.0
    for l in range(len_of):
       
        arr_of_distance[l] = np.min([calc(datapoints_array[l], centroid) for centroid in Centroids_array])
        
        arr_of_distance[len_of] = arr_of_distance[len_of] + arr_of_distance[l]


def calc(arr1, arr2):  
    sub = np.subtract(np.array(arr1), np.array(arr2))
    ret = np.sum(np.multiply(sub, sub))
    return ret


def call_fit_ex2(len_of, K, dimension, data_points, centroids, goal):
  
    try:
        final_centroids = spkmeans_module.fit(len_of, K, dimension, data_points, goal.value, centroids)
        return final_centroids
    except:
        handle_errors()

def print_last(matrix, i):
    print(",".join(["%.4f" % float(i) for i in matrix[i]]), end="")
    
def print_between(matrix, i):
    print(",".join(["%.4f" % float(i) for i in matrix[i]]))

def print_mat(matrix):
    for i in range(len(matrix)):
        if i == (len(matrix)-1):
            print_last(matrix, i)
        else:
            print_between(matrix, i)
def print_err():
    print("Invalid Input!")
        
def init_centroids(k, vector_list, vector_list_ind, vector_num):
    np.random.seed(0)
    if (not (k < vector_num)):
        print_err()
        assert(k < vector_num)
    dist = [0 for i in range(vector_num)]
    second_centroids = [0 for i in range(k)] 
    centroids_index = [0 for i in range(k)]
    rand_index = np.random.choice(vector_num)
    second_centroids[0] = vector_list[rand_index]
    centroids_index[0] = vector_list_ind[rand_index]
    z = 1
    while z < k:
        for i in range(vector_num):
            if z == 1:
                dist[i] = distance(vector_list[i], second_centroids[0])
            else:
                dist[i] = min_dist_centroid(
                    vector_list, i, second_centroids, z, dist)
        sum = np.sum(dist)
        prob = dist/sum
        chosen_ind = np.random.choice(vector_num, p=prob)
        second_centroids[z] = vector_list[int(chosen_ind)]
        centroids_index[z] = vector_list_ind[int(chosen_ind)]
        z += 1
    for i in range(k):
        second_centroids[i] = second_centroids[i].tolist()
    return second_centroids, centroids_index


def readfile(filename):
    file = pd.read_csv(filename, header=None)
    ret = file
    return ret

def distance(vector_1, vector_2):
    return np.sum((vector_1-vector_2)**2)


def min_dist_centroid(vector_list, i, second_centroids, z, dist):
    return min(distance(vector_list[i], second_centroids[z-1]), dist[i])


def fix_mat_to_zeros(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if ((matrix[i][j] > -0.00005) and (matrix[i][j] < 0)):
                matrix[i][j] = 0
    return matrix

def main(argv):

    inputs = argv
    inputs_len = len(inputs)
    K, goal = check_input(argv, inputs_len)
    data_points_array, len_of, the_g = get_goal_input(argv[3])

    if K > len_of:
        handle_errors_input()
    try:
        goal_matrix = spkmeans_module.fit(len_of, K, the_g, data_points_array.tolist(), goal.value, [])
        if(goal != Goal.spk):
            print_output(goal_matrix, len_of, goal)
        else:
            if(K == 0):
                K = len(goal_matrix[0])
            goal=Goal.second_kmeans
            centroids_index, centroids = start_kmeans(K, goal_matrix)
            the_g=K
            print_output_spk(call_fit_ex2(len_of, K, the_g,goal_matrix, centroids, goal), K, the_g, centroids_index)

    except:
        handle_errors()


if __name__ == '__main__':
    main(sys.argv)
