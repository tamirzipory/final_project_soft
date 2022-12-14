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
    kmeans_enum = 6

def print_output(mat_get_c, max_iter, goal):
    ret_out = ""
    num_rows = max_iter
    if(goal == Goal.jacobi):
        num_rows += 1

    for i in range(num_rows):
        for j in range(max_iter):
            ret_out += str('%.4f' % (mat_get_c[i][j]))
            if (j != max_iter-1):
                ret_out += ","
        if (i != num_rows-1):
            ret_out += "\n"
    print(ret_out)

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


def display_out_of_spk(mat_get_c, K, d_mat, centroids_index_list):
    ret_out = ""
    for i in range(K):
        ret_out += str(int(centroids_index_list[i]))
        if (i != K-1):
            ret_out += ","
    ret_out += "\n"

    
    for i in range(K):
        for j in range(d_mat):
            ret_out += str('%.4f' % (mat_get_c[i][j]))
            if (j != d_mat-1):
                ret_out += ","
        if ( i != K-1):
            ret_out += "\n"
    print(ret_out)


def check_input(argv, argc):
    len_of_input = argc
    if (len_of_input != 4):
        handle_errors_input()
   
    boolean_validation = argv[1].isnumeric()
    if boolean_validation:
        boolean_validation = int(argv[1]) >= 0 
    
    if boolean_validation:
        boolean_validation = argv[2] in [curr_goal.name for curr_goal in Goal]

    if not boolean_validation:
        handle_errors_input()

    return ((int)(argv[1]), Goal[argv[2]])


def kMeans_init(K, data_points_array):

    arr_of_cent_index = []  
    the_op_arr = []

    max_iter = len(data_points_array)
    ret_d_array = np.array([0.0 for i in range(max_iter+1)])
    p_arr = np.array([0.0 for i in range(max_iter)])
    arr_of_index = np.array([i for i in range(max_iter)])

    np.random.seed(0)

    index=np.random.choice(arr_of_index)
    the_op_arr.append(index)
    arr_of_cent_index.append(data_points_array[index]) 

    for i in range(1, K):  
        calc_d_array(ret_d_array, data_points_array, max_iter, arr_of_cent_index)  
        
        p_arr = np.array([(ret_d_array[l] / ret_d_array[max_iter]) for l in range(max_iter)])
        index = np.random.choice(arr_of_index, p=p_arr)
        the_op_arr.append(index)
        arr_of_cent_index.append(data_points_array[index])
    return the_op_arr, arr_of_cent_index


def calc_d_array(ret_d_array, datapoints_array, max_iter, arr_of_cent_index):
    ret_d_array[max_iter] = 0.0
    for l in range(max_iter):
       
        ret_d_array[l] = np.min([calc(datapoints_array[l], centroid) for centroid in arr_of_cent_index])
        
        ret_d_array[max_iter] = ret_d_array[max_iter] + ret_d_array[l]


def calc(x, y):  
    z = np.subtract(np.array(x), np.array(y))
    return np.sum(np.multiply(z, z))


def triget_method2(max_iter, K, dimension, data_points, centroids, goal):
  
    try:
        final_centroids = spkmeans_module.fit(max_iter, K, dimension, data_points, goal.value, centroids)
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
    data_points_array, max_iter, d_arr = get_goal_input(argv[3])

    if K > max_iter:
        handle_errors_input()
    try:
        goal_matrix = spkmeans_module.fit(max_iter, K, d_arr, data_points_array.tolist(), goal.value, [])
        if(goal != Goal.spk):
            print_output(goal_matrix, max_iter, goal)
        else:
            if(K == 0):
                K = len(goal_matrix[0])
            goal=Goal.kmeans_enum
            centroids_index, centroids = kMeans_init(K, goal_matrix)
            d_arr=K
            display_out_of_spk(triget_method2(max_iter, K, d_arr,goal_matrix, centroids, goal), K, d_arr, centroids_index)

    except:
        handle_errors()


if __name__ == '__main__':
    main(sys.argv)
