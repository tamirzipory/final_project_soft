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

def print_output(mat_get_c, N, goal):
    ret_out = ""
    num_rows = N
    if(goal == Goal.jacobi):
        num_rows += 1

    for i in range(num_rows):
        for j in range(N):
            ret_out += str('%.4f' % (mat_get_c[i][j]))
            if (j != N-1):
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

    N = len(data_points_array)
    D_array = np.array([0.0 for i in range(N+1)])
    p_arr = np.array([0.0 for i in range(N)])
    arr_of_index = np.array([i for i in range(N)])

    np.random.seed(0)

    index=np.random.choice(arr_of_index)
    the_op_arr.append(index)
    arr_of_cent_index.append(data_points_array[index]) 

    for i in range(1, K):  
        find_D(D_array, data_points_array, N, arr_of_cent_index)  
        
        p_arr = np.array([(D_array[l] / D_array[N]) for l in range(N)])
        index = np.random.choice(arr_of_index, p=p_arr)
        the_op_arr.append(index)
        arr_of_cent_index.append(data_points_array[index])
    return the_op_arr,arr_of_cent_index


def find_D(D_array, datapoints_array, N, arr_of_cent_index):
    D_array[N] = 0.0
    for l in range(N):
       
        D_array[l] = np.min([calc(datapoints_array[l], centroid) for centroid in arr_of_cent_index])
        
        D_array[N] = D_array[N] + D_array[l]


def calc(x, y):  
    z = np.subtract(np.array(x), np.array(y))
    return np.sum(np.multiply(z, z))


def call_fit_ex2(N, K, dimension, data_points, centroids, goal):
  
    try:
        final_centroids = spkmeans_module.fit(N, K, dimension, data_points, goal.value, centroids)
        return final_centroids
    except:
        handle_errors()



def main(argv):

    inputs = argv
    inputs_len = len(inputs)
    K, goal = check_input(argv, inputs_len)
    data_points_array, N, D = get_goal_input(argv[3])

    if K > N:
        handle_errors_input()
    try:
        goal_matrix = spkmeans_module.fit(N, K, D, data_points_array.tolist(), goal.value, [])
        if(goal != Goal.spk):
            print_output(goal_matrix, N, goal)
        else:
            if(K == 0):
                K = len(goal_matrix[0])
            goal=Goal.kmeans_enum
            centroids_index, centroids = kMeans_init(K, goal_matrix)
            D=K
            display_out_of_spk(call_fit_ex2(N, K, D,goal_matrix, centroids, goal), K, D, centroids_index)

    except:
        handle_errors()


if __name__ == '__main__':
    main(sys.argv)
