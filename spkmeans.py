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
    spk_ex2 = 6


def print_output(mat_from_fit, N, goal):
    output_res = ""
    num_rows = N

    if(goal == Goal.jacobi):
        num_rows += 1

    for i in range(num_rows):
        for j in range(N):
            output_res += str('%.4f' % (mat_from_fit[i][j]))
            if (j != N-1):
                output_res += ","
        if (i != num_rows-1):
            output_res += "\n"
    print(output_res)

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


def print_output_spk(mat_from_fit, K, D, centroids_index_list):
    output_res = ""
    for i in range(K):
        output_res += str(int(centroids_index_list[i]))
        if (i != K-1):
            output_res += ","
    output_res += "\n"

    
    for i in range(K):
        for j in range(D):
            output_res += str('%.4f' % (mat_from_fit[i][j]))
            if (j != D-1):
                output_res += ","
        if ( i != K-1):
            output_res += "\n"
    print(output_res)


def check_input(given_input, argLen):
    if (argLen != 4):
        handle_errors_input()
   
    is_valid = given_input[1].isnumeric()
    if is_valid:
        is_valid = int(given_input[1]) >= 0 
    
    if is_valid:
        is_valid = given_input[2] in [curr_goal.name for curr_goal in Goal]

    if not is_valid:
        handle_errors_input()

    return ((int)(given_input[1]), Goal[given_input[2]])


def kMeans_init(K, data_points_array):

    Centroids_array = []  
    Centroids_index_array = []

    N = len(data_points_array)
    D_array = np.array([0.0 for i in range(N+1)])
    Pr_array = np.array([0.0 for i in range(N)])
    Index_array = np.array([i for i in range(N)])

    np.random.seed(0)

    index=np.random.choice(Index_array)
    Centroids_index_array.append(index)
    Centroids_array.append(data_points_array[index]) 

    for i in range(1, K):  
        find_D(D_array, data_points_array, N, Centroids_array)  
        
        Pr_array = np.array([(D_array[l] / D_array[N]) for l in range(N)])
        index = np.random.choice(Index_array, p=Pr_array)
        Centroids_index_array.append(index)
        Centroids_array.append(data_points_array[index])
    return Centroids_index_array,Centroids_array


def find_D(D_array, datapoints_array, N, Centroids_array):
    D_array[N] = 0.0
    for l in range(N):
       
        D_array[l] = np.min([calc(datapoints_array[l], centroid) for centroid in Centroids_array])
        
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
            goal=Goal.spk_ex2
            centroids_index, centroids = kMeans_init(K, goal_matrix)
            D=K
            print_output_spk(call_fit_ex2(N, K, D,goal_matrix, centroids, goal), K, D, centroids_index)

    except:
        handle_errors()


if __name__ == '__main__':
    main(sys.argv)
