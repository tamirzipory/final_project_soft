import sys
from enum import Enum
import pandas as pd
import numpy as np
from mykmeanssp import *

# kmeans++ random seed
np.random.seed(0)

# goal
class Goal(Enum):
    SPK = "spk"
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"

# constants
MAX_ITER = 300
EPSILON = 0

def check_k(k):
    """
    checks if the parameter k is valid and converts it to the right type.
    :param k
    :return: the value of the parameter, None in case of an error
    """
    if not k.isnumeric(): # checks if the string contains only digits
        return None
    k = int(k) # k==0 is valid
    return k

def check_goal(goal):
    """
    checks if the parameter goal is valid and converts it to the right type.
    :param goal
    :return: the value of the parameter, None in case of an error
    """
    try:
        goal = Goal(goal)
    except:
        return None
    return goal

def get_parameters():
    """
    gets the parameters from cmd and checks if they are valid.
    :return: tuple containing the parameters.
    """
    num_params = len(sys.argv) - 1 # without program name
    k = 0
    eps = EPSILON
    max_iter = MAX_ITER
    input_file_name = ""
    goal = None

    # get parameters from cmd
    if num_params != 3:
        invalid_input()
    else:
        k = sys.argv[1]
        goal = sys.argv[2]
        input_file_name = sys.argv[3]

    # convert parameters to the right types and check validity
    k = check_k(k)
    goal = check_goal(goal)
    parameters = (k, max_iter, eps, goal, input_file_name)
    if None in parameters: # at least one of the parameters is not valid
        invalid_input()
    return parameters

def read_input_file(input_file_name):
    """
    gets the input file and returns the data points as a pandas Data Frame.
    the index of the data points is from 0 to (n-1).
    :param input_file_name: string that represents the input file name.
    :return: data_points_df, a pandas Data Frame representing the data points.
    """
    data_points_df = pd.read_csv(input_file_name, header=None) # will not fail, beacuse file is assumes to be in right foramt
    return data_points_df

def kmeans_pp_centroids(data_points_df, k):
    """
    this function chooses the random initial centroids from the data points according to the given k.
    :param data_points_df: the output of "read_input_file"
    :param k: the desired number of clusters
    :return: a tuple containing the initial centroids and the initial centroids indexes.
    """
    initial_centroids = []
    initial_centroids_indexes = []
    # choose first centroid
    mu_index = np.random.choice(data_points_df.index)
    mu = list(data_points_df.iloc[mu_index])
    initial_centroids.append(mu)
    initial_centroids_indexes.append(mu_index)
    # choose the rest of the centroids
    for i in range(1, k):
        # calculate Dls
        Dl_results = data_points_df.apply(lambda data_point: min_euclidean_norm_squared(data_point, initial_centroids), axis=1).to_numpy()  # axis 1 - apply on rows
        # calculate probabilities
        P_results = Dl_results / np.sum(Dl_results) # sum can not be 0: k > 1, N > 1
        # choose next centroid
        mu_index = np.random.choice(data_points_df.index, p=P_results)
        mu = list(data_points_df.iloc[mu_index])
        initial_centroids.append(mu)
        initial_centroids_indexes.append(mu_index)

    return initial_centroids, initial_centroids_indexes

def euclidean_norm_squared(vec1, vec2):
    """
    gets two vectors (meaning points) and calculates their euclidean norm. returns the squared result.
    :param vec1: point 1.
    :param vec2: point 2.
    :return: (euclidean norm of point 1 and point 2)^2
    """
    return sum([(vec1[i] - vec2[i])**2 for i in range(len(vec1))])

def min_euclidean_norm_squared(vec, vec_list):
    """
    returns the minimal value of euclidean_norm_squared(vec,x) from all vectors x in vec_list.
    :param vec: some vector.
    :param vec_list: list of vectors of the same dimensions.
    :return: min euclidean_norm_squared.
    """
    return min([euclidean_norm_squared(vec, x) for x in vec_list])

def print_matrix(matrix, format_zero = False):
    """
    gets a matrix (meaning, a 2 dimensional list) and prints it in the desired format.
    :param matrix: a 2 dimensional floats list.
    :param format_zero: boolean flag, removes the "-" from "-0.0000" if on.
    :return: None.
    """
    for row in matrix:
        for i in range(len(row)):
            if i < len(row) - 1:
                out_str = "%.4f," % (row[i])
            else: # last coordinate in row without comma, and with new line
                out_str = "%.4f\n" % (row[i])
            if format_zero:  # if the format_zero flag is on, remove the "-" from zeros
                out_str = out_str.replace("-0.0000", "0.0000")
            print(out_str, end="")

def print_final_centroids(initial_centroids_indexes, final_centroids):
    """
    gets the final centroids from the kmeans algorithm output and the initial chosen centroids.
    prints the data in the desired format.
    :param initial_centroids_indexes: the chosen centroids indexes from kmeans_pp_centroids.
    :param final_centroids: the final centroids after calling the c module.
    :return: None.
    """
    # print initial centroid indexes
    print(",".join([str(x) for x in initial_centroids_indexes]))
    # print final centroids
    print_matrix(final_centroids)

# error messages
def invalid_input():
    print("Invalid Input!", end="")
    sys.exit()

def error():
    print("An Error Has Occurred", end="")
    sys.exit()

# main
def main():
    k, max_iter, eps, goal, input_file_name = get_parameters()
    data_points_df = read_input_file(input_file_name)
    N = data_points_df.shape[0]
    d = data_points_df.shape[1]
    data_points = data_points_df.values.tolist()  # convert data points to 2D list

    # kmeans option
    if goal == Goal.SPK:
        if k == 1 or k >= N: # check k for SPK option
            invalid_input()
        try:
            t_matrix = np.array(fit_nsc(data_points, N, d, k)) # calculate t matrix and find k if k=0
        except:
            error()
        # change the k and d according to the t matrix that was calculated
        k = t_matrix.shape[1]
        d = k
        # the data points are the rows of t matrix
        data_points_df = pd.DataFrame(t_matrix)
        data_points = data_points_df.values.tolist()

        # find the initial centroids for the kmeans algorithm
        initial_centroids, initial_centroids_indexes = kmeans_pp_centroids(data_points_df, k)
        try:
            final_centroids = np.array(fit_kmeanspp(data_points, initial_centroids, N, d, k, max_iter, eps))
        except:
            error()
        print_final_centroids(initial_centroids_indexes, final_centroids)
            
    # jacobi option
    elif goal == Goal.JACOBI:
        if N != d: # check that data_points is a square matrix
            invalid_input()
        try:
            out = calculate_jacobi(data_points, N, d) # in this case data_points is assumed to be a symmetric matrix
        except:
            error()
        print_matrix(out[0], format_zero=True) # prints the eigenvalues without the "-" before zeros
        print_matrix(out[1])

    # matrix options
    # goal in [Goal.WAM, goal.DDG, goal.LNORM]
    else:
        try:
            out = np.array(calculate_matrix(data_points, N, d, goal.value))
        except:
            error()
        print_matrix(out)

if __name__ == "__main__":
    main()
