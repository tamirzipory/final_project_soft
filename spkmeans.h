import sys
import numpy as np
import numpy.random
import spkmeansmodule as kp

np.random.seed(0)


def euclidean_dist(vector1, vector2):
    """
    Computes the euclidean distance between two given vectors.
    """
    dist = (vector1 - vector2) ** 2
    return np.sqrt(sum(dist))


def print_result(result, index, goal, n):
    """
    Prints the result according to the goal.
    :param index: case goal == spk: an array of the chosen indices of the observations
                                    chosen as the initial centroids.
                  case goal == jacobi: an array of the eigenvalues of the matrix.
                  other cases: None.
    :param n: the row dimension of the result matrix.
    """
    if goal == "spk":
        print(*index, sep=",")
    elif goal == "jacobi":
        index = [format(round(index[i], 4) + 0, '.4f') for i in range(len(index))]
        print(*index, sep=",")

    # formatting the output as the assignment requirements
    result = ['%.4f' % result[i] for i in range(len(result))]
    result = [result[i:i + n] for i in range(0, len(result), n)]
    for i in range(len(result)):
        print(*(result[i]), sep=",")


def kmeans_pp(vectors, k):
    """
    Implantation of the k-means++ algorithm, as detailed in HW2.
    :return: clusters - 2D matrix of the observations chosen as the initial centroids.
             index - the indices of the chosen observations as initial centroids.
    """
    index = []  # contains the indices of the chosen observations

    # choosing a random observation as the first centroid
    clusters = np.zeros((k, len(vectors[0])), dtype="float")
    first_cluster_index = np.random.choice(len(vectors))
    index.append(first_cluster_index)
    clusters[0] = np.copy(vectors[first_cluster_index])

    # randomly choosing observations with a weighted probability,
    # favoring a distant observation from the centroids we already chose
    distances = np.full(len(vectors), sys.float_info.max)
    i = 1
    while i < k:
        for t in range(len(vectors)):
            for j in range(i):
                dist = euclidean_dist(vectors[t], clusters[j])
                if dist < distances[t]:
                    distances[t] = dist
        probabilities = distances / sum(distances)  # compute probabilities
        chosen_index = int(np.random.choice(len(vectors), 1, True, probabilities))
        clusters[i] = np.copy(vectors[chosen_index])
        index.append(chosen_index)
        i += 1
    return clusters, index


def k_legal(k):
    """
    Validation of the k given via line argument is in correct format.
    """
    if k - int(k) > 0 or (k != 0 and k <= 1):
        print("Invalid Input!")
        return False
    return True


def legal_goal(goal):
    """
    Validation of the goal given via line argument is in correct format.
    """
    goal_list = ["wam", "ddg", "lnorm", "jacobi", "spk"]
    return goal in goal_list


def execute_by_goal(vectors, goal):
    """
    Execution of goals wam, ddg, lnorm via calling the C extension
    """
    goal_to_function = {
        "wam": kp.wam,
        "ddg": kp.ddg,
        "lnorm": kp.lnorm,
    }
    return goal_to_function[goal](len(vectors), len(vectors[0]), vectors.flatten().tolist())


def create_vectors_from_file(file):
    """
    Form a 2D matrix containing the n given points in the file as rows
    """
    try:
        vectors = np.genfromtxt(file, delimiter=",")
        return vectors
    except IOError:
        print("Invalid Input!")
        return None


def goal_is_spk(vectors, k):
    """
    Implantation of the spk algorithm via calling the C extension.
    :param k: case k == 0 - first determine k via the Eigengap Heuristic,
    :return: clusters - 2D matrix of the k final clusters
            chosen_index - the indices of the chosen observations as initial centroids.
            n - the row dimension of clusters.
    """
    n = len(vectors[0])
    if k == 0:
        k, vectors = kp.spk(len(vectors), len(vectors[0]), vectors.flatten().tolist())
        if vectors is None or k <= 1:
            return None, None, None
        vectors = np.array(vectors).reshape(int(len(vectors) / k), k)
        n = int(k)
    k = int(k)
    init_clusters, chosen_index = kmeans_pp(vectors, k)
    # implement the k-means algorithm
    clusters = kp.fit(k, len(vectors), len(vectors[0]), vectors.flatten().tolist(), init_clusters.flatten().tolist())
    return clusters, chosen_index, n


def main():
    # Reading user CMD arguments and validation that are in correct format.
    cmd_input = sys.argv
    if len(cmd_input) != 4:
        print("Invalid Input!")
        return 1
    k = float(cmd_input[1])
    goal = cmd_input[2]
    file_name = cmd_input[3]
    if not legal_goal(goal):
        print("Invalid Input!")
        return 1

    vectors = create_vectors_from_file(file_name)
    if vectors is None:
        return 1

    # in case of 'spk' will contain the chosen indices of the observations points.
    # in case of 'jacobi' will contain the eigenvalues of the matrix.
    # other cases will stay None.
    chosen_index = None

    # the row dimension of the result matrix. may change in case of 'spk'.
    n = len(vectors)

    # Execution by goal.
    if goal == "spk":
        if not k_legal(k):
            return 1
        result, chosen_index, n = goal_is_spk(vectors, k)
    elif goal == "jacobi":
        result, chosen_index = kp.jacobi(n, len(vectors[0]), vectors.flatten().tolist())
    else:
        result = execute_by_goal(vectors, goal)

    if result is None:
        print("An Error Has Occurred", end ="")
        return 1

    print_result(result, chosen_index, goal, n)
    return 0


if __name__ == '__main__':
    main()
