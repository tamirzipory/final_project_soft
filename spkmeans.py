import sys
import pandas as pd
import numpy as np
import spkmeans_c as spk
import kmeans_c

# from sklearn.cluster import kmeans_plusplus, KMeans

goal_enum = ['spk', 'wam', 'ddg', 'lnorm', 'jacobi']


def read_from_file(file_name):
    data = pd.read_csv(file_name, header=None)
    return data


def read_params():
    k, goal, file_name = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    observations = read_from_file(file_name)

    if not (goal in goal_enum and 0 <= k <= len(observations) and k != 1):
        raise Exception("Invalid Input!")

    return k, goal_enum.index(goal), file_name, observations


def extract_inputs(file_name):
    data = pd.read_csv(file_name, sep=',', header=None)
    return data.to_numpy()


def kmeans_pp(k, datapoints):
    # centers, indices = kmeans_plusplus(datapoints, n_clusters=k, random_state=0, n_local_trials=1)
    # kmeans = KMeans(n_clusters=k, random_state=0, init=centers,algorithm="full", tol=0, max_iter=300).fit(datapoints)
  
    # print(",".join(["%d" % node for node in indices]))
    # for row in kmeans.cluster_centers_:
    #     print(",".join(["%.4f" % edge if "%.4f" % edge != "-0.0000" else "0.0000" for edge in row]))
    # print("---------")
    
    # setup
    np.random.seed(0)
    i = 0
    centeroids = np.zeros([k, datapoints.shape[1]])
    centroid_ids = []
    # randomize initial centeroid
    j = np.random.choice(datapoints.shape[0])

    centeroids[i] = datapoints[j]
    centroid_ids.append(str(j))

    # loop until all centeroids are chosen
    while i < k - 1:
        # build list of min((xl-uj)^2) for 0<=j<=i for every xl
        D_l = np.sum(np.power(
            np.swapaxes(datapoints * np.ones([i + 1, *datapoints.shape]), 0, 1)
            -
            centeroids[:i + 1]
            , 2), axis=2)
        D_l = np.min(D_l, axis=1)
        # create likelihood vector from D_l
        p = D_l / np.sum(D_l)
        # add a new randomized centeroid with p
        j = np.random.choice(datapoints.shape[0], p=p)
        i += 1
        centeroids[i] = datapoints[j]
        centroid_ids.append(str(j))

    print(",".join(centroid_ids))

    return centeroids


def main():
    k, goal, file_name, observations = read_params()
    try:
        datapoints = extract_inputs(file_name)
    except:
        raise RuntimeError("Invalid Input!")

    if goal_enum[goal] == 'spk':
        if k == 0:
            # Running spkmeans in order to set the number of k
            matrix = spk.fit(goal, file_name, observations.shape[0], observations.shape[1])
            # print(matrix)
            k = matrix['k']
            datapoints = np.array(matrix['edges'])

        # Running kmeans++ with the k above
        centroids = kmeans_pp(k, datapoints)
        centroids = kmeans_c.fit(k, len(datapoints), len(datapoints[0]), 300, 0.0, datapoints.tolist(),
                                 centroids.tolist())

        for row in centroids['edges']:
            print(",".join(["%.4f" % edge if "%.4f" % edge != "-0.0000" else "0.0000" for edge in row]))
    else:
        matrix = spk.fit(goal, file_name, observations.shape[0], observations.shape[1])

        if (goal_enum[goal] == 'jacobi'):
            print(",".join(["%.4f" % node for node in matrix['nodes']]))
        for row in matrix['edges']:
            print(",".join(["%.4f" % edge if "%.4f" % edge != "-0.0000" else "0.0000" for edge in row]))


if __name__ == "__main__":
    main()
