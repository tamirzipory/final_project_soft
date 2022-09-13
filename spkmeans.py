import sys
import numpy as np
import pandas as pd
import mykmeanssp

# --------------- Auxiliary functions --------------- #

def euclid(x, y):  # actually euc^2
    sumUp = 0
    dim = len(x)
    coordinateSub = [0 for init in range(dim)]
    for coordinate in range(dim):
        coordinateSub[coordinate] = x[coordinate] - y[coordinate]
        coordinateSub[coordinate] = coordinateSub[coordinate] * coordinateSub[coordinate]
        sumUp += coordinateSub[coordinate]
    res = (sumUp)

    return res


def getRandCentroids(points, k, indexToKey):

    # what we need:

    n = len(points)
    dim = len(points[1])

    BoolPoints = []
    for i in range(0, n):
        BoolPoints.append(False)  # initialize to false values

    Yindexes = []
    for i in range(0, k):
        Yindexes.append(-1)  # initialize to -1 values

    Dpoints = []
    for i in range(0, n):
        Dpoints.append(-1)  # initialize to -1 values

    Ppoints = []
    for i in range(0, n):
        Ppoints.append(0)  # initialize to 0 values
    toprint = ""

    # actual function:

    SumD = 0
    np.random.seed(0)

    for z in range(0, k):
        SelectedIndex = 0
        if (z == 0):
            SelectedIndex = np.random.choice(n, 1)[0]  # select a random number Rand from 0 to n-1
            toprint = toprint + str(indexToKey[SelectedIndex])
            Yindexes[z] = SelectedIndex
            BoolPoints[SelectedIndex] = True

        else:  # z!=0
            SelectedIndex = np.random.choice(n, 1, p=Ppoints)[0]
            toprint = toprint + "," + str(indexToKey[SelectedIndex])
            Yindexes[z] = SelectedIndex
            BoolPoints[SelectedIndex] = True

            if (z == k - 1):  # no need to update Ppoints again
                break
        for i in range(0, n):  # reset Dpoints, Ppoints
            Dpoints[i] = -1
            Ppoints[i] = 0  # seems like shouldnt reset P
        SumD = 0

        # Yindexes[z] = SelectedIndex
        # BoolPoints[SelectedIndex] = True


        for i in range(0, n):
            for j in range(0, z + 1):
                Yj = points[Yindexes[j]]
                CurrDi = euclid(points[i], Yj)
                # CurrDi = CurrDi * CurrDi
                if (Dpoints[i] == -1 or CurrDi < Dpoints[i]):
                    Dpoints[i] = CurrDi
        SumD = sum(Dpoints)

        for index in range(0, n):
            Ppoints[index] = Dpoints[index] / SumD

    CentArray = []
    for i in range(0, k):
        CentArray.append(points[Yindexes[i]])

    # print("toprint:", toprint)
    return CentArray, Yindexes


# print the recieved matrix
def printResult(ArrayToPrint, rowDim, colDim):
    ArrayToPrint = np.round(ArrayToPrint, decimals=4)
    for i in range(rowDim):
        lineString = ""
        for j in range(colDim):
            if j > 0:
                lineString += ","
            if abs(ArrayToPrint[i * colDim + j]) < 0.00005:
                lineString += "%.4f" % abs(ArrayToPrint[i * colDim + j])
            else:
                #lineString += str(ArrayToPrint[i * colDim + j])
                lineString += "%.4f" % (ArrayToPrint[i * colDim + j])
        print(lineString.replace("[", "").replace("]", ""))


# --------------- K-means-PP algorithm --------------- #

def kmeansPP(points, dim, n, k):

    # ----- verify length's correctness ----- #

    # K must be bigger than zero and integer
    if k < 1 or type(k) != int:
        sys.exit("An Error Has Occured")
    if n == 0:  # At least one point is needed
        sys.exit("An Error Has Occured")
    if n <= k:  # K must be smaller than N and not equal to N
        sys.exit("An Error Has Occured")
    if dim == 0:  # Dim of points is zero
        sys.exit("An Error Has Occured")

    # ----- indexing centroids  ----- #
    
    points = pd.DataFrame(points)
    
    indexToKey = []  # indexToKey[i] is the i-th smallest key.
    for i in range(len(points.index)):
        indexToKey.append(int(points.index[i]))

    # ----- choosing random initial centroids  ----- #

    regPoints = np.array(points).tolist()
    _, centroidsIndexes = getRandCentroids(regPoints, k, indexToKey)

    # ----- print random initial centroids  ----- #

    print(",".join([str(index) for index in centroidsIndexes]))

    # ----- create arrays for kmeans algo in c  ----- #

    centroids_vals = []
    for cent_index in centroidsIndexes:
        centroids_vals.append(points.loc[cent_index])
    centroids = pd.DataFrame(centroids_vals)

    centroidsResult = np.zeros(dim * k)

    # ----- calls kmeans in c through fit function in API  ----- #
    # maxIter = 300

    mykmeanssp.fit(points, centroidsResult, centroids, dim, k, 300)

    # ----- printing results ----- #
    printResult(centroidsResult, k, dim)


# --------------- main  --------------- #

if len(sys.argv) != 4:
    sys.exit("An Error Has Occured")

# ----- calling SPKmeans in c ----- #
kStr = sys.argv[1]
goalStr = sys.argv[2]
fileNameStr = sys.argv[3]

n = 0
with open(fileNameStr, 'r') as fin:
  for line in fin:
    n += 1

# TODO: get n,k from file!!!
#n = 10
k = int(kStr)
x,y = n+1, n
#results = np.zeros((n+1, n))
results = np.zeros((x, y))
k_ans = np.zeros((1, 1))

mykmeanssp.fitToSPKmeans(kStr, goalStr, fileNameStr, results, k_ans)

# ----- taking k as int ----- #
k = int(kStr)

if goalStr == "spk":
  k_ans = int(k_ans[0][0])
  results.shape = (x*y, )
  results = results[0:n*k_ans]
  results.shape = (n, k_ans)
  kmeansPP(results, k_ans, n, k_ans)

elif goalStr == "jacobi":   # (n+1) * n
  results.shape = (x*y, 1)
  results = results[0:(n+1)*n]
  printResult(results, n+1, n)

elif goalStr == "wam":
  results.shape = (x*y, )
  results = results[0:n*n]
  printResult(results, n, n)

elif goalStr == "ddg":
  results.shape = (x*y, )
  results = results[0:n*n]
  printResult(results, n, n)

elif goalStr == "lnorm":
  results.shape = (x*y, )
  results = results[0:n*n]
  printResult(results, n, n)
