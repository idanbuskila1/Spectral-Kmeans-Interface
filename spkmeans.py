import sys
import math
import numpy as np
import pandas as pd
import mykmeanssp as kmeans

# initializes K,N,epsilon, datapoints and validates input
def initialize():
    argc = len(sys.argv)
    K = -1
    dimension = 0
    goal_index = 1

    if argc == 4:  # K was passed, read it
        goal_index = 2
        try:
            K = int(sys.argv[1])
        except ValueError:
            print("Invalid number of clusters!")
            raise SystemExit

    elif argc < 3:  # wrong number of args
        print("An Error Has Occured")
        raise SystemExit

    try:  # read goal
        goal = sys.argv[goal_index]
    except ValueError:
        print("Invalid goal!")
        raise SystemExit

    # initialize datapoints/matrix from file
    datapoints = pd.read_csv(sys.argv[goal_index+1])
    header = ["d"+str(i) for i in range(datapoints.shape[1])]
    header[0] = "index"
    datapoints = pd.read_csv(sys.argv[goal_index+1], names=header)
    N = datapoints.shape[0]
    dimension = datapoints.shape[1]-1

    return (K, datapoints, N, dimension,goal)


def distance(p, q):
    x = np.subtract(p, q)
    return np.sqrt(np.sum(np.multiply(x, x)))


def updateDistances(cur_distances, centroids_indexes, datapoints):
    index = centroids_indexes[-1]
    new_centroid = np.array(datapoints.loc[datapoints.index == index])[0, 1:]
    for i in range(cur_distances.shape[0]):
        j = cur_distances[i, 1]
        point = np.array(datapoints.loc[datapoints.index == j])[0, 1:]
        dist = distance(point, new_centroid)
        if cur_distances[i, 0] > dist:
            cur_distances[i, 0] = dist


def addCentroid(min_distances, datapoints):
    dist_sum = np.sum(min_distances[:, 0])
    distribution = np.array([float(min_distances[i, 0] / dist_sum)
                            for i in range(min_distances.shape[0])])
    indexes = np.array(datapoints.index)
    centroid_index = np.random.choice(indexes, p=distribution)
    return centroid_index


def initialize_centroids(K, datapoints):
    centroids_indexes = np.array([])
    min_distances = np.array([[math.inf, datapoints.iloc[i, 0]]
                             for i in range(datapoints.shape[0])])
    # randomly choose first centroid
    rand_index = np.random.choice(np.array(datapoints.iloc[:, 0]))
    centroids_indexes = np.append(centroids_indexes, rand_index)
    updateDistances(min_distances, centroids_indexes, datapoints)

    # choose another K-1 centroids by kmeans++ algorithm
    for i in range(1, K):
        centroids_indexes = np.append(centroids_indexes, addCentroid(min_distances, datapoints))
        updateDistances(min_distances, centroids_indexes, datapoints)
    # return centroids
    centroids = np.array(
        datapoints.loc[datapoints.index == centroids_indexes[0]])[0, 1:]
    for i in range(1, K):
        centroid = np.array(
            datapoints.loc[datapoints.index == centroids_indexes[i]])[0, 1:]
        centroids = np.vstack([centroids, centroid])
    return centroids, centroids_indexes


def printRow(row, N, format):
    for i in range(N-1):
        print(format % row[i], end=",")
    print(format % row[N-1])


### MAIN ###
np.random.seed(0)
K, data, N, dimension, goal = initialize()
datapoints = data.values.tolist()

if goal=="spk":
    gl = kmeans.gl(datapoints)
    eigns = kmeans.jacobi(gl)
    U = kmeans.getKLowestEignvectors(K,eigns)
    K = len(U[0])#if K was passed its not changing.
    col=np.zeros((N,1))
    for i in range(N):col[i,0]=i
    U=np.hstack((col,U))
    centroids, indexes = initialize_centroids(K, pd.DataFrame(U))
    U = U[:,1:].tolist()
    centroids = np.array(kmeans.Kmeans(K, N, K, U, centroids.tolist()))
    printRow(indexes,K,"%d")
    for i in range(centroids.shape[0]): printRow(centroids[i], centroids.shape[1], "%.4f")
else:
    if goal == "wam":
        ret = np.array(kmeans.wam(datapoints))
    elif goal == "ddg":
        ret = np.array(kmeans.ddg(datapoints))
    elif goal == "gl":
        ret = np.array(kmeans.gl(datapoints))
    elif goal == "jacobi":
        ret = np.array(kmeans.jacobi(datapoints))
    for i in range(ret.shape[0]): printRow(ret[i], ret.shape[1], "%.4f")
