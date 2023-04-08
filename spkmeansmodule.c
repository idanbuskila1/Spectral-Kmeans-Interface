#include "spkmeans.h"
#define PY_SSIZE_T_CLEAN
#include <Python.h>

/*Gets python list of lists, updates N and Dimension,and returns c array of arrays of the same data*/
double **PythonListsToMatrix(PyObject *list, int *N, int *dimension)
{
    *N = PyObject_Length(list);
    *dimension = PyObject_Length(PyList_GetItem(list, 0));
    double **ret = (double **)malloc((*N) * sizeof(double *));
    int i, j;
    PyObject *InnerList, *Item;

    for (i = 0; i < (*N); i++)
    {
        InnerList = PyList_GetItem(list, i);
        ret[i] = (double *)calloc(*dimension, sizeof(double));
        for (j = 0; j < (*dimension); j++)
        {
            Item = PyList_GetItem(InnerList, j);
            ret[i][j] = PyFloat_AsDouble(Item);
        }
    }
    return ret;
}
/*gets n arrays of size m, returns n python lists of size m with the same data*/
PyObject *MatrixToPythonLists(double **matrix, int n, int m)
{
    PyObject *retList = PyList_New(n);
    PyObject *Item, *InnerList;
    int i, j;

    for (i = 0; i < n; i++)
    {
        InnerList = PyList_New(m);
        for (j = 0; j < m; j++)
        {
            Item = Py_BuildValue("f", matrix[i][j]);
            PyList_SetItem(InnerList, j, Item);
        }
        PyList_SetItem(retList, i, InnerList);
    }
    return retList;
}

/*gets python object pointer to list of N lists of size dimension.
returns N connected vectors with |dimension| cords each, which represents the list of lists
if an error occured return NULL */
vector *ListOfListsToVectors(PyObject *lol, unsigned N, unsigned dimension)
{
    PyObject *PyPoint;
    PyObject *item;
    vector *curVec = NULL, *headVec = NULL, *tmpVec;
    cord *headCord = NULL, *curCord;
    unsigned i, j;
    double val;

    for (i = 0; i < N; i++)
    {
        tmpVec = (vector *)malloc(sizeof(vector *));
        curCord = (cord *)malloc(sizeof(cord *));
        if (tmpVec == NULL || curCord == NULL)
        {                                /*if memory error occured*/
            freeVectorGroup(headVec, 1); /*free vectors already malloced(if exist)*/
            if (tmpVec != NULL)
                free(tmpVec);
            if (curCord != NULL)
                free(curCord);
            return NULL;
        }
        PyPoint = PyList_GetItem(lol, i);
        for (j = 0; j < dimension; j++)
        {
            item = PyList_GetItem(PyPoint, j);
            val = PyFloat_AsDouble(item);
            curCord->value = val;
            if (j == 0)
                headCord = curCord; /* first cord in vector */
            if (j != dimension - 1)
            {
                curCord->next = (cord *)malloc(sizeof(cord *));
                curCord = curCord->next;
                if (curCord == NULL)
                {                        /*if memory error occured */
                    freeCords(headCord); /*free cords already malloced(if exist)*/
                    return NULL;
                }
            }
            else
                curCord->next = NULL; /* last cord in vector */
        }
        tmpVec->cords = headCord;
        tmpVec->next = NULL;
        if (i == 0) /* first vector in list */
        {
            headVec = tmpVec;
            curVec = headVec;
        }
        else
        {
            curVec->next = tmpVec;
            curVec = curVec->next;
        }
    }
    return headVec;
}

/*gets array with K pointers to vectors with |dimension| cords.
returns python object pointer to a python list of K lists of size |dimension|"*/
PyObject *VectorArrayToListOfLists(vector **array, unsigned K, unsigned dimension)
{
    PyObject *retList = PyList_New(K);
    unsigned i = 0, j = 0;
    PyObject *python_float;
    cord *tmp;

    for (; i < K; i++)
    {
        PyObject *innerList = PyList_New(dimension);
        tmp = array[i]->cords;
        for (; j < dimension; j++)
        {
            python_float = Py_BuildValue("f", tmp->value);
            PyList_SetItem(innerList, j, python_float);
            tmp = tmp->next;
        }
        PyList_SetItem(retList, i, innerList);
        j = 0;
    }
    return retList;
}
// GETS: data points from python
// RETURNS: weighted adjacancy matrix of them
static PyObject *wam(PyObject *self, PyObject *args)
{
    PyObject *PyDataset, *ret;
    double **CDataset, **wam;
    int N, dimension;

    /*parse variables from python*/
    if (!PyArg_ParseTuple(args, "O", &PyDataset))
        return NULL;

    CDataset = PythonListsToMatrix(PyDataset, &N, &dimension);
    wam = createWAM(N, dimension, CDataset);

    ret = MatrixToPythonLists(wam, N, N);

    freeMatrix(CDataset, N);
    freeMatrix(wam, N);
    return ret;
}
// GETS:data points from python
// RETURNS: diagonal degree matrix of them
static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *PyDataset, *ret;
    double **CDataset, **wam, **ddg;
    int N, dimension;

    /*parse variables from python*/
    if (!PyArg_ParseTuple(args, "O", &PyDataset))
        return NULL;

    CDataset = PythonListsToMatrix(PyDataset, &N, &dimension);
    wam = createWAM(N, dimension, CDataset);
    ddg = createDegreeMatrix(wam, N);
    ret = MatrixToPythonLists(ddg, N, N);

    freeMatrix(CDataset, N);
    freeMatrix(wam, N);
    freeMatrix(ddg, N);
    return ret;
}
// GETS:data points from python
// RETURNS: graph laplacian of them
static PyObject *gl(PyObject *self, PyObject *args)
{
    PyObject *PyDataset, *ret;
    double **CDataset, **wam, **ddg, **gl;
    int N, dimension;

    /*parse variables from python*/
    if (!PyArg_ParseTuple(args, "O", &PyDataset))
        return NULL;

    CDataset = PythonListsToMatrix(PyDataset, &N, &dimension);
    wam = createWAM(N, dimension, CDataset);
    ddg = createDegreeMatrix(wam, N);
    gl = createGraphLaplacian(wam, ddg, N);
    ret = MatrixToPythonLists(gl, N, N);

    freeMatrix(CDataset, N);
    freeMatrix(gl, N);
    freeMatrix(ddg, N);
    return ret;
}
// GETS:symmetric matrix from python of size NxN
// RETURNS: first row eignvalues, N lower rows are eignvetors as columns matrix
static PyObject *jacobi(PyObject *self, PyObject *args)
{
    PyObject *PyDataset, *ret;
    double **CDataset, **eignvectors;
    int N, dimension;

    /*parse variables from python*/
    if (!PyArg_ParseTuple(args, "O", &PyDataset))
        return NULL;

    CDataset = PythonListsToMatrix(PyDataset, &N, &dimension);
    eignvectors = jacobiProcess(N, CDataset); // frees CDataset
    ret = MatrixToPythonLists(eignvectors, N + 1, N);

    freeMatrix(eignvectors, N);
    return ret;
}
// GETS:desired K, N+1xN matrix with first row eignvalues and lower rows are eignvetors as columns matrix
// RETURNS: NxK matrix of k eignvectors corresponding with lowest k eignvalues as columns
static PyObject *getKLowestEignvectors(PyObject *self, PyObject *args)
{
    int K, i, N, dimension;
    PyObject *data, *retPy;
    double **eignvectors, **retC;
    pair **SortPremutation;

    if (!PyArg_ParseTuple(args, "IO", &K, &data))
        return NULL;

    eignvectors = PythonListsToMatrix(data, &N, &dimension);
    N -= 1; // not including eignvalues row
    SortPremutation = (pair **)malloc(N * sizeof(pair *));
    /*get sort premutation of the eignvalues*/
    for (i = 0; i < N; i++)
    {
        SortPremutation[i] = (pair *)malloc(sizeof(pair));
        SortPremutation[i]->index = i;
        SortPremutation[i]->value = eignvectors[0][i]; /*row 0 is the eignvalues*/
    }
    qsort((void **)SortPremutation, N, sizeof(pair *), PairComparator);
    if (K == -1)
        K = EignvaluesGapHuristic(SortPremutation, N); /*if K not given, determine it by the huristic*/

    /*store K lowest eignvectors as columns in retC and convert to python */
    retC = createKEignvectorsMatrix(N, K, (eignvectors + 1), SortPremutation);
    retPy = MatrixToPythonLists(retC, N, K);
    /*free Data*/
    freeMatrix(eignvectors, N + 1);
    freeMatrix(retC, N);
    for (i = 0; i < N; i++)
        free(SortPremutation[i]);
    free(SortPremutation);
    return retPy;
}
// GETS:N datapoints of the dimension, initial centroids
// RETURNS: final centroids clusteering the data, after kmeans++ algorithm
static PyObject *Kmeans(PyObject *self, PyObject *args)
{
    unsigned K, N, dimension, i = 0;
    PyObject *dataset;
    PyObject *PyCentroids, *retList;
    vector *datapoints, **centroids, *tmpVec, *tmpCentroids;

    /*parse variables from python*/
    if (!PyArg_ParseTuple(args, "IIIOO", &K, &N, &dimension, &dataset, &PyCentroids))
        return NULL;
    centroids = (vector **)malloc(sizeof(vector *) * K); /* array of K pointers*/
    tmpCentroids = ListOfListsToVectors(PyCentroids, K, dimension);
    datapoints = ListOfListsToVectors(dataset, N, dimension);
    tmpVec = tmpCentroids;
    if (datapoints != NULL && tmpCentroids != NULL && centroids != NULL)
    { /* all mallocs were succesful: */
        for (; i < K; i++)
        { /*store centroids adresses in the designated pointers array */
            centroids[i] = tmpVec;
            tmpVec = tmpVec->next;
        }
        /*find k means with initial centroids received from kmeans_pp.py*/
        centroids = FindKmeans(K, dimension, datapoints, centroids);
        if (centroids != NULL)
        {
            /*build PyObject list of lists of the calculated centroids*/
            retList = VectorArrayToListOfLists(centroids, K, dimension);
            /* free all memory */
            i = 0;
            for (; i < K; i++)
            {
                if (centroids[i] != NULL)
                {
                    if (centroids[i]->cords != NULL)
                        freeCords(centroids[i]->cords);
                    free(centroids[i]);
                }
            }
            free(centroids);
            return retList;
        }
    }
    if (centroids != NULL)
        free(centroids);
    PyErr_SetString(PyExc_RuntimeError, "An error has occured");
    return NULL;
}
static PyMethodDef SPKmeansMethods[] = {
    {"wam",
     (PyCFunction)wam,
     METH_VARARGS,
     PyDoc_STR("Gets N datapoints, Returns their Weighted Adjacancy Matrix")},
    {"ddg",
     (PyCFunction)ddg,
     METH_VARARGS,
     PyDoc_STR("Gets N datapoints, Returns their Diagonal Degree Matrix")},
    {"gl",
     (PyCFunction)gl,
     METH_VARARGS,
     PyDoc_STR("Gets N datapoints, Returns their Graph Laplacian Matrix")},
    {"jacobi",
     (PyCFunction)jacobi,
     METH_VARARGS,
     PyDoc_STR("Gets NxN symmetric matrix, Returns its N eignvalues and coresponding eignvectors")},
    {"getKLowestEignvectors",
     (PyCFunction)getKLowestEignvectors,
     METH_VARARGS,
     PyDoc_STR("Gets N eignvectors & eignvalues and K, if K null calculates optimal K and Returns matrix with K lowest eignvectors as columns")},
    {"Kmeans",
     (PyCFunction)Kmeans,
     METH_VARARGS,
     PyDoc_STR("Gets N,K,dimension, K initial centroids and N datapoints of the dimension. Returns K centroids after Kmeans")},
    {NULL, NULL, 0, NULL}

};
static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",                                                                                                                             /* name of module */
    "implements spectral kmeans algorithm, along with jacobi algorithm, and finding weighted adjacancy, degree and graph laplacian matrices", /* module documentation, may be NULL */
    -1,                                                                                                                                       /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    SPKmeansMethods                                                                                                                           /* the PyMethodDef array from before containing the methods of the extension */
};
PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}
