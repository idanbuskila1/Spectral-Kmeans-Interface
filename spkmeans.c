#include "spkmeans.h"

/*GETS: WAM, DDG of the points
RETURNS: Graph Laplacian of the points, inplace on argument wam*/
double **createGraphLaplacian(double **wam, double **ddg, int N)
{
    int i = 0, j;
    for (; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i == j)
                wam[i][i] = ddg[i][i];
            else
                wam[i][j] *= (-1);
        }
    }
    return wam;
}
/*GETS:n dimensional vector
RETURNS: vector sum*/
double sum(double *row, int n)
{
    double sum = 0.0;
    int i = 0;
    for (; i < n; i++)
        sum += row[i];
    return sum;
}
/* GETS:WAM of points, N matrix size
 RETURNS: Diagonal Degree Matrix of points*/
double **createDegreeMatrix(double **wam, int N)
{
    double **ddg = (double **)malloc(sizeof(double *) * N);
    int i = 0;
    if (!ddg)
    {
        printf("An Error has occured");
        exit(1);
    }
    for (; i < N; i++)
    {
        ddg[i] = (double *)calloc(N, sizeof(double));
        if (!ddg[i])
        {
            printf("An Error has occured");
            exit(1);
        }
        ddg[i][i] = sum(wam[i], N);
    }
    return ddg;
}
/* GETS: 2 vectors of the dimension x,y
 RETURNS: sum((x-y)^2)*/
double squaredDistance(double *x, double *y, int dimension)
{
    int i = 0;
    double dist = 0;
    for (; i < dimension; i++)
    {
        dist += pow(x[i] - y[i], 2);
    }
    return dist;
}
/* GETS: pointer to N data points of the dimension
 RETURNS: Weighted Adjacancy Matrix of the points*/
double **createWAM(int N, int dimension, double **data)
{
    int i = 0, j;
    double exponent;
    double **wam = (double **)malloc(sizeof(double *) * N);
    if (!wam)
    {
        printf("An Error has occured");
        exit(1);
    }
    /*initialize matrix*/
    for (; i < N; i++)
    {
        wam[i] = (double *)malloc(sizeof(double) * N);
        if (!wam[i])
        {
            printf("An Error has occured");
            exit(1);
        }
        wam[i][i] = 0.0;
    }
    /*calculate weights*/
    for (i = 0; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            exponent = squaredDistance(data[i], data[j], dimension) * (-0.5);
            exponent = exp(exponent);
            wam[i][j] = exponent;
            wam[j][i] = wam[i][j];
        }
    }

    return wam;
}
/* GETS: string of vector values saparated by comma
 RETURNS: array represenation of the vector*/
double *split(char *line, int dimension)
{
    int i = 1;
    char *token = strtok(line, ",");
    double *ret = (double *)malloc(sizeof(double) * dimension);
    if (!ret)
    {
        printf("An Error has occured");
        exit(1);
    }

    ret[0] = atof(token);
    for (; i < dimension; i++)
        ret[i] = atof(strtok(NULL, ","));
    return ret;
}
/* GETS: file name with N points of the dimension
 RETURNS: N x dimension matrix of the data*/
double **ReadFile(char *fileName, int *N, int *dimension)
{
    FILE *fp;
    char *line = NULL, *token;
    size_t len = 0;
    ssize_t read;
    int i = 0;
    double **ret;

    fp = fopen(fileName, "r");
    if (fp == NULL)
    {
        printf("An Error has occured");
        exit(1);
    }
    /* get number of datapoints/matrix rows */
    while ((read = getline(&line, &len, fp)) != -1)
        (*N)++;
    rewind(fp);
    /*get dimension of datapoints/ matrix columns*/
    read = getline(&line, &len, fp);
    token = strtok(line, ",");
    while (token != NULL)
    {
        (*dimension)++;
        token = strtok(NULL, ",");
    }
    rewind(fp);
    /*get data into ret 2 dimensional array*/
    ret = (double **)malloc(sizeof(double *) * (*N));
    if (!ret)
    {
        printf("An Error has occured");
        exit(1);
    }
    for (; i < *N; i++)
    {
        read = getline(&line, &len, fp);
        ret[i] = split(line, *dimension);
    }

    fclose(fp);
    if (line)
        free(line);
    return ret;
}
/* GETS: matrix with n rows and m columns
RETURNS: void. prints the matrix line by line, values saparated with comma*/
void printMatrix(int n, int m, double **matrix)
{
    int i = 0, j;
    for (; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (j != m - 1)
                printf("%.4f,", matrix[i][j]);
            else
                printf("%.4f\n", matrix[i][j]);
        }
    }
}
/* GETS: pointers to 2 pointers of pairs of int index and double value
 RETURNS: comparison of value*/
int PairComparator(const void *p, const void *q)
{
    double val1 = (*((pair **)p))->value;
    double val2 = (*((pair **)q))->value;
    if (val1 > val2)
        return 1;
    if (val1 < val2)
        return -1;
    return 0;
}
/* GETS: pairs of N eignvalues sorted
 RETURNS: the i index of biggest jump in values of eignvalues i,i+1*/
int EignvaluesGapHuristic(pair **eignvalues, int N)
{
    int K = -1, i = 0;
    double maxGap = -1.0, temp;

    for (; i < ((int)(N / 2)); i++)
    {
        temp = fabs((eignvalues[i]->value) - (eignvalues[i + 1]->value));
        if (temp > maxGap)
        {
            maxGap = temp;
            K = i;
        }
    }
    /*counting from 1 and not zero*/
    return K + 1;
}
/* GETS: N eignvectors as columns, pairs of sorted eignvalue and original index, desired K
 RETURNS:matrix with K columns of the eignvectors corresponding with the K min eignvalues*/
double **createKEignvectorsMatrix(int N, int K, double **eignvectors, pair **SortPremutaion)
{
    int i, j;
    double **ret = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
    {
        ret[i] = (double *)malloc(K * sizeof(double));
        for (j = 0; j < K; j++)
        {
            ret[i][j] = eignvectors[i][SortPremutaion[j]->index];
        }
    }
    return ret;
}

/*C INTERFACE*/
int main(int argc, char *argv[])
{
    int N = 0, dimension = 0;
    char *goal = argv[1];
    double **data = ReadFile(argv[argc - 1], &N, &dimension);
    double **wam, **ddg, **gl, **result;

    if (strcmp("jacobi", goal) == 0)
    {
        result = jacobiProcess(N, data);
        printMatrix(N + 1, N, result);
        freeMatrix(result, N + 1);
        return 0;
    }
    wam = createWAM(N, dimension, data);
    if (strcmp("wam", goal) == 0)
    {
        printMatrix(N, N, wam);
        freeMatrix(wam, N);
        freeMatrix(data, N);
        return 0;
    }
    ddg = createDegreeMatrix(wam, N);
    if (strcmp("ddg", goal) == 0)
    {
        printMatrix(N, N, ddg);
        freeMatrix(ddg, N);
        freeMatrix(wam, N);
        freeMatrix(data, N);
        return 0;
    }
    gl = createGraphLaplacian(wam, ddg, N);
    if (strcmp("gl", goal) == 0)
    {
        printMatrix(N, N, gl);
        freeMatrix(gl, N);
        freeMatrix(ddg, N);
        freeMatrix(data, N);
        return 0;
    }
    return 1;
}
