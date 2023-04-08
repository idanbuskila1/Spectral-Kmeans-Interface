#ifndef KMEANS_H_
#define KMEANS_H_
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
/*structs def*/
typedef struct cord
{
    double value;
    struct cord *next;
} cord;
typedef struct vector
{
    struct vector *next;
    struct cord *cords;
} vector;
typedef struct pair
{
    int index;
    double value;
} pair;
/*prototypes*/
/*spkmeans.c*/
double **createGraphLaplacian(double **wam, double **ddg, int N);
double sum(double *row, int n);
double **createDegreeMatrix(double **wam, int N);
double squaredDistance(double *x, double *y, int dimension);
double **createWAM(int N, int dimension, double **data);
double *split(char *line, int dimension);
double **ReadFile(char *fileName, int *N, int *dimension);
void printMatrix(int n, int m, double **matrix);
int PairComparator(const void *p, const void *q);
int EignvaluesGapHuristic(pair **eignvalues, int N);
double **createKEignvectorsMatrix(int N, int K, double **eignvectors, pair **SortPremutaion);
/*jacobi.c*/
double **jacobiProcess(int N, double **A);
double **calculateRotationMatrix(int N, double **A, double *c, double *s, int *i, int *j);
void FindPivot(int *i, int *j, int N, double **A);
double **jacobiStep(double **A, int N, double c, double s, int i, int j);
int isConverged(double **A, double **A2, int N);
void freeMatrix(double **A, int N);
double **MultuplyMatrices(double **A, double **B, int N);
/*kmeans.c*/
void freeCords(cord *cor);
void freeVectorGroup(vector *vec, int includeCords);
vector **FindKmeans(unsigned K, unsigned dimension, vector *datapoints, vector **centroids);
#endif
