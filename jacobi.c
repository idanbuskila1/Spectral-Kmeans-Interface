#include "spkmeans.h"
/* GETS: 2 square Matrices of NxN
 RETURNS: C=AB*/
double **MultuplyMatrices(double **A, double **B, int N)
{
    double **C = (double **)malloc(N * sizeof(double *));
    int i = 0, j, k;
    for (; i < N; i++)
    {
        C[i] = (double *)malloc(N * sizeof(double));
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
    return C;
}
/* GETS: Matrix A with n rows
 RETURNS: void. frees it*/
void freeMatrix(double **A, int N)
{
    int i = 0;
    for (; i < N; i++)
        free(A[i]);
    free(A);
}
/* GETS: Matrix A, Matrix A2 after jacobi step, rows number N
 RETURNS: 1 if they are converged after the step, 0 otherwise*/
int isConverged(double **A, double **A2, int N)
{
    double offA = 0.0, offA2 = 0.0;
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i != j)
            { /*sum of diagonal values*/
                offA += pow(A[i][j], 2);
                offA2 += pow(A2[i][j], 2);
            }
        }
    }
    if (fabs(offA - offA2) <= 0.00001) /*convergence rule*/
        return 1;
    return 0;
}
/* GETS:square Matrix A NxN, pivot index (i,j), step parameters s,c
 RETURNS:A2 = jacobi algorithm step on A*/
double **jacobiStep(double **A, int N, double c, double s, int i, int j)
{
    int k, l;
    double val;
    double **A2 = (double **)malloc(N * sizeof(double *));
    if (!A2)
    {
        printf("An Error Has Occured!");
        exit(1);
    }
    /*duplicate A*/
    for (k = 0; k < N; k++)
    {
        A2[k] = (double *)malloc(N * sizeof(double));
        for (l = 0; l < N; l++)
        {
            A2[k][l] = A[k][l];
        }
    }
    /*change values in i'th column and row after rotation*/
    for (k = 0; k < N; k++)
    {
        val = (c * A[k][i]) - (s * A[k][j]);
        A2[k][i] = val;
        A2[i][k] = val;
    }
    /*change values in j'th column and row after rotation*/
    for (k = 0; k < N; k++)
    {
        val = (c * A[k][j]) + (s * A[k][i]);
        A2[k][j] = val;
        A2[j][k] = val;
    }
    /*change rest of values after rotation*/
    A2[i][i] = (c * c * A[i][i]) + (s * s * A[j][j]) - (2 * s * c * A[i][j]);
    A2[j][j] = (s * s * A[i][i]) + (c * c * A[j][j]) + (2 * s * c * A[i][j]);
    A2[i][j] = 0;
    A2[j][i] = 0;

    return A2;
}
/* GETS:square Matrix A NxN, empty pointers i,j
 RETURNS:updates i,j to pivot values of A*/
void FindPivot(int *i, int *j, int N, double **A)
{
    double max = 0.0;
    int k = 0, l = 0;

    for (; k < N; k++)
    {
        for (l = 0; l < N; l++)
        {
            if (k != l && max < fabs(A[k][l]))
            { /*pivot condition*/
                max = fabs(A[k][l]);
                *i = k;
                *j = l;
            }
        }
    }
}
/* GETS: square Matrix A NxN, empty pointers i,j,s,c
 RETURNS: finds parameters s,c, pivot (i,j), and returns according P rotation matrix*/
double **calculateRotationMatrix(int N, double **A, double *c, double *s, int *i, int *j)
{
    int k = 0;
    double t, theta;
    double **P = (double **)malloc(N * sizeof(double *));
    if (!P)
    {
        printf("An Error has occured");
        exit(1);
    }

    /*initialize P as Identity matrix*/
    for (; k < N; k++)
    {
        P[k] = (double *)calloc(N, sizeof(double));
        P[k][k] = 1;
    }

    /*calculate s,c,t*/
    FindPivot(i, j, N, A);
    theta = (A[*j][*j] - A[*i][*i]) / (2 * A[*i][*j]);
    t = 1 / (fabs(theta) + sqrt((theta * theta) + 1));
    if (theta < 0)
        t = t * (-1);
    *c = 1 / sqrt((t * t) + 1);
    *s = t * (*c);

    /*assign values to P*/
    P[*i][*i] = *c;
    P[*j][*j] = *c;
    P[*i][*j] = *s;
    P[*j][*i] = (-1) * (*s);

    return P;
}
/* GETS: square Matrix A NxN
 RETURNS:FREES A, returns eignvalues and corresponding eignvectors, calculated with jacobi algorithm*/
double **jacobiProcess(int N, double **A)
{
    double c, s;
    double **P2, **temp;
    int i, j, k = 0, iterations = 1;
    double **P = calculateRotationMatrix(N, A, &c, &s, &i, &j);
    double **A2 = jacobiStep(A, N, c, s, i, j);
    int flag = isConverged(A, A2, N);
    double **ret = (double **)malloc((N + 1) * sizeof(double *));
    if (!ret)
    {
        printf("An Error has occured");
        exit(1);
    }
    freeMatrix(A, N);
    A = A2;

    while (!flag && iterations < 100)
    {
        P2 = calculateRotationMatrix(N, A, &c, &s, &i, &j);
        A2 = jacobiStep(A, N, c, s, i, j);
        temp = P;
        P = MultuplyMatrices(P, P2, N);
        freeMatrix(temp, N);
        freeMatrix(P2, N);
        flag = isConverged(A, A2, N);
        freeMatrix(A, N);
        A = A2;
        iterations++;
    }

    ret[0] = (double *)malloc(N * sizeof(double));
    for (; k < N; k++)
    {
        ret[0][k] = A[k][k];
        ret[k + 1] = P[k];
    }
    freeMatrix(A, N);
    free(P);
    return ret;
}
