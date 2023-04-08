#include "spkmeans.h"

/*gets: 2 datapoints.
returns: distance between datapoints*/
double d(vector *p, vector *q)
{
    double sum = 0.0;
    cord *curP = p->cords;
    cord *curQ = q->cords;
    while (curP != NULL)
    {
        sum += pow(curP->value - curQ->value, 2);
        curP = curP->next;
        curQ = curQ->next;
    }
    return sqrt(sum);
}
/*gets two vectors arrays and arrays len K
returns true if d(prev_centroids[i], centroids[i])<epsilon for every i. else return false*/
int isConverge(vector *prev_centroids[], vector *centroids[], unsigned K)
{
    unsigned i = 0;
    for (; i < K; i++)
    {
        if (d(prev_centroids[i], centroids[i]) >= 0.0)
        {
            return 0;
        }
    }
    return 1;
}
/*gets: centroids array, NULL array datapointsByCentroids, all datapoints, K
after: each datapoint is is in datapointsByCentroids[i], when centroid[i] is the closest centroid to datapoint
returns 0 on succesful run, otherwise 1 */
int assignDatapointsToCentroid(vector *centroids[], vector *dataPointsByCentroids[], vector *datapoints, unsigned K)
{
    vector *curDatapoint = datapoints;
    vector **curVectors = malloc(K * sizeof(vector *)); /* K pointers for the building process of dataByCentroids*/
    if (curVectors == NULL)
    {
        return 1;
    }

    while (curDatapoint != NULL)
    {
        double minDist = 20000.0;
        unsigned i = 0, minCentroid = K + 1;
        /*find centroid with minimal distance from curDatapoint */
        for (; i < K; i++)
        {
            double dist = d(curDatapoint, centroids[i]);
            if (dist < minDist)
            {
                minCentroid = i;
                minDist = dist;
            }
        }
        /* assign curDatapoint to closest centroid*/
        if (dataPointsByCentroids[minCentroid] == NULL)
        { /* first point assigned to this centroid*/
            dataPointsByCentroids[minCentroid] = (vector *)malloc(sizeof(vector));
            if (dataPointsByCentroids[minCentroid] == NULL)
            {
                return 1;
            }
            dataPointsByCentroids[minCentroid]->cords = curDatapoint->cords;
            dataPointsByCentroids[minCentroid]->next = NULL;
            curVectors[minCentroid] = dataPointsByCentroids[minCentroid];
        }
        else
        {
            curVectors[minCentroid]->next = (vector *)malloc(sizeof(vector));
            if (curVectors[minCentroid]->next == NULL)
            {
                return 1;
            }
            curVectors[minCentroid] = curVectors[minCentroid]->next;
            curVectors[minCentroid]->cords = curDatapoint->cords;
            curVectors[minCentroid]->next = NULL;
        }
        curDatapoint = curDatapoint->next;
    }
    free(curVectors);
    return 0;
}
/*gets:K length array centroids,
        K lenght array datapointsBycentroids with all datapoints allocated to centroid i in index i
        arrays len K
        datapoints dimension
calculate new K centroids according to datapointsByCentroids[] and stores in centroids[]
return 0 on succesful run, otherwise 1 */
int updateCentroids(vector *centroids[], vector *datapointsByCentroids[], unsigned K, int dimension)
{
    unsigned i = 0;
    for (; i < K; i++)
    {
        /*initialize new centroid*/
        int j = 1;
        cord *curCord;
        vector *newCentroid = (vector *)malloc(sizeof(vector));
        if (newCentroid == NULL)
        {
            return 1;
        }
        newCentroid->cords = (cord *)malloc(sizeof(cord));
        newCentroid->next = NULL;
        if (newCentroid->cords == NULL)
        {
            free(newCentroid);
            return 1;
        }

        newCentroid->cords->value = 0.0;
        curCord = newCentroid->cords;
        for (; j < dimension; j++)
        {
            curCord->next = (cord *)malloc(sizeof(cord));
            curCord = curCord->next;
            if (curCord == NULL)
            {
                freeCords(newCentroid->cords);
                free(newCentroid);
                return 1;
            }
            curCord->value = 0.0;
        }
        curCord->next = NULL;
        /* if no datapoints assigned to this centroid*/
        if (datapointsByCentroids[i] == NULL)
        {
            centroids[i] = newCentroid;
        }
        else
        {
            /* loop over all datapoints assigned to centroid i and sum cords*/
            vector *curDatapoint = datapointsByCentroids[i];
            int pointCount = 0;
            while (curDatapoint != NULL)
            {
                cord *tmpCord = curDatapoint->cords;
                pointCount++;
                curCord = newCentroid->cords;
                while (tmpCord != NULL)
                { /* loop over cords of new centroid and current datapoint*/
                    curCord->value += tmpCord->value;
                    curCord = curCord->next;
                    tmpCord = tmpCord->next;
                }
                curDatapoint = curDatapoint->next;
            }
            /* divide each cord by pointCount to get average*/
            curCord = newCentroid->cords;
            while (curCord != NULL)
            {
                curCord->value = curCord->value / pointCount;
                curCord = curCord->next;
            }
            centroids[i] = newCentroid;
        }
    }
    return 0;
}
/*gets: datapoints array, K array length
prints datapoints*/
void printCentroids(vector *centroids[], unsigned K)
{
    unsigned i = 0;
    for (; i < K; i++)
    {
        cord *curCord = centroids[i]->cords;
        while (curCord->next != NULL)
        {
            printf("%.4f,", curCord->value);
            curCord = curCord->next;
        }
        printf("%.4f\n", curCord->value);
    }
}
/*gets: pointer to cord
frees allocated memory of this cord and all next cords.*/
void freeCords(cord *cor)
{
    if (cor != NULL && cor->next != NULL)
    {
        freeCords(cor->next);
    }
    if (cor != NULL)
        free(cor);
}
/*gets:vector pointer and includeCords flag
if includeCords==0 frees all memory allocated to vector chain.
else- frees all memory of vectors AND their cords. */
void freeVectorGroup(vector *vec, int includeCords)
{
    if (vec != NULL && vec->next != NULL)
    {
        freeVectorGroup(vec->next, includeCords);
    }
    if (vec != NULL && includeCords)
    {
        freeCords(vec->cords);
    }
    if (vec != NULL)
        free(vec);
}
/*returns k centroids by kmeans algorithm, starting with the centroids recieved as arguments*/
vector **FindKmeans(unsigned K, unsigned dimension, vector *datapoints, vector **centroids)
{
    /*varaible declaration*/
    unsigned iter = 300;
    int mallocError = 0;
    unsigned i = 0, iteration_count = 0;
    vector **prev_centroids = malloc(K * sizeof(vector *));
    vector **datapointsByCentroids = malloc(K * sizeof(vector *));

    /*check mallocs were valid*/
    if (prev_centroids == NULL || datapointsByCentroids == NULL)
        mallocError = 1;

    if (!mallocError)
    {
        /*initialize pointers array to null*/
        for (; i < K; i++)
            datapointsByCentroids[i] = NULL;
        /*calculate K centroids if no mem errors occured so far: */
        do
        {
            unsigned i = 0;
            /*assign each datapoint to current closest centroid*/
            mallocError = assignDatapointsToCentroid(centroids, datapointsByCentroids, datapoints, K);
            if (mallocError)
            { /*unsuccesful malloc in method. free datapoints allocated before error*/
                for (; i < K; i++)
                    if (datapointsByCentroids[i] != NULL)
                        freeVectorGroup(datapointsByCentroids[i], 0);
                break;
            }
            /*save current centroids in prev_centroid, and free allocated data in prev_centroid*/
            for (; i < K; i++)
            {
                if (iteration_count != 0)
                { /*free allocated data if not first iteration*/
                    freeCords(prev_centroids[i]->cords);
                    free(prev_centroids[i]);
                }
                prev_centroids[i] = centroids[i];
            }
            /*calculate new centroids with thr new datapointsByCentroids*/
            mallocError = updateCentroids(centroids, datapointsByCentroids, K, dimension);
            if (mallocError) /*error occured in allocation memory for new centroids*/
            {
                for (; i < K; i++) /*free dataPointsByCentroids that were already allocated in this iteration*/
                    freeVectorGroup(datapointsByCentroids[i], 0);
                break;
            }
            i = 0;
            /*free allocated data in datapointsByCentroid and initialize to NULL array*/
            for (; i < K; i++)
            {
                freeVectorGroup(datapointsByCentroids[i], 0);
                datapointsByCentroids[i] = NULL;
            }
            /*count iterations*/
            iteration_count++;
        } while (!isConverge(prev_centroids, centroids, K) && iteration_count < iter);
    }
    /*free all memory allocated.it is needed in both cases of mallocError*/
    i = 0;
    for (; i < K; i++)
    { /*free prev_centroids */
        if (prev_centroids[i] != NULL)
        {
            if (prev_centroids[i]->cords != NULL)
                freeCords(prev_centroids[i]->cords);
            free(prev_centroids[i]);
        }
    }
    /*free rest of memory, excepts for returned centroids*/
    if (prev_centroids != NULL)
        free(prev_centroids);
    if (datapointsByCentroids != NULL)
        free(datapointsByCentroids);
    if (datapoints != NULL)
        freeVectorGroup(datapoints, 1);
    /*return program value in dependency on errors occured*/
    if (!mallocError)
        return centroids;
    else
        return NULL;
}
