#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifndef KMEANS_UTILS_H_
#define KMEANS_UTILS_H_

typedef struct {
    float x;
    float y;
    int group;
} POINT;

inline double seconds() {
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

#define dist2(a, b) (((a)->x - (b)->x) * ((a)->x - (b)->x) + ((a)->y - (b)->y) * ((a)->y - (b)->y))

void printDebug(POINT* centroids, int num_clusters);
void genRandomPoints(POINT* pts, int num_pts, float radius);
int findNumRows(char* filename);
void readPointsFromFile(POINT* pts, char* filename);
int nearest(POINT* pt, POINT* cent, int n_cluster);
int bisectionSearch(float* x, int n, float v);
void initClustersRandom(POINT* pts, int num_pts, int num_clusters);
void initClusters(POINT* pts, int num_pts, POINT* centroids, int num_clusters);
void writePointsToFile(char* filename, POINT* pts, int num_pts, POINT* centroids, int num_clusters);
void printResult(POINT* pts, int num_pts, POINT* centroids, int num_clusters);

#endif