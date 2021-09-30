
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifndef KMEANS_UTILS_H_
#define KMEANS_UTILS_H_

typedef struct {
    float* x;
    float* y;
    int* group;
} POINT;

#define dist2(ax, ay, bx, by) (((ax) - (bx)) * ((ax) - (bx)) + ((ay) - (by)) * ((ay) - (by)))

void mallocPS(POINT* ps, int n);
void freePS(POINT* ps);
void gen_xy(POINT pts, int num_pts, float radius);
int nearest(POINT pts, int ptPos, POINT centroids, int n_cluster);
int bisectionSearch(float* x, int n, float v);
void initClusters(POINT pts, int num_pts, POINT centroids, int num_clusters);
void printResult(POINT pts, int num_pts, POINT centroids, int num_clusters);
void printDebug(POINT centroids, int num_clusters);

#endif