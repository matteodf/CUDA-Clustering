
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef KMEANS_UTILS_H_
#define KMEANS_UTILS_H_

#define dist2(a, b) (((a)->x - (b)->x) * ((a)->x - (b)->x) + ((a)->y - (b)->y) * ((a)->y - (b)->y)) 

typedef struct {
    float x;
    float y;
    int group;
} POINT;

void gen_xy(POINT* pts, int num_pts, float radius);
int bisectionSearch(float* x, int n, float v);
void initClusters(POINT* pts, int num_pts, POINT* centroids, int num_clusters);
void printResult(POINT* pts, int num_pts, POINT* centroids, int num_clusters);

#endif