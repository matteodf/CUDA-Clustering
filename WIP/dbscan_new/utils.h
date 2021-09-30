#ifndef DBSCAN_UTILS_H_
#define DBSCAN_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <cmath>

#include "Points.h"

#define dist(ax, ay, bx, by) sqrt(((ax) - (bx)) * ((ax) - (bx)) + ((ay) - (by)) * ((ay) - (by)))

const int NOISE = -2;
const int NOT_CLASSIFIED = -1;

inline double seconds() {
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

void printResult(Points pts, int num_pts, int num_clusters);

#endif