#ifndef DBSCAN_H
#define DBSCAN_H

#include <math.h>

#include "InputReader.h"
#include "Points.h"
#include "utils.h"

#define log2(x) (log(x) / log(2));

#define PRINT_RESULTS false

class DBSCAN {
   public:
    int n, minPts;
    double eps;
    Points points;
    int size;
    int* adjPoints;
    int adjOffset;
    std::vector<std::vector<int>> cluster;
    int clusterIdx;

    DBSCAN(double eps, int minPts, Points points);
    void run();
    void dfs(int now, int c);
    void checkNearPoints();
    bool isCoreObject(int idx);
    std::vector<std::vector<int>> getCluster();
};

#endif