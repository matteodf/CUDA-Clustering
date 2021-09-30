#include "dbscan.h"

DBSCAN::DBSCAN(double eps, int minPts, Points points) {
    this->eps = eps;
    this->minPts = minPts;
    this->points = points;
    this->size = points.getSize();
    adjOffset = minPts + log2(size);                             //questa è la dimensione predefinita per le adiacenze di diciascun punto
    adjPoints = (int*)malloc(sizeof(int) * (size * adjOffset));  // l'array di adiacenze linearizzato
    this->clusterIdx = -1;
}
void DBSCAN::run() {
    double start, stop;

    start = seconds();
    checkNearPoints();
    stop = seconds();
    if (!PRINT_RESULTS) printf("checkNearPoints time: %f (sec)\n", stop - start);

    start = seconds();
    for (int i = 0; i < size; i++) {
        if (points.cluster[i] != NOT_CLASSIFIED) continue;

        if (isCoreObject(i)) {
            dfs(i, ++clusterIdx);
        } else {
            points.cluster[i] = NOISE;
        }
    }
    stop = seconds();
    if (!PRINT_RESULTS) printf("loop with dfs call time: %f (sec)\n", stop - start);

    start = seconds();
    cluster.resize(clusterIdx + 1);
    n = clusterIdx + 1;
    for (int i = 0; i < size; i++) {
        if (points.cluster[i] != NOISE) {
            cluster[points.cluster[i]].push_back(i);
        }
    }
    stop = seconds();
    if (!PRINT_RESULTS) printf("last loop time: %f (sec)\n", stop - start);
}

void DBSCAN::dfs(int now, int c) {
    points.cluster[now] = c;
    if (!isCoreObject(now)) return;

    //auto& next : adjPoints[now]
    for (int i = 0; i < points.ptsCnt[now]; i++) {
        int next = adjPoints[now * adjOffset + i];
        if (points.cluster[next] != NOT_CLASSIFIED && points.cluster[next] != NOISE) continue;
        dfs(next, c);
    }
}

void DBSCAN::checkNearPoints() {
    int k;
    for (int i = 0; i < size; i++) {
        k = 0;
        for (int j = 0; j < size; j++) {
            if (i == j) continue;
            if (dist(points.x[i], points.y[i], points.x[j], points.y[j]) <= eps) {
                points.ptsCnt[i]++;
                if (k >= adjOffset) {
                    printf("ERROR too many adjPoints\n");
                    exit(1);
                }
                adjPoints[i * adjOffset + k] = j;
                k++;
            }
        }
    }
}
// is idx'th point core object?
bool DBSCAN::isCoreObject(int idx) {
    return points.ptsCnt[idx] >= minPts;
}

std::vector<std::vector<int>> DBSCAN::getCluster() {
    return cluster;
}