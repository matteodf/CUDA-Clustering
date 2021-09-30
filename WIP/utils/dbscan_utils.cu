#include <sys/time.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <sys/stat.h>
//#include <sys/types.h>
//#include <algorithm>
//#include <map>
//#include <string>
//include <bits/stdc++.h>

using namespace std;

const int NOISE = -2;
const int NOT_CLASSIFIED = -1;

class Point {
   public:
    float x, y;
    int ptsCnt, cluster;
    vector<int> adjPoints;
    float getDis(const Point& ot) {
        return sqrt((x - ot.x) * (x - ot.x) + (y - ot.y) * (y - ot.y));
    }
};

class DBSCAN {
   public:
    int minPts;
    float eps;
    vector<Point> points;
    int size;
    vector<bool> visited;
    vector<vector<int> > cluster;
    int clusterIdx;

    DBSCAN(float eps, int minPts, vector<Point> points) {
        this->eps = eps;
        this->minPts = minPts;
        this->points = points;
        this->size = (int)points.size();
        this->clusterIdx = -1;
    }
    void run() {

        int nBytesPts = sizeof(Point) * size;
        Point* dev_pts;
        int threads = 32;
        int blocks = size / threads;

        CHECK(cudaMalloc((void**) &dev_pts, nBytesPts));
        
        CHECK(cudaMemcpy(dev_pts, &points[0], nBytesPts, cudaMemcpyHostToDevice));

        checkNearPoints<<<block, threads>>>();
        classifyPoints<<<block, threads>>>();

        for (int i = 0; i < size; i++) {
            if (points[i].cluster != NOT_CLASSIFIED) continue;

            if (isCoreObject(i)) {
                dfs(i, ++clusterIdx);
            } else {
                points[i].cluster = NOISE;
            }
        }

        cluster.resize(clusterIdx + 1);
        for (int i = 0; i < size; i++) {
            if (points[i].cluster != NOISE) {
                cluster[points[i].cluster].push_back(i);
            }
        }
    };

    __global__ void classifyPoints(){
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (dev_pts[tid].cluster != NOT_CLASSIFIED) continue;
        if (dev_pts[tid].ptsCnt >= minPts;){

        } else {
            dev_pts[tid].cluster = NOISE            
        }
    };

    void dfs(int now, int c) {
        points[now].cluster = c;
        if (!isCoreObject(now)) return;

        for (auto& next : adjPoints[now]) {
            if (points[next].cluster != NOT_CLASSIFIED && points[next].cluster != NOISE) continue;
            dfs(next, c);
        }
    }

    __global__ void checkNearPoints() {
        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        for (int j = 0; j < size; j++) {
            if (tid == j) continue;
            if (dev_pts[tid].getDis(dev_pts[j]) <= eps) {
                dev_pts[tid].ptsCnt++;
                dev_pts[tid].adjPoints.push_back(j);
            }
        }
    }
    // is idx'th point core object?
    bool isCoreObject(int idx) {
        return points[idx].ptsCnt >= minPts;
    }

    vector<vector<int> > getCluster() {
        return cluster;
    }
};

class InputReader {
   private:
    ifstream fin;
    vector<Point> points;

   public:
    InputReader(string filename) {
        filename = "../../input_data/" + filename;
        fin.open(filename);
        if (!fin) {
            cout << filename << " file could not be opened\n";
            exit(0);
        }
        parse();
    }
    void parse() {
        int idx;
        float x, y;
        while (!fin.eof()) {
            fin >> idx >> x >> y;
            points.push_back({x, y, 0, NOT_CLASSIFIED});
        }
        points.pop_back();
    }
    vector<Point> getPoints() {
        return points;
    }
};

double seconds();
void printResults(DBSCAN dbScan, string filename);