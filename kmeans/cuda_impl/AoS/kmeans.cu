
#define NUMBER_OF_POINTS 1024 * 1024
#define NUMBER_OF_CLUSTERS 32
#define MAXIMUM_ITERATIONS 100
#define RADIUS 100.0

#define THREADSxBLOCK 512

#include "common.h"
#include "kmeans_utils.h"

__constant__ POINT const_cnt[NUMBER_OF_CLUSTERS];

void printDebug(POINT* pts, int num_pts, POINT* centroids, int num_clusters) {
    printf("centroids: -----------------\n");
    for (int i = 0; i < num_clusters; i++) {
        printf("\t%f %f, %d\n", centroids[i].x, centroids[i].y, centroids[i].group);
    }
}

__global__ void addUpClusters(POINT* pts, int numPts, POINT* centroids) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < numPts) {
        int g = pts[tid].group;
        atomicAdd(&centroids[g].group, 1);
        atomicAdd(&centroids[g].x, pts[tid].x);
        atomicAdd(&centroids[g].y, pts[tid].y);
    }
}

__global__ void computeMean(POINT* centroids, int numCentroids) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < numCentroids) {
        centroids[tid].x /= centroids[tid].group;
        centroids[tid].y /= centroids[tid].group;
    }
}

void updateCentroids(POINT* d_pts, int numPts, POINT* d_centroids, int numCentroids) {
    int threads = THREADSxBLOCK;
    int blocks = (numPts + threads - 1) / threads;
    CHECK(cudaMemset(d_centroids, 0, sizeof(POINT) * numCentroids));
    addUpClusters<<<blocks, threads>>>(d_pts, numPts, d_centroids);
    CHECK(cudaDeviceSynchronize());
    computeMean<<<(numCentroids + threads - 1) / threads, threads>>>(d_centroids, numCentroids);
    CHECK(cudaDeviceSynchronize());
}

__global__ void updatePoints(int* d_changes, POINT* pts, int num_pts, POINT* centroids, int num_clusters) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float d, min_d;
    int j, clusterIndex;
    POINT p;

    min_d = HUGE_VAL;
    if (tid < num_pts) {
        clusterIndex = pts[tid].group;
        p = pts[tid];
        for (j = 0; j < num_clusters; j++) {
            d = dist2(&const_cnt[j], &p);
            if (d < min_d) {
                min_d = d;
                clusterIndex = j;
            }
        }
        if (clusterIndex != pts[tid].group) {
            pts[tid].group = clusterIndex;
            atomicAdd(d_changes, 1);
        }
    }
}

void kmeans(POINT* pts, int num_pts, POINT* centroids, int num_clusters, int maxTimes) {
    int acceptable = num_pts / 1000;

    int nBytesPts = sizeof(POINT) * num_pts;
    int nBytesCentroids = sizeof(POINT) * num_clusters;

    POINT* dev_pts;
    POINT* dev_centroids;
    CHECK(cudaMalloc((void**)&dev_pts, nBytesPts));
    CHECK(cudaMalloc((void**)&dev_centroids, nBytesCentroids));

    int h_changes;
    int* d_changes;
    CHECK(cudaMalloc((void**)&d_changes, sizeof(int)));

    int threads = THREADSxBLOCK;
    int blocks = (num_pts + threads - 1) / threads;

    CHECK(cudaMemcpy(dev_pts, pts, nBytesPts, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dev_centroids, centroids, nBytesCentroids, cudaMemcpyHostToDevice));

    do {
        updateCentroids(dev_pts, num_pts, dev_centroids, num_clusters);

        CHECK(cudaMemcpyToSymbol(const_cnt, dev_centroids, nBytesCentroids, 0, cudaMemcpyDeviceToDevice));

        h_changes = 0;
        CHECK(cudaMemcpy(d_changes, &h_changes, sizeof(int), cudaMemcpyHostToDevice));
        updatePoints<<<blocks, threads>>>(d_changes, dev_pts, num_pts, dev_centroids, num_clusters);
        CHECK(cudaMemcpy(&h_changes, d_changes, sizeof(int), cudaMemcpyDeviceToHost));
        maxTimes--;
    } while ((h_changes > acceptable) && (maxTimes > 0));

    CHECK(cudaMemcpy(pts, dev_pts, nBytesPts, cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(centroids, dev_centroids, nBytesCentroids, cudaMemcpyDeviceToHost));

    for (int i = 0; i < num_clusters; i++)
        centroids[i].group = i;

    CHECK(cudaFree(dev_pts));
    CHECK(cudaFree(dev_centroids));
}

int main() {
    int num_pts = NUMBER_OF_POINTS;
    int num_clusters = NUMBER_OF_CLUSTERS;
    int maxTimes = MAXIMUM_ITERATIONS;
    float radius = RADIUS;
    double start, time;

    if (num_clusters == 1 || num_pts <= 0 || num_clusters > num_pts) {
        printf("Error\n");
        return 1;
    }
    if (maxTimes < 1)
        maxTimes = 1;

    POINT* pts = (POINT*)malloc(sizeof(POINT) * num_pts);
    POINT* centroids = (POINT*)malloc(sizeof(POINT) * num_clusters);

    gen_xy(pts, num_pts, radius);
    initClusters(pts, num_pts, centroids, num_clusters);

    start = seconds();
    kmeans(pts, num_pts, centroids, num_clusters, maxTimes);
    time = seconds() - start;
    printf("kmeans total time: %f (sec)\n", time);

    //printDebug(pts, num_pts, centroids, num_clusters);
    //printResult(pts, num_pts, centroids, num_clusters);

    free(pts);
    free(centroids);
    return 0;
}