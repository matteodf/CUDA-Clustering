
#include "common.h"
#include "kmeans_utils.h"
#include "msort/msort.h"
#include "scan/scan.h"

#include <string>

#define NUMBER_OF_POINTS (1024 * 1024)
#define NUMBER_OF_CLUSTERS 32
#define MAXIMUM_ITERATIONS 100
#define RADIUS 100.0

#define READ_FROM_FILE true
#define PRINT_RESULTS false

#define THREADSxBLOCK 512

__constant__ float const_cx[NUMBER_OF_CLUSTERS];
__constant__ float const_cy[NUMBER_OF_CLUSTERS];

__global__ void updateCentroids(POINT pts, int numPts, POINT centroids, int *counts, int numCentroids) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int i;

    extern __shared__ int sharedMem[];
    float *scx = (float *)sharedMem;
    float *scy = (float *)&scx[numCentroids];
    int *scg = (int *)&scy[numCentroids];

    if (tid < numCentroids) {
        scg[tid] = 0;
        scx[tid] = 0;
        scy[tid] = 0;

        int group_dim = (tid == (numCentroids - 1) ? (numPts - counts[tid]) : (counts[tid + 1] - counts[tid]));
        int group_start = counts[tid];
        scg[tid] = group_dim;

        for (i = group_start; i < (group_start + group_dim); i++) {
            scx[tid] += pts.x[i];
            scy[tid] += pts.y[i];
        }

        scx[tid] /= scg[tid];
        scy[tid] /= scg[tid];

        centroids.x[tid] = scx[tid];
        centroids.y[tid] = scy[tid];
        centroids.group[tid] = scg[tid];
    }
}

__global__ void updatePoints(int *d_changes, POINT pts, int numPts, POINT centroids, int numCentroids) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float d, min_d, px, py;
    int j, clusterIndex;

    min_d = HUGE_VAL;
    if (tid < numPts) {
        clusterIndex = pts.group[tid];
        px = pts.x[tid];
        py = pts.y[tid];
        for (j = 0; j < numCentroids; j++) {
            d = dist2(const_cx[j], const_cy[j], px, py);
            if (d < min_d) {
                min_d = d;
                clusterIndex = j;
            }
        }
        if (clusterIndex != pts.group[tid]) {
            pts.group[tid] = clusterIndex;
            atomicAdd(d_changes, 1);
        }
    }
}

__global__ void countElements(int *pts_group, int numPts, int *counts) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid < numPts) {
        atomicAdd(&counts[pts_group[tid]], 1);
    }
}

void calcPSum(POINT dev_pts, int numPts, int *&d_counts, int numCentroids) {
    CHECK(cudaMemset(d_counts, 0, sizeof(int) * numCentroids));

    countElements<<<(numPts / THREADSxBLOCK), THREADSxBLOCK>>>(dev_pts.group, numPts, d_counts);
    CHECK(cudaDeviceSynchronize());

    int *d_counts_cumul;
    CHECK(cudaMalloc((void **)&d_counts_cumul, NUMBER_OF_CLUSTERS * sizeof(int)));
    scan(d_counts_cumul, d_counts, numCentroids);

    CHECK(cudaFree(d_counts));
    d_counts = d_counts_cumul;
}

int kmeans(POINT pts, int numPts, POINT centroids, int numCentroids, int maxTimes) {
    int acceptable = numPts / 1000;

    POINT dev_pts;
    POINT dev_centroids;

    CHECK(cudaMalloc((void **)&(dev_pts.x), sizeof(float) * numPts));
    CHECK(cudaMalloc((void **)&(dev_pts.y), sizeof(float) * numPts));
    CHECK(cudaMalloc((void **)&(dev_pts.group), sizeof(int) * numPts));

    CHECK(cudaMalloc((void **)&(dev_centroids.x), sizeof(float) * numCentroids));
    CHECK(cudaMalloc((void **)&(dev_centroids.y), sizeof(float) * numCentroids));
    CHECK(cudaMalloc((void **)&(dev_centroids.group), sizeof(int) * numCentroids));

    int h_changes;
    int *d_changes;
    CHECK(cudaMalloc((void **)&d_changes, sizeof(int)));

    int threadsUpdtPts = THREADSxBLOCK;
    int blocksUpdtPts = ((numPts / threadsUpdtPts) == 0 ? 1 : (numPts / threadsUpdtPts));

    int threadsUpdtCds = THREADSxBLOCK;
    int blocksUpdtCds = ((numCentroids / threadsUpdtCds) == 0 ? 1 : (numPts / threadsUpdtCds));

    cudaMemcpyPS(dev_pts, pts, numPts, cudaMemcpyHostToDevice);
    cudaMemcpyPS(dev_centroids, centroids, numCentroids, cudaMemcpyHostToDevice);

    int *d_counts;
    CHECK(cudaMalloc((void **)&d_counts, numCentroids * sizeof(int)));

    int sharedMemSize = sizeof(float) * 3 * numCentroids;

    do {
        runMergeSort(dev_pts, numPts);
        calcPSum(dev_pts, numPts, d_counts, numCentroids);
        updateCentroids<<<blocksUpdtCds, threadsUpdtCds, sharedMemSize>>>(dev_pts, numPts, dev_centroids, d_counts, numCentroids);
        h_changes = 0;
        CHECK(cudaMemcpy(d_changes, &h_changes, sizeof(int), cudaMemcpyHostToDevice));

        CHECK(cudaMemcpyToSymbol(const_cx, dev_centroids.x, sizeof(float) * numCentroids, 0, cudaMemcpyDeviceToDevice));
        CHECK(cudaMemcpyToSymbol(const_cy, dev_centroids.y, sizeof(float) * numCentroids, 0, cudaMemcpyDeviceToDevice));

        updatePoints<<<blocksUpdtPts, threadsUpdtPts>>>(d_changes, dev_pts, numPts, dev_centroids, numCentroids);
        CHECK(cudaMemcpy(&h_changes, d_changes, sizeof(int), cudaMemcpyDeviceToHost));
        maxTimes--;
    } while ((h_changes > acceptable) && (maxTimes > 0));

    cudaMemcpyPS(pts, dev_pts, numPts, cudaMemcpyDeviceToHost);
    cudaMemcpyPS(centroids, dev_centroids, numCentroids, cudaMemcpyDeviceToHost);

    for (int i = 0; i < numCentroids; i++)
        centroids.group[i] = i;

    cudaFreePS(dev_pts);
    cudaFreePS(dev_centroids);
    CHECK(cudaDeviceReset());
    return (MAXIMUM_ITERATIONS - maxTimes);
}

void test(char* inputFile, bool initRand,int testIter){
    int numPts, numCentroids, maxTimes, numIter;
    double start, stop, timeInit, timeKmeans;
    
    numPts = findNumRows(inputFile);
    numCentroids = NUMBER_OF_CLUSTERS;
    maxTimes = MAXIMUM_ITERATIONS;

    POINT pts, centroids;
    mallocPS(&pts, numPts);
    mallocPS(&centroids, numCentroids);

    if (numCentroids == 1 || numPts <= 0 || numCentroids > numPts) {
        printf("Error, wrong parameters\n");
        exit(1);
    }
    if (maxTimes < 1) maxTimes = 1;
    
    readPointsFromFile(&pts, inputFile);
    
    for (int i=0; i < testIter; i++){
        start = seconds();
        if (initRand) initClustersRandom(pts, numPts, numCentroids);
        if (!initRand) initClusters(pts, numPts, centroids, numCentroids);
        stop = seconds();
        timeInit = stop - start;

        start = seconds();
        numIter = kmeans(pts, numPts, centroids, numCentroids, maxTimes);
        stop = seconds();
        timeKmeans = stop - start;

        printf("%f %f %d\n", timeInit, timeKmeans, numIter);
    }

    char outputFile[] = "../../../../output/output_impl2.txt";
    writePointsToFile(outputFile, pts, numPts, centroids, numCentroids);

    freePS(&pts);
    freePS(&centroids);
}

int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Please follow this format: ./app [intputFile] [initRand = t / f] [numTestIterations]\n");
        return 0;
    }
    //char filename[] = "../../../../input/input_rand2.txt";
    char *filename = argv[1];
    bool initRand = (((std::string)argv[2]) == (std::string)"t");
    int testIter = std::stoi(argv[3]);

    test(filename, initRand, testIter);
    return 0;
}