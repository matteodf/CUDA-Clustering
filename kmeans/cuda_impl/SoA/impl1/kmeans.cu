
#include "kmeans_utils.h"
#include <string>

#define NUMBER_OF_POINTS (1024*1024)
#define NUMBER_OF_CLUSTERS 32
#define MAXIMUM_ITERATIONS 100
#define RADIUS 100.0

#define PRINT_RESULTS false
#define READ_FROM_FILE true

#define THREADSxBLOCK 512

__constant__ float const_cx[NUMBER_OF_CLUSTERS];
__constant__ float const_cy[NUMBER_OF_CLUSTERS];

__global__ void addUpClusters(POINT pts, int numPts, POINT centroids) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < numPts) {
        int g = pts.group[tid];
        atomicAdd(&centroids.group[g], 1);
        atomicAdd(&centroids.x[g], pts.x[tid]);
        atomicAdd(&centroids.y[g], pts.y[tid]);
    }
}

__global__ void computeMean(POINT centroids, int numCentroids) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < numCentroids) {
        centroids.x[tid] /= centroids.group[tid];
        centroids.y[tid] /= centroids.group[tid];
    }
}

void updateCentroids(POINT d_pts, int numPts, POINT d_centroids, int numCentroids) {
    int threads = THREADSxBLOCK;
    int blocks = (numPts + threads - 1) / threads;
    CHECK(cudaMemset(d_centroids.group, 0, sizeof(int) * numCentroids));
    CHECK(cudaMemset(d_centroids.x, 0, sizeof(float) * numCentroids));
    CHECK(cudaMemset(d_centroids.y, 0, sizeof(float) * numCentroids));
    addUpClusters<<<blocks, threads>>>(d_pts, numPts, d_centroids);
    CHECK(cudaDeviceSynchronize());
    computeMean<<<(numCentroids + threads - 1) / threads, threads>>>(d_centroids, numCentroids);
    CHECK(cudaDeviceSynchronize());
}

__global__ void updatePoints(int *changes, POINT pts, int numPts, int numCentroids) {
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
            atomicAdd(changes, 1);
        }
    }
}

int kmeans(POINT pts, int numPts, POINT centroids, int numCentroids, int maxTimes) {
    int acceptable = numPts / 1000;

    POINT d_pts;
    POINT d_centroids;

    CHECK(cudaMalloc((void **)&(d_pts.x), sizeof(float) * numPts));
    CHECK(cudaMalloc((void **)&(d_pts.y), sizeof(float) * numPts));
    CHECK(cudaMalloc((void **)&(d_pts.group), sizeof(int) * numPts));

    CHECK(cudaMalloc((void **)&(d_centroids.x), sizeof(float) * numCentroids));
    CHECK(cudaMalloc((void **)&(d_centroids.y), sizeof(float) * numCentroids));
    CHECK(cudaMalloc((void **)&(d_centroids.group), sizeof(int) * numCentroids));

    int h_changes;
    int *d_changes;
    CHECK(cudaMalloc((void **)&d_changes, sizeof(int)));

    int threads = THREADSxBLOCK;
    int blocks = (numPts + threads - 1) / threads;

    cudaMemcpyPS(d_pts, pts, numPts, cudaMemcpyHostToDevice);
    cudaMemcpyPS(d_centroids, centroids, numCentroids, cudaMemcpyHostToDevice);

    do {
        updateCentroids(d_pts, numPts, d_centroids, numCentroids);
        h_changes = 0;
        CHECK(cudaMemcpy(d_changes, &h_changes, sizeof(int), cudaMemcpyHostToDevice));

        CHECK(cudaMemcpyToSymbol(const_cx, d_centroids.x, sizeof(int) * numCentroids, 0, cudaMemcpyDeviceToDevice));
        CHECK(cudaMemcpyToSymbol(const_cy, d_centroids.y, sizeof(int) * numCentroids, 0, cudaMemcpyDeviceToDevice));

        updatePoints<<<blocks, threads>>>(d_changes, d_pts, numPts, numCentroids);
        CHECK(cudaMemcpy(&h_changes, d_changes, sizeof(int), cudaMemcpyDeviceToHost));
        maxTimes--;
    } while ((h_changes > acceptable) && (maxTimes > 0));

    cudaMemcpyPS(pts, d_pts, numPts, cudaMemcpyDeviceToHost);
    cudaMemcpyPS(centroids, d_centroids, numCentroids, cudaMemcpyDeviceToHost);

    for (int i = 0; i < numCentroids; i++)
        centroids.group[i] = i;

    cudaFreePS(d_pts);
    cudaFreePS(d_centroids);
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

    char outputFile[] = "../../../../output/output_impl1.txt";
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