
#include "kmeans_utils.h"

#include <string>

#define NUMBER_OF_POINTS (1024 * 1024)
#define NUMBER_OF_CLUSTERS 32
#define MAXIMUM_ITERATIONS 100
#define RADIUS 100.0

#define READ_FROM_FILE true
#define PRINT_RESULTS true

int kmeans(POINT* pts, int num_pts, POINT* centroids, int num_clusters, int maxTimes) {
    int i, clusterIndex;
    int changes;
    int acceptable = num_pts / 1000;

    do {
        for (i = 0; i < num_clusters; i++) {
            centroids[i].group = 0;
            centroids[i].x = 0;
            centroids[i].y = 0;
        }

        for (i = 0; i < num_pts; i++) {
            clusterIndex = pts[i].group;
            centroids[clusterIndex].group++;
            centroids[clusterIndex].x += pts[i].x;
            centroids[clusterIndex].y += pts[i].y;
        }

        for (i = 0; i < num_clusters; i++) {
            centroids[i].x /= centroids[i].group;
            centroids[i].y /= centroids[i].group;
        }

        changes = 0;
        for (i = 0; i < num_pts; i++) {
            clusterIndex = nearest(&pts[i], centroids, num_clusters);
            if (clusterIndex != pts[i].group) {
                pts[i].group = clusterIndex;
                changes++;
            }
        }

        maxTimes--;
    } while ((changes > acceptable) && (maxTimes > 0));

    for (i = 0; i < num_clusters; i++)
        centroids[i].group = i;
    return MAXIMUM_ITERATIONS - maxTimes;
}

void test(char* inputFile, bool initRand, int testIter){
    int numPts, numCentroids, maxTimes, numIter;
    double start, stop, timeInit, timeKmeans;
    
    numPts = findNumRows(inputFile);
    numCentroids = NUMBER_OF_CLUSTERS;
    maxTimes = MAXIMUM_ITERATIONS;

    POINT* pts = (POINT*)malloc(sizeof(POINT) * numPts);
    POINT* centroids = (POINT*)malloc(sizeof(POINT) * numCentroids);

    if (numCentroids == 1 || numPts <= 0 || numCentroids > numPts) {
        printf("Error, wrong parameters\n");
        exit(1);
    }
    if (maxTimes < 1) maxTimes = 1;
    
    readPointsFromFile(pts, inputFile);
    
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

    char outputFile[] = "../../output/output_seq.txt";
    writePointsToFile(outputFile, pts, numPts, centroids, numCentroids);

    free(pts);
    free(centroids);
}

int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Please follow this format: ./app [intputFile] [initRand = t / f] [numTestIterations]\n");
        return 0;
    }
    //char filename[] = "../../input/input_rand2.txt";
    char *filename = argv[1];
    bool initRand = (((std::string)argv[2]) == (std::string)"t");
    int testIter = std::stoi(argv[3]);

    test(filename, initRand, testIter);
    return 0;
}