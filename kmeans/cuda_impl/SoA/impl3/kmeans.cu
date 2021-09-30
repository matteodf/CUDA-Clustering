#include "common.h"
#include "kmeans_utils.h"
#include "scan/scan.h"

#define NUMBER_OF_POINTS (1024 * 1024)
#define NUMBER_OF_CLUSTERS 32
#define MAXIMUM_ITERATIONS 100
#define RADIUS 100.0

#define THREADSxBLOCK 32

void cudaMemcpyPS(POINT dst, POINT src, int n, enum cudaMemcpyKind kind) {
    CHECK(cudaMemcpy(dst.x, src.x, sizeof(float) * n, kind));
    CHECK(cudaMemcpy(dst.y, src.y, sizeof(float) * n, kind));
    CHECK(cudaMemcpy(dst.group, src.group, sizeof(int) * n, kind));
}

void cudaFreePS(POINT ps) {
    CHECK(cudaFree(ps.x));
    CHECK(cudaFree(ps.y));
    CHECK(cudaFree(ps.group));
}

__global__ void cudaPrintArray(int *arr, int n) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid == 0) {
        for (int i = 0; i < n; i++) {
            printf("%d\t", arr[i]);
        }
        printf("\n");
    }
}
__global__ void cudaPrintArray(float *arr, int n) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid == 0) {
        for (int i = 0; i < n; i++) {
            printf("%.2f\t", arr[i]);
        }
        printf("\n");
    }
}

__global__ void updateCentroids(POINT pts, int num_pts, POINT centroids, int num_clusters) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int i;

    extern __shared__ int sharedMem[];
    float *scx = (float *)sharedMem;
    float *scy = (float *)&scx[num_clusters];
    int *scg = (int *)&scy[num_clusters];

    if (tid < num_clusters) {
        scg[tid] = 0;
        scx[tid] = 0;
        scy[tid] = 0;

        for (i = 0; i < num_pts; i++) {
            if (pts.group[i] == tid) {
                scg[tid]++;
                scx[tid] += pts.x[i];
                scy[tid] += pts.y[i];
            }
        }

        scx[tid] /= scg[tid];
        scy[tid] /= scg[tid];

        centroids.x[tid] = scx[tid];
        centroids.y[tid] = scy[tid];
        centroids.group[tid] = scg[tid];
    }
}

__global__ void updatePoints(int *d_changes, POINT pts, int num_pts, POINT centroids, int num_clusters) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float d, min_d;
    int j, clusterIndex;

    //TODO usare sharedMem per centroids?

    min_d = HUGE_VAL;
    if (tid < num_pts) {
        clusterIndex = pts.group[tid];
        for (j = 0; j < num_clusters; j++) {
            d = dist2(centroids.x[j], centroids.y[j], pts.x[tid], pts.y[tid]);
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

__global__ void countElements(int *pts_group, int num_pts, int *counts) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid < num_pts) {
        atomicAdd(&counts[pts_group[tid]], 1);
    }
}

__global__ void findIndexes(int *pts_group, float *indexes, int num_pts, int *counts, int *counts_cumul, int num_clusters) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid < num_pts) {
        int group = pts_group[tid];
        //int group_dim = (tid == (num_pts - 1) ? (num_pts - counts_cumul[group]) : (counts_cumul[group + 1] - counts_cumul[group]));
        int i = counts_cumul[group] + atomicAdd(&counts[group], -1) - 1;
        indexes[i] = tid;
    }
}

__global__ void sortGroup(int *pts_group, int num_pts, int *counts, int num_clusters) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int high = (tid == num_clusters - 1) ? num_pts : counts[tid + 1];
    for (int i = counts[tid]; i < high; i++) {
        pts_group[i] = tid;
    }
}

__global__ void sortXY(POINT pts, float *indexes, int num_pts) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid < num_pts) {
        // per fare lo swap dei valori di x e y viene sfruttato l'array indexes per salvare
        // i valori ordinati di x e x viene usato per salvare i valori ordinati di y
        int index = indexes[tid];
        // // __syncthreads();  //deve essere su tutta la griglia
        indexes[index] = pts.x[tid];
        // // __syncthreads();  //deve essere su tutta la griglia
        pts.x[index] = pts.y[tid];
    }
}

void countingSort(POINT dev_pts, int num_pts, int *&d_counts, int num_clusters) {
    //step 1
    CHECK(cudaMemset(d_counts, 0, sizeof(int) * num_clusters));

    cudaPrintArray<<<1, 1>>>(dev_pts.group, num_pts);
    CHECK(cudaDeviceSynchronize());
    cudaPrintArray<<<1, 1>>>(dev_pts.x, num_pts);
    CHECK(cudaDeviceSynchronize());
    printf("\n");

    //step 2
    countElements<<<(num_pts / THREADSxBLOCK), THREADSxBLOCK>>>(dev_pts.group, num_pts, d_counts);
    CHECK(cudaDeviceSynchronize());

    //step 3
    int *d_counts_cumul;
    CHECK(cudaMalloc((void **)&d_counts_cumul, NUMBER_OF_CLUSTERS * sizeof(int)));
    scan(d_counts_cumul, d_counts, num_clusters);

    //step 3.1 -> alloca e riempi array degli indici per pts.x e pts.y
    float *d_indexes;
    CHECK(cudaMalloc((void **)&d_indexes, num_pts * sizeof(float)));
    findIndexes<<<(num_pts / THREADSxBLOCK), THREADSxBLOCK>>>(dev_pts.group, d_indexes, num_pts, d_counts, d_counts_cumul, num_clusters);
    CHECK(cudaDeviceSynchronize());

    //step 3.2
    CHECK(cudaFree(d_counts));
    d_counts = d_counts_cumul;

    //step 4.1
    sortGroup<<<(num_clusters / THREADSxBLOCK), THREADSxBLOCK>>>(dev_pts.group, num_pts, d_counts, num_clusters);
    CHECK(cudaDeviceSynchronize());

    // cudaPrintArray<<<1, 1>>>(dev_pts.group, num_pts);
    // CHECK(cudaDeviceSynchronize());
    // printf("\n");

    //step 4.2
    sortXY<<<(num_pts / THREADSxBLOCK), THREADSxBLOCK>>>(dev_pts, d_indexes, num_pts);
    CHECK(cudaDeviceSynchronize());

    //step 4.3 -> sistema gli array che sono stati spostati nel precedente passaggio
    // dev_pts.y = dev_pts.x;
    // CHECK(cudaFree(dev_pts.x));
    // dev_pts.x = d_indexes;

    cudaPrintArray<<<1, 1>>>(dev_pts.group, num_pts);
    CHECK(cudaDeviceSynchronize());
    cudaPrintArray<<<1, 1>>>(d_indexes, num_pts);
    CHECK(cudaDeviceSynchronize());
    printf("\n");
}

void kmeans(POINT pts, int num_pts, POINT centroids, int num_clusters, int maxTimes) {
    int acceptable = num_pts / 1000;

    // for (int i = 0; i < num_pts; i++)
    //     printf("(%d, %.2f)\t", pts.group[i], pts.x[i]);
    // printf("\n");
    // printf("\n");

    // for (int i = 0; i < num_pts; i++)
    //     printf("%d\t", pts.group[i]);
    // printf("\n");
    // printf("\n");

    POINT dev_pts;
    POINT dev_centroids;

    CHECK(cudaMalloc((void **)&(dev_pts.x), sizeof(float) * num_pts));
    CHECK(cudaMalloc((void **)&(dev_pts.y), sizeof(float) * num_pts));
    CHECK(cudaMalloc((void **)&(dev_pts.group), sizeof(int) * num_pts));

    CHECK(cudaMalloc((void **)&(dev_centroids.x), sizeof(float) * num_clusters));
    CHECK(cudaMalloc((void **)&(dev_centroids.y), sizeof(float) * num_clusters));
    CHECK(cudaMalloc((void **)&(dev_centroids.group), sizeof(int) * num_clusters));

    int h_changes;
    int *d_changes;
    CHECK(cudaMalloc((void **)&d_changes, sizeof(int)));

    int threads = THREADSxBLOCK;
    int blocks = num_pts / threads;

    cudaMemcpyPS(dev_pts, pts, num_pts, cudaMemcpyHostToDevice);
    cudaMemcpyPS(dev_centroids, centroids, num_clusters, cudaMemcpyHostToDevice);

    int *d_counts;
    CHECK(cudaMalloc((void **)&d_counts, num_clusters * sizeof(int)));

    int sharedMemSize = sizeof(float) * 3 * num_clusters;

    do {
        updateCentroids<<<(num_clusters / threads), threads, sharedMemSize>>>(dev_pts, num_pts, dev_centroids, num_clusters);
        h_changes = 0;
        CHECK(cudaMemcpy(d_changes, &h_changes, sizeof(int), cudaMemcpyHostToDevice));
        updatePoints<<<blocks, threads>>>(d_changes, dev_pts, num_pts, dev_centroids, num_clusters);
        CHECK(cudaMemcpy(&h_changes, d_changes, sizeof(int), cudaMemcpyDeviceToHost));
        maxTimes--;
    } while ((h_changes > acceptable) && (maxTimes > 0));

    cudaMemcpyPS(pts, dev_pts, num_pts, cudaMemcpyDeviceToHost);
    cudaMemcpyPS(centroids, dev_centroids, num_clusters, cudaMemcpyDeviceToHost);

    // for (int i = 0; i < num_pts; i++)
    //     printf("(%d, %.2f)\t", pts.group[i], pts.x[i]);
    // printf("\n");
    // printf("\n");

    // for (int i = 0; i < num_pts; i++)
    //     printf("%d\t", pts.group[i]);
    // printf("\n");
    // printf("\n");

    for (int i = 0; i < num_clusters; i++)
        centroids.group[i] = i;

    cudaFreePS(dev_pts);
    cudaFreePS(dev_centroids);
}

int main() {
    int num_pts = NUMBER_OF_POINTS;
    int num_clusters = NUMBER_OF_CLUSTERS;
    int maxTimes = MAXIMUM_ITERATIONS;
    float radius = RADIUS;
    double start, time;

    if (num_clusters == 1 || num_pts <= 0 || num_clusters > num_pts) {
        printf("Error, wrong parameters\n");
        return 1;
    }
    if (maxTimes < 1)
        maxTimes = 1;

    POINT pts, centroids;
    mallocPS(&pts, num_pts);
    mallocPS(&centroids, num_clusters);

    gen_xy(pts, num_pts, radius);
    initClusters(pts, num_pts, centroids, num_clusters);

    start = seconds();
    kmeans(pts, num_pts, centroids, num_clusters, maxTimes);
    time = seconds() - start;
    printf("kmeans total time: %f (sec)\n", time);

    printDebug(centroids, num_clusters);
    //printResult(pts, num_pts, centroids, num_clusters);

    freePS(&pts);
    freePS(&centroids);
    return 0;
}