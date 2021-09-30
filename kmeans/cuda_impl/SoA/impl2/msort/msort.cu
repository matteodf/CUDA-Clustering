#include "../common.h"
#include "../kmeans_utils.h"

void swap(int *&a, int *&b) {
    int *temp = a;
    a = b;
    b = temp;
}
void swap(float *&a, float *&b) {
    float *temp = a;
    a = b;
    b = temp;
}

void check_up_sorting(int array[], unsigned size) {
    bool flag = true;
    for (int i = 0; i < size - 1; i++)
        if (array[i] > array[i + 1]) {
            printf("Sorting error! array[%d]=%d array[%d]=%d\n", i, array[i], i + 1,
                   array[i + 1]);
            flag = false;
            break;
        }
    if (flag)
        printf("   Sorting OK!\n");
}

__device__ int binarySearch(int arr[], int x, int k, bool UP) {
    int l = 0, r = k;

    while (l < r) {
        int m = (l + r) / 2;
        if (UP) {  // for upper chunk B
            if (arr[m] <= x)
                l = m + 1;
            else
                r = m;
        } else {  // for lower chunk A
            if (arr[m] < x)
                l = m + 1;
            else
                r = m;
        }
    }

    return l;
}

__global__ void cudaMergeSortMulti(POINT ps, POINT sorted, int n, int k) {
    // k = 1,2,4,8,16,..., 2^m chunk dims
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int j = tid % k;
    int l = (tid - j) * 2;  // first element of the fisrt chunk
    int i = l + j;          // A[i] first chunk   [][][*][] and  B[i+k] [][][*][]

    if (k == 1) {
        l = 2 * tid;
        i = l;
    }

    // find the relative position of x within B[*]
    int x = ps.group[i];
    int p = binarySearch(ps.group + l + k, x, k, 1);
    sorted.group[i + p] = x;
    sorted.x[i + p] = ps.x[i];
    sorted.y[i + p] = ps.y[i];

    // find the relative position of y within A[*]
    int y = ps.group[i + k];
    p = binarySearch(ps.group + l, y, k, 0);
    sorted.group[i + p] = y;
    sorted.x[i + p] = ps.x[i + k];
    sorted.y[i + p] = ps.y[i + k];
}

void runMergeSort(POINT &dev_pts, int N) {
    int BLOCK_SIZE = 32;
    bool array2sorted = false;
    int nThreads = N / 2;

    POINT d_sorted;
    CHECK(cudaMalloc((void **)&(d_sorted.x), sizeof(float) * N));
    CHECK(cudaMalloc((void **)&(d_sorted.y), sizeof(float) * N));
    CHECK(cudaMalloc((void **)&(d_sorted.group), sizeof(int) * N));

    dim3 block(min(nThreads, BLOCK_SIZE));
    dim3 grid((nThreads + block.x - 1) / block.x);

    for (int chunk = 1; chunk <= N / 2; chunk *= 2) {
        array2sorted = !array2sorted;
        if (array2sorted)
            cudaMergeSortMulti<<<grid, block>>>(dev_pts, d_sorted, N, chunk);
        else
            cudaMergeSortMulti<<<grid, block>>>(d_sorted, dev_pts, N, chunk);
    }
    CHECK(cudaDeviceSynchronize());

    if (array2sorted) {
        swap(dev_pts.x, d_sorted.x);
        swap(dev_pts.y, d_sorted.y);
        swap(dev_pts.group, d_sorted.group);
    }
}