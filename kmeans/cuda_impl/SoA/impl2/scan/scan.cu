#include "scan.h"

void scan(int *d_out, int *d_in, int length) {
    if (length > ELEMENTS_PER_BLOCK) {
        scanLargeDeviceArray(d_out, d_in, length);
    } else {
        scanSmallDeviceArray(d_out, d_in, length);
    }

    return;
}

void scanLargeDeviceArray(int *d_out, int *d_in, int length) {
    int remainder = length % (ELEMENTS_PER_BLOCK);
    if (remainder == 0) {
        scanLargeEvenDeviceArray(d_out, d_in, length);
    } else {
        // perform a large scan on a compatible multiple of elements
        int lengthMultiple = length - remainder;
        scanLargeEvenDeviceArray(d_out, d_in, lengthMultiple);

        // scan the remaining elements and add the (inclusive) last element of the large scan to this
        int *startOfOutputArray = &(d_out[lengthMultiple]);
        scanSmallDeviceArray(startOfOutputArray, &(d_in[lengthMultiple]), remainder);

        add<<<1, remainder>>>(startOfOutputArray, remainder, &(d_in[lengthMultiple - 1]), &(d_out[lengthMultiple - 1]));
    }
}

void scanSmallDeviceArray(int *d_out, int *d_in, int length) {
    int powerOfTwo = nextPowerOfTwo(length);
    prescan_arbitrary<<<1, (length + 1) / 2, 2 * powerOfTwo * sizeof(int)>>>(d_out, d_in, length, powerOfTwo);
}

void scanLargeEvenDeviceArray(int *d_out, int *d_in, int length) {
    const int blocks = length / ELEMENTS_PER_BLOCK;
    const int sharedMemArraySize = ELEMENTS_PER_BLOCK * sizeof(int);

    int *d_sums, *d_incr;
    cudaMalloc((void **)&d_sums, blocks * sizeof(int));
    cudaMalloc((void **)&d_incr, blocks * sizeof(int));

    prescan_large<<<blocks, THREADS_PER_BLOCK, 2 * sharedMemArraySize>>>(d_out, d_in, ELEMENTS_PER_BLOCK, d_sums);

    const int sumsArrThreadsNeeded = (blocks + 1) / 2;
    if (sumsArrThreadsNeeded > THREADS_PER_BLOCK) {
        // perform a large scan on the sums arr
        scanLargeDeviceArray(d_incr, d_sums, blocks);
    } else {
        // only need one block to scan sums arr so can use small scan
        scanSmallDeviceArray(d_incr, d_sums, blocks);
    }

    add<<<blocks, ELEMENTS_PER_BLOCK>>>(d_out, ELEMENTS_PER_BLOCK, d_incr);

    cudaFree(d_sums);
    cudaFree(d_incr);
}

__global__ void prescan_arbitrary(int *output, int *input, int n, int powerOfTwo) {
    extern __shared__ int temp[];  // allocated on invocation
    int threadID = threadIdx.x;

    int ai = threadID;
    int bi = threadID + (n / 2);
    int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
    int bankOffsetB = CONFLICT_FREE_OFFSET(bi);

    if (threadID < n) {
        temp[ai + bankOffsetA] = input[ai];
        temp[bi + bankOffsetB] = input[bi];
    } else {
        temp[ai + bankOffsetA] = 0;
        temp[bi + bankOffsetB] = 0;
    }

    int offset = 1;
    for (int d = powerOfTwo >> 1; d > 0; d >>= 1)  // build sum in place up the tree
    {
        __syncthreads();
        if (threadID < d) {
            int ai = offset * (2 * threadID + 1) - 1;
            int bi = offset * (2 * threadID + 2) - 1;
            ai += CONFLICT_FREE_OFFSET(ai);
            bi += CONFLICT_FREE_OFFSET(bi);

            temp[bi] += temp[ai];
        }
        offset *= 2;
    }

    if (threadID == 0) {
        temp[powerOfTwo - 1 + CONFLICT_FREE_OFFSET(powerOfTwo - 1)] = 0;  // clear the last element
    }

    for (int d = 1; d < powerOfTwo; d *= 2)  // traverse down tree & build scan
    {
        offset >>= 1;
        __syncthreads();
        if (threadID < d) {
            int ai = offset * (2 * threadID + 1) - 1;
            int bi = offset * (2 * threadID + 2) - 1;
            ai += CONFLICT_FREE_OFFSET(ai);
            bi += CONFLICT_FREE_OFFSET(bi);

            int t = temp[ai];
            temp[ai] = temp[bi];
            temp[bi] += t;
        }
    }
    __syncthreads();

    if (threadID < n) {
        output[ai] = temp[ai + bankOffsetA];
        output[bi] = temp[bi + bankOffsetB];
    }
}

__global__ void prescan_large(int *output, int *input, int n, int *sums) {
    extern __shared__ int temp[];

    int blockID = blockIdx.x;
    int threadID = threadIdx.x;
    int blockOffset = blockID * n;

    int ai = threadID;
    int bi = threadID + (n / 2);
    int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
    int bankOffsetB = CONFLICT_FREE_OFFSET(bi);
    temp[ai + bankOffsetA] = input[blockOffset + ai];
    temp[bi + bankOffsetB] = input[blockOffset + bi];

    int offset = 1;
    for (int d = n >> 1; d > 0; d >>= 1)  // build sum in place up the tree
    {
        __syncthreads();
        if (threadID < d) {
            int ai = offset * (2 * threadID + 1) - 1;
            int bi = offset * (2 * threadID + 2) - 1;
            ai += CONFLICT_FREE_OFFSET(ai);
            bi += CONFLICT_FREE_OFFSET(bi);

            temp[bi] += temp[ai];
        }
        offset *= 2;
    }
    __syncthreads();

    if (threadID == 0) {
        sums[blockID] = temp[n - 1 + CONFLICT_FREE_OFFSET(n - 1)];
        temp[n - 1 + CONFLICT_FREE_OFFSET(n - 1)] = 0;
    }

    for (int d = 1; d < n; d *= 2)  // traverse down tree & build scan
    {
        offset >>= 1;
        __syncthreads();
        if (threadID < d) {
            int ai = offset * (2 * threadID + 1) - 1;
            int bi = offset * (2 * threadID + 2) - 1;
            ai += CONFLICT_FREE_OFFSET(ai);
            bi += CONFLICT_FREE_OFFSET(bi);

            int t = temp[ai];
            temp[ai] = temp[bi];
            temp[bi] += t;
        }
    }
    __syncthreads();

    output[blockOffset + ai] = temp[ai + bankOffsetA];
    output[blockOffset + bi] = temp[bi + bankOffsetB];
}

__global__ void add(int *output, int length, int *n) {
    int blockID = blockIdx.x;
    int threadID = threadIdx.x;
    int blockOffset = blockID * length;

    output[blockOffset + threadID] += n[blockID];
}

__global__ void add(int *output, int length, int *n1, int *n2) {
    int blockID = blockIdx.x;
    int threadID = threadIdx.x;
    int blockOffset = blockID * length;

    output[blockOffset + threadID] += n1[blockID] + n2[blockID];
}

int nextPowerOfTwo(int x) {
    int power = 1;
    while (power < x) {
        power *= 2;
    }
    return power;
}