#include <stdio.h>
#include <stdlib.h>

#ifndef SCAN_H_
#define SCAN_H_

#define THREADS_PER_BLOCK 512
#define ELEMENTS_PER_BLOCK (THREADS_PER_BLOCK * 2)

#define SHARED_MEMORY_BANKS 32
#define LOG_MEM_BANKS 5
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_MEM_BANKS)

void scan(int *d_out, int *d_in, int length);
void scanLargeDeviceArray(int *d_out, int *d_in, int length);
void scanSmallDeviceArray(int *d_out, int *d_in, int length);
void scanLargeEvenDeviceArray(int *d_out, int *d_in, int length);

__global__ void prescan_arbitrary(int *output, int *input, int n, int powerOfTwo);
__global__ void prescan_large(int *output, int *input, int n, int *sums);
__global__ void add(int *output, int length, int *n);
__global__ void add(int *output, int length, int *n1, int *n2);

int nextPowerOfTwo(int x);

#endif