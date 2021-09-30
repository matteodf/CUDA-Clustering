#ifndef POINTS_H
#define POINTS_H

#include <stdlib.h>

#include <vector>

class Points {
   public:
    float* x;
    float* y;
    int* ptsCnt;
    int* cluster;
    int size;
    void allocate(int size);
    int getSize();
};

#endif