#include "Points.h"

void Points::allocate(int size) {
    this->size = size;
    this->x = (float*)malloc(size * sizeof(float));
    this->y = (float*)malloc(size * sizeof(float));
    this->ptsCnt = (int*)malloc(size * sizeof(int));
    this->cluster = (int*)malloc(size * sizeof(int));
}
int Points::getSize() {
    return size;
}