#ifndef FINAL_LAGRANGE_COMMON_H
#define FINAL_LAGRANGE_COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "globals.h"

typedef struct Points3D {
    float x;
    float y;
    float z;
} Point3D;

typedef struct Points3D_int {
    int x;
    int y;
    int z;
} Point3D_int;

typedef struct Points2D {
    float x;
    float y;
} Point2D;

bool load_UVW(float* U, float* V, float* W, int t);
double time_elapsed(clock_t start);
void setup_points_cpu(Point3D* pt_lagrange);

#endif //FINAL_LAGRANGE_COMMON_H
