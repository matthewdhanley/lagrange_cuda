//
// Created by matt on 5/12/19.
//

#ifndef FINAL_LAGRANGE_CPU_H
#define FINAL_LAGRANGE_CPU_H
#include "common.h"
#include "globals.h"
#include "lagrange_cuda.cuh"

void cpu_streamlines(Point3D* pt_sav, float* x_axis, float* y_axis, float* z_axis, int t, float* U, float* V, float* W, float* Q_u, float* Q_v, float* Q_w);
#endif //FINAL_LAGRANGE_CPU_H
