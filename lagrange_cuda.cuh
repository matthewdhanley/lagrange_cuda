#ifndef LAGRANGE_CUDA_CUH_
#define LAGRANGE_CUDA_CUH_

#include <cuda.h>
#include "globals.h"
#include "common.h"

__host__ __device__ void get_Q(float* Q, float* M, Point3D_int low, Point3D_int high);
__host__ __device__ Point3D get_next_rk4_k(Point3D pt, float* x_axis, float* y_axis, float* z_axis,
                                           float* U, float* V, float* W,  float* Q_u, float* Q_v, float* Q_w);
__host__ __device__ float trilinear(float* Q, Point3D pt1, Point3D pt2, Point3D pt);
__host__ __device__ float bilinear(float* Q, Point2D pt1, Point2D pt2, Point2D pt);
__host__ __device__ Point3D lagrange3_step(Point3D pt, float* x_axis, float* y_axis, float* z_axis, float* U, float* V, float* W, float* Q_u, float* Q_v, float* Q_w);
__global__ void streamlines(int t, float* x_axis, float* y_axis, float* z_axis,
                            float* U, float* V, float* W,
                            float* Q_u, float* Q_v, float* Q_w,
                            Point3D* pt_sav);
__global__ void setup_axes(float* x_axis, float* y_axis, float* z_axis);
__global__ void setup_pt_lagrange(Point3D* pt_lagrange, int npoints);
//__global__ void setup_sav(Point3D* pt_sav, Point3D* pt_lagrange);


#endif // LAGRANGE_CUDA_CUH