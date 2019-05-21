#include "lagrange_cuda.cuh"


__host__ __device__ float bilinear(float* Q, Point2D pt1, Point2D pt2, Point2D pt){

    if (fabsf(pt1.x - pt.x) < ZERO_THRESH || fabsf(pt2.x - pt.x) < ZERO_THRESH){
        return (pt2.y + pt1.y)/2.f;
    }
    else if (fabsf(pt1.y - pt.y) < ZERO_THRESH || fabsf(pt2.y - pt.y) < ZERO_THRESH){
        return (pt2.x + pt1.x)/2.f;
    }
    else {

        float first_term = 1 / ((pt2.x - pt1.x) * (pt2.y - pt1.y));
        return first_term * (Q[0] * (pt1.x - pt.x) * (pt2.y - pt.y) +
                             Q[1] * (pt.x - pt1.x) * (pt2.y - pt.y) +
                             Q[2] * (pt2.x - pt.x) * (pt.y - pt1.y) +
                             Q[3] * (pt.x - pt1.x) * (pt.y - pt1.y));
    }
}


__host__ __device__ float trilinear(float* Q, Point3D pt1, Point3D pt2, Point3D pt){
    /**
     * Q - 3d array
    **/
    Point2D bpt1 = {pt1.x, pt1.y};
    Point2D bpt2 = {pt2.x, pt2.y};
    Point2D bpt = {pt.x,  pt.y};
    float* Q1 = &Q[0];
    float* Q2 = &Q[4];
    if (pt1.z - pt.z == 0){
        return bilinear(Q1, bpt1, bpt2, bpt);
    }
    else if (pt2.z - pt.z == 0){
        return bilinear(Q2, bpt1, bpt2, bpt);
    }
    else if (pt2.z == pt1.z){
        return bilinear(Q1, bpt1, bpt2, bpt);
    }
    else {
        float zd = (pt.z - pt1.z) / (pt2.z - pt1.z);

        float b1 = bilinear(Q1, bpt1, bpt2, bpt);
        float b2 = bilinear(Q2, bpt1, bpt2, bpt);

        return b1 * (1.f - zd) + b2 * zd;
    }
}


__host__ __device__ void get_Q(float* Q, float* M, Point3D_int low, Point3D_int high){
//    printf("%d ",low.z  * (NX * NY) + low.y  * NY  + low.x );
    Q[0] = M[low.z  * (NX * NY) + low.y  * NY  + low.x ];
    Q[1] = M[low.z  * (NX * NY) + low.y  * NY  + high.x];
    Q[2] = M[low.z  * (NX * NY) + high.y * NY  + low.x ];
    Q[3] = M[low.z  * (NX * NY) + high.y * NY  + high.x];
    Q[4] = M[high.z * (NX * NY) + low.y  * NY  + low.x ];
    Q[5] = M[high.z * (NX * NY) + low.y  * NY  + high.x];
    Q[6] = M[high.z * (NX * NY) + low.y  * high.y + low.x ];
    Q[7] = M[high.z * (NX * NY) + low.y  * high.y + high.x];
}


__host__ __device__ Point3D get_next_rk4_k(Point3D pt, float* x_axis, float* y_axis, float* z_axis,
                                           float* U, float* V, float* W, float* Q_u, float* Q_v, float* Q_w){

    Point3D k = {
            pt.x / (float) RESOLUTION,
            pt.y / (float) RESOLUTION,
            pt.z / (float) RESOLUTION
    };

    Point3D_int k_high = {
            (int) ceilf(k.x),
            (int) ceilf(k.y),
            (int) ceilf(k.z)
    };

    Point3D_int k_low = {
            (int) floorf(k.x),
            (int) floorf(k.y),
            (int) floorf(k.z)
    };
//    printf("%f %f %f\n", pt.x, pt.y, pt.z);
//    printf("%d %d %d\n", k_high.x, k_high.y, k_high.z);
    get_Q(Q_u, U, k_low, k_high);
    get_Q(Q_v, V, k_low, k_high);
    get_Q(Q_w, W, k_low, k_high);

    Point3D low = {
            x_axis[k_low.x],
            y_axis[k_low.y],
            z_axis[k_low.z]
    };

    Point3D high = {
            x_axis[k_high.x],
            y_axis[k_high.y],
            z_axis[k_high.z]
    };


    Point3D d_dt = {
            trilinear(Q_u, low, high, pt),
            trilinear(Q_v, low, high, pt),
            trilinear(Q_w, low, high, pt)
    };


    Point3D rk4_k = {
            (float) DT * d_dt.x,
            (float) DT * d_dt.y,
            (float) DT * d_dt.z
    };
    return rk4_k;
}


__host__ __device__ Point3D lagrange3_step(Point3D pt, float* x_axis, float* y_axis, float* z_axis, float* U, float* V,
                                           float* W, float* Q_u, float* Q_v, float* Q_w){
    Point3D nan_point = {NAN, NAN, NAN};

    Point3D k1 = get_next_rk4_k(pt, x_axis, y_axis, z_axis, U, V, W, Q_u, Q_v, Q_w);
    Point3D tmp = {
            pt.x + k1.x/2.f,
            pt.y + k1.y/2.f,
            pt.z + k1.z/2.f
    };
    if (tmp.x < (float) MINVAL || tmp.y < (float) MINVAL || tmp.z < (float) MINVAL ||
        tmp.x > (float) (MAXVALX-RESOLUTION) || tmp.y > (float) (MAXVALY-RESOLUTION) || tmp.z > (float) (MAXVALZ-RESOLUTION)
        || tmp.x != tmp.x || tmp.y != tmp.y || tmp.z != tmp.z ){
        return nan_point;
    }
//    printf("from step %f %f %f\n", tmp.x, tmp.y, tmp.z);
    Point3D k2 = get_next_rk4_k(tmp, x_axis, y_axis, z_axis, U, V, W, Q_u, Q_v, Q_w);
    tmp = (Point3D){
            pt.x + k2.x/2.f,
            pt.y + k2.y/2.f,
            pt.z + k2.z/2.f
    };
    if (tmp.x < (float) MINVAL || tmp.y < (float) MINVAL || tmp.z < (float) MINVAL ||
        tmp.x > (float) (MAXVALX-RESOLUTION) || tmp.y > (float) (MAXVALY-RESOLUTION) || tmp.z > (float) (MAXVALZ-RESOLUTION)
        || tmp.x != tmp.x || tmp.y != tmp.y || tmp.z != tmp.z ){
        return nan_point;
    }
    Point3D k3 = get_next_rk4_k(tmp, x_axis, y_axis, z_axis, U, V, W, Q_u, Q_v, Q_w);
    tmp = (Point3D){
            pt.x + k3.x/2.f,
            pt.y + k3.y/2.f,
            pt.z + k3.z/2.f
    };

    if (tmp.x < (float) MINVAL || tmp.y < (float) MINVAL || tmp.z < (float) MINVAL ||
        tmp.x > (float) (MAXVALX-RESOLUTION) || tmp.y > (float) (MAXVALY-RESOLUTION) || tmp.z > (float) (MAXVALZ-RESOLUTION)
        || tmp.x != tmp.x || tmp.y != tmp.y || tmp.z != tmp.z ){
        return nan_point;
    }
    Point3D k4 = get_next_rk4_k(tmp, x_axis, y_axis, z_axis, U, V, W, Q_u, Q_v, Q_w);
    Point3D result = {
            pt.x + 1/6.f * (k1.x + 2*k2.x + 2*k3.x + k4.x),
            pt.y + 1/6.f * (k1.y + 2*k2.y + 2*k3.y + k4.y),
            pt.z + 1/6.f * (k1.z + 2*k2.z + 2*k3.z + k4.z)
    };
    return result;
}


__global__ void streamlines(int t, float* x_axis, float* y_axis, float* z_axis,
                            float* U, float* V, float* W,
                            float* Q_u, float* Q_v, float* Q_w,
                            Point3D* pt_sav){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = tid;
    while (stride < N_PARTICLES){
        if (pt_sav[stride].x < 0){
            stride += gridDim.x * blockDim.x;
            continue;
        }

        Point3D pt_i = pt_sav[stride];
        Point3D lagrange = lagrange3_step(pt_i, x_axis, y_axis, z_axis, U, V, W, &Q_u[8*tid], &Q_v[8*tid], &Q_w[8*tid]);

        if (lagrange.x >= MAXVALX || lagrange.y >= MAXVALY || lagrange.z >= MAXVALZ ||
            lagrange.x <= MINVAL || lagrange.y <= MINVAL || lagrange.z <= MINVAL ||
            lagrange.x == NAN || lagrange.y == NAN || lagrange.z == NAN){
            pt_sav[stride].x = -1;
            stride += gridDim.x * blockDim.x;
            continue;
        }

        pt_sav[stride] = lagrange;
        stride += gridDim.x * blockDim.x;
    }

}

__global__ void setup_axes(float* x_axis, float* y_axis, float* z_axis){
    // assign x axis values
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = tid;

    while (stride < (int)(NX)){
        float val = ((float) stride) * RESOLUTION;
        x_axis[stride] = val;
        y_axis[stride] = val;
        stride += gridDim.x * blockDim.x;
    }
    stride = tid;
    while (stride < (int)(NZ)){
        float val = ((float) stride) * RESOLUTION;
        z_axis[stride] = val;
        stride += gridDim.x * blockDim.x;
    }
}


__global__ void setup_pt_lagrange(Point3D* pt_lagrange, int npoints){
    // initial points.
    int tidx = threadIdx.x + blockIdx.x * blockDim.x;
    int tidy = threadIdx.y + blockIdx.y * blockDim.y;
    int tidz = threadIdx.z + blockIdx.z * blockDim.z;
    int tid = tidz * blockDim.x * gridDim.x * blockDim.y * gridDim.y + tidy * blockDim.y * gridDim.y + tidx;
    int stride = tid;

    while (stride < npoints / 8){
        pt_lagrange[stride].x = (float) tidx * RESOLUTION * 8;
        pt_lagrange[stride].y = (float) tidy * RESOLUTION * 8;
        pt_lagrange[stride].z = (float) tidz * RESOLUTION * 8;
        stride += blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
    }
}


//__global__ void setup_sav(Point3D* pt_sav, Point3D* pt_lagrange) {
//    int tid = threadIdx.x + blockIdx.x * blockDim.x;
//    int stride = tid;
//    while (stride < N_PARTICLES) {
//        pt_sav[stride] = pt_lagrange[stride];
//        stride += gridDim.x * blockDim.x;
//
//    }
//}


//
//
//void streamlines_driver(int t, float* x_axis, float* y_axis, float* z_axis,
//                        float* U, float* V, float* W,
//                        float* Q_u, float* Q_v, float* Q_w, Point3D* pt_sav){
//    Point3D nan_point = {NAN, NAN, NAN};
//
//}
//

