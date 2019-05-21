#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "common.h"
#include "globals.h"
#include "cpu.h"
//#include "World.h"
#include "lagrange_cuda.cuh"
#include "gpu_err_check.h"

int main(){
    cudaFree(0);
    /** Set up Axes **/
    // initialize pointers to the axes' arrays
    float* x_axis; float* y_axis; float* z_axis;

    // allocate space for axes
    x_axis = (float*) malloc((size_t)(int)((NX) * sizeof(float)));
    y_axis = (float*) malloc((size_t)(int)((NY) * sizeof(float)));
    z_axis = (float*) malloc((size_t)(int)((NZ) * sizeof(float)));

    // assign values to x, y, and z axes
    for (int i = 0; i < NX; i++){
        x_axis[i] = ((float) i) * RESOLUTION;
    }

    for (int i = 0; i < NY; i++){
        y_axis[i] = ((float) i) * RESOLUTION;
    }

    for (int i = 0; i < NZ; i++){
        z_axis[i] = ((float) i) * RESOLUTION;
    }

    /** Initialize points that will be tracked **/

    Point3D* pt_sav = (Point3D*) malloc(N_PARTICLES * sizeof(Point3D));
    Point3D* pt_sav_out = (Point3D*) malloc(N_PARTICLES * sizeof(Point3D));

    /** Setup Gradient Variables **/
    float* U = (float*) malloc((size_t) NX*NY*NZ * sizeof(float));
    float* V = (float*) malloc((size_t) NX*NY*NZ * sizeof(float));
    float* W = (float*) malloc((size_t) NX*NY*NZ * sizeof(float));

    float* Q_u = (float*) malloc((size_t) 8 * sizeof(float));
    float* Q_v = (float*) malloc((size_t) 8 * sizeof(float));
    float* Q_w = (float*) malloc((size_t) 8 * sizeof(float));

    /** GPU VARIABLES **/
    float* d_x_axis;
    float* d_y_axis;
    float* d_z_axis;
    cudaSafeCall( cudaMalloc((void**) &d_x_axis, (int)(NX * sizeof(float)) ));
    cudaSafeCall( cudaMalloc((void**) &d_y_axis, (int)(NY * sizeof(float)) ));
    cudaSafeCall( cudaMalloc((void**) &d_z_axis, (int)(NZ * sizeof(float)) ));
    setup_axes<<<8, 32>>>(d_x_axis, d_y_axis, d_z_axis);
    cudaDeviceSynchronize();

    float* d_U;
    float* d_V;
    float* d_W;
    cudaSafeCall( cudaMalloc((void**) &d_U, NX*NY*NZ*sizeof(float)) );
    cudaSafeCall( cudaMalloc((void**) &d_V, NX*NY*NZ*sizeof(float)) );
    cudaSafeCall( cudaMalloc((void**) &d_W, NX*NY*NZ*sizeof(float)) );

    float *d_Q_u;
    float *d_Q_v;
    float *d_Q_w;
    cudaSafeCall( cudaMalloc((void **) &d_Q_u, NBLOCKS * THREADS_PER_BLOCK * 8 * sizeof(float)) );
    cudaSafeCall( cudaMalloc((void **) &d_Q_v, NBLOCKS * THREADS_PER_BLOCK * 8 * sizeof(float)) );
    cudaSafeCall( cudaMalloc((void **) &d_Q_w, NBLOCKS * THREADS_PER_BLOCK * 8 * sizeof(float)) );

//    dim3 threads_per_block(8,8,8);
//    dim3 blocks(32,32,256);
//    setup_pt_lagrange <<< blocks, threads_per_block >>> (d_sav, N_PARTICLES);
//    cudaDeviceSynchronize();
//    cudaCheckError();


    Point3D *d_pt_sav;
    cudaSafeCall(cudaMalloc((void **) &d_pt_sav, (N_PARTICLES) * sizeof(Point3D)));
    cudaDeviceSynchronize();

    setup_points_cpu(pt_sav);
    cudaSafeCall( cudaMemcpy(d_pt_sav, pt_sav, N_PARTICLES * sizeof(Point3D), cudaMemcpyHostToDevice) );
    cudaDeviceSynchronize();

    /** LOOP: for each timestep of the gradient data **/
    clock_t start;
    clock_t cpu_time;
    clock_t gpu_time;
    clock_t write_time;

    double gpu_time_elapsed;
    double cpu_time_elapsed;
    double load_time_elapsed;
    double write_time_elapsed;

    int t = 1;

    if (load_UVW(U, V, W, t) == 0) {
        printf("Error with loading.\n");
        exit(-1);
    }

    FILE* outfile;
    outfile = fopen("/media/matt/UNTITLED/output_val.bin", "wb");

    for (int i = 0; i<=TSTEPS; i++){
        /** Load gradients from file **/

        start = clock();

//        // load memory
//        if ((t-1) % 5 == 0) {
//            if (load_UVW(U, V, W, t) == 0) {
//                printf("Error with loading.\n");
//                fclose(outfile);
//                exit(-1);
//            }
//        }

        load_time_elapsed = time_elapsed(start);

        /** ON THE CPU **/
        cpu_time = clock();
        if (CPU and ~GPU){
            if ((t-1) % 5 == 0) {
                if (load_UVW(U, V, W, t) == 0) {
                    printf("Error with loading.\n");
                    fclose(outfile);
                    exit(-1);
                }
            }
        }
        if (CPU) cpu_streamlines(pt_sav, x_axis, y_axis, z_axis, t, U, V, W, Q_u, Q_v, Q_w);
        cpu_time_elapsed = time_elapsed(cpu_time);

        /** ON THE GPU **/
        if (GPU) {
            gpu_time = clock();

            // copy gradients to the device
            if ((t-1) % 5 == 0) {
                cudaSafeCall(cudaMemcpy(d_U, U, NX * NY * NZ * sizeof(float), cudaMemcpyHostToDevice));
                cudaSafeCall(cudaMemcpy(d_V, V, NX * NY * NZ * sizeof(float), cudaMemcpyHostToDevice));
                cudaSafeCall(cudaMemcpy(d_W, W, NX * NY * NZ * sizeof(float), cudaMemcpyHostToDevice));
            }

            // run GPU Analysis
            streamlines << < NBLOCKS, THREADS_PER_BLOCK >> >
                                      (t, d_x_axis, d_y_axis, d_z_axis, d_U, d_V, d_W, d_Q_u, d_Q_v, d_Q_w, d_pt_sav);
            if ((t) % 5 == 0) {
                if (load_UVW(U, V, W, t+1) == 0) {
                    printf("Error with loading.\n");
                    fclose(outfile);
                    exit(-1);
                }
            }
            cudaDeviceSynchronize();
            cudaCheckError();

            // copy memory back
            cudaSafeCall(cudaMemcpy(pt_sav_out, d_pt_sav, N_PARTICLES * sizeof(Point3D), cudaMemcpyDeviceToHost));
            cudaDeviceSynchronize();

            printf("X: %1.16f\t", pt_sav[1000].x);
            printf("Y: %1.16f\t", pt_sav[1000].y);
            printf("Z: %1.16f\n", pt_sav[1000].z);
        }
        gpu_time_elapsed = time_elapsed(gpu_time);
        write_time = clock();

        fwrite(pt_sav, sizeof(Point3D), N_PARTICLES, outfile);
        write_time_elapsed = time_elapsed(write_time);
        t++;

        // print timing metrics.
        printf("Iteration: %3d | Load: %f | Write: %f | GPU: %f | CPU: %f\n", i, load_time_elapsed, write_time_elapsed, gpu_time_elapsed, cpu_time_elapsed);
        time_t rawtime;
        struct tm * timeinfo;

        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        printf ( "Current local time and date: %s", asctime (timeinfo) );
    }
    fclose(outfile);


    return 0;
}