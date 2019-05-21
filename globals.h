#ifndef FINAL_LAGRANGE_GLOBALS_H
#define FINAL_LAGRANGE_GLOBALS_H

#include "common.h"

#define GPU 1
#define CPU 0

#define N_PARTICLES 4096
#define SPARCITY 4
//#define N_PARTICLES 8192
//#define SPARCITY 128
//#define N_PARTICLES 32768
//#define SPARCITY 64
//#define N_PARTICLES 131072
//#define SPARCITY 32
//#define N_PARTICLES 524288
//#define SPARCITY 16
//#define N_PARTICLES 2097152
//#define SPARCITY 8
//#define N_PARTICLES 8388608
//#define SPARCITY 4
//#define N_PARTICLES 33554432
//#define SPARCITY 2
//#define N_PARTICLES 134217728
//#define SPARCITY 1
#define NX 256
#define NY 256
#define NZ 2048
#define RESOLUTION 0.0027358225
#define TSTEPS 1280  // this is defined by the number of files of gradients given.
//#define TSTEPS 3  // this is defined by the number of files of gradients given.
#define ZERO_THRESH 1e-10
#define DT 0.00000117402
#define MINVAL 0.0f + RESOLUTION
#define MAXVALX (float)NX*RESOLUTION - RESOLUTION
#define MAXVALY (float)NY*RESOLUTION - RESOLUTION
#define MAXVALZ (float)NZ*RESOLUTION - RESOLUTION
//#define THREADS_PER_BLOCK 512
//#define NBLOCKS 15000
#define THREADS_PER_BLOCK 512
#define NBLOCKS 150

#endif //FINAL_LAGRANGE_GLOBALS_H
