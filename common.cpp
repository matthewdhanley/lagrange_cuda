#include "common.h"

bool load_UVW(float* U, float* V, float* W, int t) {
    FILE *fid;
    t = t-1;
    char *fbase = (char *) malloc(49 * sizeof(char));

    sprintf(fbase, "/media/matt/UNTITLED/Velocity%d/Velocity%d_%04d.bin", 1, 1, t);
//    sprintf(fbase, "/media/matt/Windows/Users/mhanl/Velocity%d_0000.bin", 1);

    printf("/media/matt/UNTITLED/Velocity%d/Velocity%d_%04d.bin\n", 1, 1, t);
//    printf("/media/matt/Windows/Users/mhanl/Velocity%d_%04d.bin\n", 1, 1, t);
    fid = fopen(fbase, "rb");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file 1\n");
        return 0;
    }
    fread(U, sizeof(float), NX * NY * NZ, fid);
    fclose(fid);
    sprintf(fbase, "/media/matt/UNTITLED/Velocity%d/Velocity%d_%04d.bin", 2, 2, t);
//    sprintf(fbase, "/media/matt/Windows/Users/mhanl/Velocity%d_0000.bin", 2);
    fid = fopen(fbase, "rb");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file 2\n");
        return 0;
    }
    fread(V, sizeof(float), NX * NY * NZ, fid);
    fclose(fid);
    sprintf(fbase, "/media/matt/UNTITLED/Velocity%d/Velocity%d_%04d.bin", 3, 3, t);
//    sprintf(fbase, "/media/matt/Windows/Users/mhanl/Velocity%d_0000.bin", 3);

    fid = fopen(fbase, "rb");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file 3\n");
        return 0;
    }
    fread(W, sizeof(float), NX * NY * NZ, fid);
    fclose(fid);

    // FOR TESTING WITHOUT HARD DRIVE.
//    sync();
//    int fd = open("/proc/sys/vm/drop_caches", O_WRONLY);
//    if(fd < 0)
//    {
//        printf("Error clearing filesystem cache. Try running as sudo.\n");
//    }

    return 1;
}



double time_elapsed(clock_t start){
    return ((float)clock()-start)/CLOCKS_PER_SEC;
}


void setup_points_cpu(Point3D* pt_lagrange){
    if (((NX * NY * NZ) % N_PARTICLES)){
        printf("Particles cannot be evenly divided into volume.\n");
        exit(1);
    }

//    int sparcity = NX * NY * NZ / N_PARTICLES;

    int pt_counter = 0;
    for (int i = 0; i < NX; i += SPARCITY){
        for (int j = 0; j < NY; j+= SPARCITY){
//            for (int k=0; k < NZ; k+= SPARCITY) {
//                pt_lagrange[pt_counter].x = .07;
//                pt_lagrange[pt_counter].y = .07;
//                pt_lagrange[pt_counter++].z = .63;
                pt_lagrange[pt_counter].x = RESOLUTION/2 + ((float) j) * RESOLUTION;
                pt_lagrange[pt_counter].y = RESOLUTION/2 + ((float) i) * RESOLUTION;
                pt_lagrange[pt_counter++].z = RESOLUTION * 1024;
//            }
        }
    }
    printf("pt_counter = %d\n",pt_counter);
}