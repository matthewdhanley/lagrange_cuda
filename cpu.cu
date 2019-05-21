#include "cpu.h"

void cpu_streamlines(Point3D* pt_sav, float* x_axis, float* y_axis, float* z_axis, int t, float* U, float* V, float* W, float* Q_u, float* Q_v, float* Q_w){

    for (int j=0; j<N_PARTICLES; j++){

        // skip if the point is no longer relevant
        if (pt_sav[j].x < 0){
            continue;
        }
        // load up point at previous time step
        Point3D pt = pt_sav[j];

        Point3D lagrange = lagrange3_step(pt, x_axis, y_axis, z_axis, U, V, W, Q_u, Q_v, Q_w);

        if (lagrange.x >= MAXVALX || lagrange.y >= MAXVALY || lagrange.z >= MAXVALZ ||
            lagrange.x <= MINVAL || lagrange.y <= MINVAL || lagrange.z <= MINVAL ||
            lagrange.x != lagrange.x || lagrange.y != lagrange.y || lagrange.z != lagrange.z){
            pt_sav[j].x = -1;
        }

        else{
            // save the new point
            pt_sav[j] = lagrange;
        }
    }
}
