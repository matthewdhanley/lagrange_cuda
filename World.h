//
// Created by matt on 5/12/19.
//

#ifndef PETER_LAGRANGE_WORLD_H
#define PETER_LAGRANGE_WORLD_H


class World {
    public:
        World(int num_particles);

    protected: // data
        int num_particles;

        float resolution;

        float* x_axis;
        float* y_axis;
        float* z_axis;

};


#endif //PETER_LAGRANGE_WORLD_H
