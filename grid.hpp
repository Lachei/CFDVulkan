#ifndef CFDLAB_GRID_H
#define CFDLAB_GRID_H
#include <vector>
#include "datastructures.hpp"
#include "vk_utils.h"
#include "gpu_images.h"

/*
 * CamelCase for function that returns or sets Objects
 * snake_case for functions that return a specific attribute like imax or jmax
 */

class Grid {
public:
    Grid(SimulationContext& sim_ctx, matrix2d<cell_flag>& cell_flags);

    matrix2d<double> velocity[2];
    matrix2d<double> pressure;
    matrix2d<double> temperature;
    matrix2d<double> density;
    matrix2d<cell_flag> cell_flagg;
    matrix2d<Neighbourhood> neighbourhood;
    std::vector<matrix2d<double>> species_concentrations;
    std::vector<Species> species;

    int _imax, _jmax;

    int i_co2 = -1;
    int i_o2 = -1;
    // Print matrices
    void print_velocity(velocity_type type);
    void print_pressure();
    void print_cell_flags();
    void print_neighbourhood();
private:
};


class GpuGrid {
public:
	GpuGrid(SimulationContext& sim_ctx, VkUtil::VkContext& vk_context, matrix2dRM<cell_flag>& flags);
    
    int _imax, _jmax;
    
    GpuImage velocity;
    GpuImage pressure;
    GpuDoubleImage temperature;
    GpuImage density;
    GpuImage cell_flags;
    GpuImage neighbourhood;
    GpuDoubleArrayImage species_concentrations;
    GpuImage fg;
    GpuImage rs;

    matrix2dRM<cell_flag> cell_flagg;
    matrix2dRM<uint8_t> neighbours;

    int i_co2 = -1;
    int i_o2 = -1;
private:
};

#endif //CFDLAB_GRID_H
