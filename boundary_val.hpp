#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"
#include "grid.hpp"
#include "cell.hpp"

#define MAXFLOAT FLT_MAX

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(SimulationContext& sim_ctx, Grid& grid);
void set_fg_boundary(Grid& g, matrix2d<double> &F, matrix2d<double> &G, int i, int j);
void set_pressure_boundary(Neighbourhood n,matrix2d<double> &P, int i, int j);
#endif
