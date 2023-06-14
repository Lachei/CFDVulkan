#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include "boundary_val.hpp"
#include "grid.hpp"
#include "par.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <iostream>
#include <mpi.h>

void calculate_temp(SimulationContext& sim_ctx, Grid& grid, matrix2d<double> &T){
	static std::vector<double> send, rec;
    matrix2d<double>& u = grid.velocity[0], &v = grid.velocity[1];
	MPI_Status mpi_s;
    T = grid.temperature;
    double min, max;

	matrix_comm(T, sim_ctx.rank_l, sim_ctx.rank_r, sim_ctx.rank_b, sim_ctx.rank_t, send, rec, &mpi_s, 0);
	
	for (int i = 2; i < T.width() - 1; ++i) {
		for (int j = 2; j < T.height() - 1; ++j) {
			double duT_dx = (1 / sim_ctx.dx) * (.5 * u(i,j) * (T(i,j) + T(i + 1,j)) - .5 * u(i - 1,j) * (T(i - 1,j) + T(i,j)))
				+ (sim_ctx.alpha / sim_ctx.dx) * (.5 * abs(u(i,j)) * (T(i,j) - T(i + 1,j)) - .5 * abs(u(i - 1,j)) * (T(i - 1,j) - T(i,j)));

			double dvT_dy = (1 / sim_ctx.dy) * (.5 * v(i,j) * (T(i,j) + T(i, j + 1)) - .5 * v(i,j - 1) * (T(i, j - 1) + T(i,j)))
				+ (sim_ctx.alpha / sim_ctx.dy) * (.5 * abs(v(i,j)) * (T(i,j) - T(i, j + 1)) - .5 * abs(v(i, j - 1)) * (T(i, j - 1) - T(i,j)));

			double d2T_dx2 = (T(i + 1,j) - 2 * T(i,j) + T(i - 1,j)) / (sim_ctx.dx * sim_ctx.dx);

			double d2T_dy2 = (T(i, j + 1) - 2 * T(i,j) + T(i, j - 1)) / (sim_ctx.dy * sim_ctx.dy);

			grid.temperature(i,j) = T(i,j) + sim_ctx.dt * ((1 / (sim_ctx.Re * sim_ctx.Pr)) * (d2T_dx2 + d2T_dy2) - duT_dx - dvT_dy);

            min = std::min(std::min(std::min(std::min(T(i, j), T(i - 1, j)), T(i + 1, j)), T(i, j - 1)), T(i, j + 1));
            max = std::max(std::max(std::max(std::max(T(i, j), T(i - 1, j)), T(i + 1, j)), T(i, j - 1)), T(i, j + 1));
            grid.temperature(i, j) = std::clamp(grid.temperature(i, j), min, max);
			//grid.temperature(i, j) = (std::max(min,grid.temperature(i, j)) == min) ? min : std::min(max, grid.temperature(i,j));
		}
	}
}

void calculate_matrix(SimulationContext& sim_ctx, Grid& grid, matrix2d<double>& M)
{
    static std::vector<double> send, rec;
    matrix2d<double>& u = grid.velocity[0], & v = grid.velocity[1];
    MPI_Status mpi_s;

    matrix_comm(M, sim_ctx.rank_l, sim_ctx.rank_r, sim_ctx.rank_b, sim_ctx.rank_t, send, rec, &mpi_s, 0);

    for (int i = 2; i < M.width() - 1; ++i) {
        for (int j = 2; j < M.height() - 1; ++j) {
            double duT_dx = (1 / sim_ctx.dx) * (.5 * u(i, j) * (M(i, j) + M(i + 1, j)) - .5 * u(i - 1, j) * (M(i - 1, j) + M(i, j)))
                + (sim_ctx.alpha / sim_ctx.dx) * (.5 * abs(u(i, j)) * (M(i, j) - M(i + 1, j)) - .5 * abs(u(i - 1, j)) * (M(i - 1, j) - M(i, j)));

            double dvT_dy = (1 / sim_ctx.dy) * (.5 * v(i, j) * (M(i, j) + M(i, j + 1)) - .5 * v(i, j - 1) * (M(i, j - 1) + M(i, j)))
                + (sim_ctx.alpha / sim_ctx.dy) * (.5 * abs(v(i, j)) * (M(i, j) - M(i, j + 1)) - .5 * abs(v(i, j - 1)) * (M(i, j - 1) - M(i, j)));

            double d2T_dx2 = (M(i + 1, j) - 2 * M(i, j) + M(i - 1, j)) / (sim_ctx.dx * sim_ctx.dx);

            double d2T_dy2 = (M(i, j + 1) - 2 * M(i, j) + M(i, j - 1)) / (sim_ctx.dy * sim_ctx.dy);

            M(i, j) = std::max(.0,M(i, j) + sim_ctx.dt * ((1 / (sim_ctx.Re * sim_ctx.Pr)) * (d2T_dx2 + d2T_dy2) - duT_dx - dvT_dy));
        }
    }
}

void normalize_species_ratios(SimulationContext& sim_ctx, Grid& grid)
{

    for (int i = 0; i < grid._imax; ++i) {
        for (int j = 0; j < grid._jmax; ++j) {
            double length = 0;
			if((int)grid.cell_flagg(i,j) > 9){
				for (int s = 0; s < grid.species_concentrations.size(); ++s) {
					if(grid.species_concentrations[s](i,j)>.0){
                		length += grid.species_concentrations[s](i, j);
					}
            	}
				length += 1e-20;
				for (int s = 0; s < grid.species_concentrations.size(); ++s) { 
					grid.species_concentrations[s](i, j) /= length;
				}
				convert_species(sim_ctx, grid, i, j, (int)grid.cell_flagg(i,j)-10); 
			}
			length=0;
            for (int s = 0; s < grid.species_concentrations.size(); ++s) {
				if(grid.species_concentrations[s](i,j)<0.0)std::cout<<grid.species_concentrations[s](i,j)<<s<<" "<< i<<" "<<j<<std::endl;
                length += grid.species_concentrations[s](i, j);
            }
            length += 1e-20;            // epsilon normalization if some error should occur
            grid.density(i, j) = 0;     // calculating density according to multispecies concentrations
            for (int s = 0; s < grid.species_concentrations.size(); ++s) {
                grid.species_concentrations[s](i, j) /= length;
                grid.density(i, j) += grid.species_concentrations[s](i, j) * grid.species[s].molar_mass * sim_ctx.PI / (GAS_CONSTANT * grid.temperature(i,j));
            }
        }
    }
}

void convert_species(SimulationContext& sim_ctx, Grid& grid, int i, int j, int c){
	std::vector<double> mass_in_cell(grid.species_concentrations.size());
	std::vector<double> reduction(grid.species_concentrations.size());
	

	for(int k=0;k<grid.species_concentrations.size();++k){
		//assert concentrations >0 (saftey check)
		if(grid.species_concentrations[k](i, j) <=.0){
			reduction[k]=.0;
			grid.species_concentrations[k](i, j)=.0;
			
		}else{
			//current mass of species k in cell
			mass_in_cell[k]=grid.species_concentrations[k](i, j) * grid.species[k].molar_mass * sim_ctx.PI / (GAS_CONSTANT * grid.temperature(i,j))*sim_ctx.dx*sim_ctx.dy;
			//reduction is either the defined reduction or all mass in cell
			if((sim_ctx.species_converters[c].reduction_per_qm_sec[k]>0)){
				if(TEST_CASE==1)reduction[k]=std::min(mass_in_cell[k], linear_c3_plant_consumption(sim_ctx.species_converters[c].reduction_per_qm_sec[k], grid.species_concentrations[k](i,j))*sim_ctx.dx*sim_ctx.dy*sim_ctx.dt);
				else if(TEST_CASE==2)reduction[k]=std::min(mass_in_cell[k], fire_O2_consumption(sim_ctx.species_converters[c].reduction_per_qm_sec[k],grid.species_concentrations[k](i,j))*sim_ctx.dx*sim_ctx.dy*sim_ctx.dt);
				else reduction[k]=std::min(mass_in_cell[k], sim_ctx.species_converters[c].reduction_per_qm_sec[k] *sim_ctx.dx*sim_ctx.dy*sim_ctx.dt);
			}else reduction[k]=0.0;
			
			if(reduction[k]>=mass_in_cell[k]){
				mass_in_cell[k]=.0;
			}else{
				mass_in_cell[k]-=reduction[k];
			}
		}
	}
	//convert consumed species into other species
	for(int k=0;k<grid.species_concentrations.size();++k){
		for(int z=0; z<grid.species_concentrations.size();++z){
			//species k gets defined ratio of reduction mass of species z
			mass_in_cell[k] += reduction[z] * sim_ctx.species_converters[c].conversion_ratio[k][z];
		}
		grid.species_concentrations[k](i, j)=mass_in_cell[k]/sim_ctx.dx/sim_ctx.dy/grid.species[k].molar_mass / sim_ctx.PI * (GAS_CONSTANT * grid.temperature(i,j));
	}
}

double linear_c3_plant_consumption(double base,double x){
	double y=(base/0.003)*x-(base/0.003*0.001); //linear consuption rate
	return std::clamp(y,.0,(base/0.003)*0.012-(base/0.003*0.001));
}

double fire_O2_consumption(double reduction, double concentration){
	return concentration > 0.16 ? reduction : 0;
}

// Determines the value of F and G
void calculate_fg(
        SimulationContext& sim_ctx,
        Grid& grid,
        matrix2d<double> &F,
        matrix2d<double> &G)
{
	matrix2d<double>& u = grid.velocity[0], &v = grid.velocity[1], &T = grid.temperature;

	//calculate f
	for (int i = 1; i < u.width() - 1; ++i) {
		for (int j = 1; j < u.height() - 1; ++j) {

			if((grid.cell_flagg(i,j)==cell_flag::FLUID || (int)grid.cell_flagg(i,j)>9)&& (grid.cell_flagg(i+1,j)==cell_flag::FLUID || grid.cell_flagg(i + 1, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i+1,j) > 9)){
				double d2u_dx2 = (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j )) / pow(sim_ctx.dx, 2);
				double d2u_dy2 = (u(i, j + 1) - 2 * u(i, j) + u(i , j  - 1)) / pow(sim_ctx.dy, 2);
				double du2_dx = (pow((u(i, j) + u(i + 1, j)) / 2, 2) - pow((u(i, j) + u(i - 1, j )) / 2, 2)) / sim_ctx.dx;
				double duv_dy = ((u(i, j) + u(i, j + 1)) * (v(i, j) + v(i + 1, j)) / 4 - (u(i , j  - 1) + u(i, j)) * (v(i , j  - 1) + v(i + 1 , j  - 1)) / 4) / sim_ctx.dy;
				double kW = (u(i, j) + u(i - 1, j )) / 2;
				double kE = (u(i, j) + u(i + 1, j)) / 2;
				double kS = (v(i, j) + v(i , j  - 1)) / 2;
				double kN = (v(i, j) + v(i, j + 1)) / 2;
				double dc_dx = sim_ctx.alpha / (2 * sim_ctx.dx) * (kE * (u(i, j) + u(i + 1, j)) - kW * (u(i - 1 , j ) + u(i, j)) + abs(kE) * (u(i, j)-u(i + 1, j)) - abs(kW) * (u(i - 1 , j ) - u(i, j))) + (1 - sim_ctx.alpha) * du2_dx;
				double dc_dy = sim_ctx.alpha / (2 * sim_ctx.dx) * (kN * (u(i, j) + u(i, j + 1)) - kS * (u(i , j  - 1) + u(i, j)) + abs(kN) * (u(i, j) - u(i, j + 1)) - abs(kS) * (u(i , j  - 1) - u(i, j))) + (1 - sim_ctx.alpha) * duv_dy;
				F(i, j) = u(i, j) + sim_ctx.dt * (1 / sim_ctx.Re * (d2u_dx2 + d2u_dy2) - dc_dx - dc_dy) + sim_ctx.dt * grid.density(i,j) * sim_ctx.GX * sim_ctx.dx * sim_ctx.dy;
			}
		}
	}
	//calculate g
    for (int i = 1; i < u.width() - 1; ++i) {
        for (int j = 1; j < u.height() - 1; ++j) {
			if((grid.cell_flagg(i,j)==cell_flag::FLUID || (int)grid.cell_flagg(i,j)>9) && (grid.cell_flagg(i, j + 1) == cell_flag::FLUID || grid.cell_flagg(i, j + 1) == cell_flag::OUTFLOW)|| (int)grid.cell_flagg(i,j+1) > 9){
				double d2v_dx2 = (v(i + 1, j) - 2 * v(i , j) + v(i - 1, j )) / pow(sim_ctx.dx, 2);
				double d2v_dy2 = (v(i, j + 1) - 2 * v(i , j) + v(i , j  - 1)) / pow(sim_ctx.dy, 2);
				double duv_dx = ((v(i , j) + v(i + 1, j)) * (u(i , j) + u(i, j + 1)) / 4 - (v(i - 1, j ) + v(i , j)) * (u(i - 1, j ) + u(i - 1,j + 1)) / 4) / sim_ctx.dx;
				double dv2_dy = (pow((v(i , j) + v(i, j + 1)) / 2, 2) - pow((v(i , j) + v(i , j  - 1)) / 2, 2)) / sim_ctx.dy;
				double kW = (u(i , j) + u(i - 1, j )) / 2;
				double kE = (u(i , j) + u(i + 1, j)) / 2;
				double kS = (v(i , j) + v(i , j  - 1)) / 2;
				double kN = (v(i , j) + v(i, j + 1)) / 2;
				double dc_dx = sim_ctx.alpha / (2 * sim_ctx.dx) * (kE * (v(i , j) + v(i + 1, j)) - kW * (v(i - 1, j ) + v(i , j)) + abs(kE) * (v(i , j) - v(i + 1, j)) - abs(kW) * (v(i - 1, j ) - v(i , j))) + (1 - sim_ctx.alpha) * duv_dx;
				double dc_dy = sim_ctx.alpha / (2 * sim_ctx.dx) * (kN * (v(i , j) + v(i, j + 1)) - kS * (v(i , j  - 1) + v(i , j)) + abs(kN) * (v(i , j) - v(i, j + 1)) - abs(kS) * (v(i , j  - 1) - v(i , j))) + (1 - sim_ctx.alpha) * dv2_dy;
				G(i , j) = v(i , j) + sim_ctx.dt * (1 / sim_ctx.Re * (d2v_dx2 + d2v_dy2) - dc_dy - dc_dx) + sim_ctx.dt * grid.density(i,j) * sim_ctx.GY * sim_ctx.dx * sim_ctx.dy;
			}
		}
	}
	//printf("sds: %f\n",G[1][5]);
	for(int i=1;i<F.width();i++){
		for(int j=1;j<F.height();j++){
			if(grid.cell_flagg(i,j)!=cell_flag::FLUID && (int)grid.cell_flagg(i,j) <= 9){
				set_fg_boundary(grid,F,G,i,j);
			}
		}
	}
	//printf("sds: %f\n",G[0][5]);
}

// This operation computes the right hand side of the pressure poisson equation.
void calculate_rs(
        SimulationContext& sim_ctx,
        matrix2d<double> &F,
        matrix2d<double> &G,
        matrix2d<double> &RS)
{
	for (int i = 1; i < F.width(); ++i) {
		for (int j = 1; j < F.height(); ++j) {
			RS(i,j) = 1 / sim_ctx.dt * ((F(i,j) - F(i - 1, j )) / sim_ctx.dx + (G(i,j) - G(i,j - 1)) / sim_ctx.dy);
		}
	}
}

// Determines the maximal time step size
void calculate_dt(SimulationContext& sim_ctx, Grid &grid)
{
        //Determine maximum absolute velocities
		double a=1/sim_ctx.Re/sim_ctx.Pr;
        double vmax=0;
        double umax=0;
        sim_ctx.alpha = 0;
        for(int i=0;i<grid._imax;i++){
                for(int j=0; j<grid._jmax;j++){
                        double utmp=abs(grid.velocity[int(velocity_type::U)](i,j));
                        double vtmp=abs(grid.velocity[int(velocity_type::V)](i,j));
                        if(vmax<vtmp)vmax=vtmp;
                        if(umax<utmp)umax=utmp;
                        double m = std::max(abs(grid.velocity[0](i, j) * sim_ctx.dt / sim_ctx.dx), abs(grid.velocity[1](i, j) * sim_ctx.dt / sim_ctx.dy));
                        if (m > sim_ctx.alpha) sim_ctx.alpha = m * 1.00000000000000000001;
                }
        }
		
        //satisfy all four conditions
        sim_ctx.dt=    sim_ctx.tau *                                                   //safty factor
                fmin(fmin(fmin(                                              //Minimum of the three conditions
                          
						.5*(sim_ctx.Pr* sim_ctx.Re* sim_ctx.dx* sim_ctx.dx* sim_ctx.dy* sim_ctx.dy)/(sim_ctx.dx* sim_ctx.dx+ sim_ctx.dy* sim_ctx.dy),			//Condition: (1/(2a)) *(1/dx^2 +1/dy^2)^-1, a = 1/Re/Pr
						.5*(sim_ctx.dx* sim_ctx.dx* sim_ctx.dy* sim_ctx.dy* sim_ctx.Re)/(sim_ctx.dx* sim_ctx.dx+ sim_ctx.dy* sim_ctx.dy)),				//Condition: (1/(2v)) *(1/dx^2 +1/dy^2)^-1
                        sim_ctx.dx/umax),                                       //CFL-Conditon
                        sim_ctx.dy/vmax                                         //CFL-Condtion
                ); 
			
		double tmp = sim_ctx.dt;
		if (MPI_Allreduce(&tmp, &sim_ctx.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS) {
			std::cout << "Error in reduction min for dt value." << std::endl;
			exit(-1);
		}
        tmp = sim_ctx.alpha;
        if (MPI_Allreduce(&tmp, &sim_ctx.alpha, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
            std::cout << "Error in reduction max for alpha value." << std::endl;
            exit(-1);
        }
}

void calculate_uv(
		SimulationContext& sim_ctx,
        Grid& grid,
        matrix2d<double> &F,
        matrix2d<double> &G)
{
	double d_old;
	matrix2d<double>& u = grid.velocity[0],& v = grid.velocity[1];

	for (int i = 1; i < u.width() - 1; i++) {
		for (int j = 1; j < u.height() - 1; j++) {
			if(grid.cell_flagg(i,j) == cell_flag::FLUID || grid.cell_flagg(i,j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i,j)>9){
				if((grid.cell_flagg(i+1,j)==cell_flag::FLUID || grid.cell_flagg(i + 1, j) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i+1,j) > 9)){
					u(i,j) = F(i,j) - sim_ctx.dt / sim_ctx.dx * (grid.pressure(i + 1,j) - grid.pressure(i, j));
				}
				if((grid.cell_flagg(i,j+1)==cell_flag::FLUID || grid.cell_flagg(i, j + 1) == cell_flag::OUTFLOW || (int)grid.cell_flagg(i,j+1) > 9)){
					v(i,j) = G(i,j) - sim_ctx.dt / sim_ctx.dy * (grid.pressure(i, j + 1) - grid.pressure(i, j));
				}
			}
		}
	}
	
	static std::vector<double> send, rec;
	MPI_Status mpi_s;
	uv_comm(u, v, sim_ctx.rank_l, sim_ctx.rank_r, sim_ctx.rank_b, sim_ctx.rank_t, send, rec, &mpi_s, 0);
}
