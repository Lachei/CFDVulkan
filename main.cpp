#include "helper.hpp"
#include "visual.hpp"
#include "init.hpp"
#include "gpu_pipelines.h"
#include "par.hpp"
#include "uvp.hpp"
#include "boundary_val.hpp"
#include "sor.hpp"
#include <cstdio>
#include <iostream>
#include <cassert>
#include <chrono> 
#include <iomanip>
#include <mpi.h>
using namespace std::chrono;

void print_matrix(matrix<double>& m) {
	for (auto& x : m) {
		for (double y : x) {
			std::cout << y << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void print_matrix(matrix2dRM<double>& m) {
    for (int x = 0; x < m.width(); ++x) {
        for (int y = 0; y < m.height(); ++y) {
            std::cout << std::setw(5) << std::setprecision(2) << m(x,y) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_matrix(matrix2dRM<cell_flag>& m) {
	for (int y = 0; y < m.height(); ++y) {
		for (int x = 0; x < m.width(); ++x) {
			std::cout << std::setprecision(2) << int(m(x, y)) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void print_matrix(matrix2dRM<float>& m) {
	for (int y = 0; y < m.height(); ++y) {
		for (int x = 0; x < m.width(); ++x) {
			std::cout << std::setw(5) << std::setprecision(2) << m(x, y) << " ";
		}
		std::cout << ";" << std::endl;
	}
	std::cout << std::endl;
}

void print_matrix(matrix2dRM<uint8_t>& m) {
	for (int y = 0; y < m.height(); ++y) {
		for (int x = 0; x < m.width(); ++x) {
			std::cout  << std::setprecision(2) << int(m(x, y)) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed. Use the predefined matrix<typename> type and give initial values in the constructor.
 * - perform the main loop
 * - at the end: destroy any memory allocated and print some useful statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two-dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop, the following steps are required (for some of the 
 * operations, a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

int main(int argn, char** args){
	// Initialize MPI
	MPI_Init(&argn, &args);

	std::string parameterFilename = "simulation_settings/indoor_fire.dat";
	char output_filename[] = "output/";

	SimulationContext sim_ctx{};
	if (argn > 1) parameterFilename = "simulation_settings/" + std::string(args[1]) + ".dat";

	int res = read_parameters(parameterFilename, sim_ctx);
	if (res != 1) {
		std::cout << "Error while trying to read " << parameterFilename << ".\n";
		return -1;
	}

	if (sim_ctx.gpu_index < 0)
		init_parallel(sim_ctx);
	else
		sim_ctx.my_rank = 0;

	if (sim_ctx.my_rank == 0) {
		if (sim_ctx.gpu_index < 0)
			std::cout << "Using cpu to simulate." << std::endl;
		else
			std::cout << "Using Gpu to simulate." << std::endl;

		print_species_information(sim_ctx);
	}

	// -------------------------------------------------
	//             Start of gpu simulation
	// -------------------------------------------------
	if (sim_ctx.gpu_index >= 0) {
		if (sim_ctx.my_rank != 0) {
			std::cout << "Gpu simulation only supported for single thread simulations!" << std::endl;
			MPI_Finalize();
			return -1;
		}
		VkUtil::VkContext vk_context{ VK_QUEUE_COMPUTE_BIT };
		std::vector<const char*> extensions{};
		VkUtil::create_vk_context(extensions, sim_ctx.gpu_index, vk_context);
		{	/* this part has to be in brackets, as all objects using vk_context attributes have to be destroyed before the vk_context is dsetroyed*/
			matrix2dRM<cell_flag> flags = read_gpu_pgm("simulation_settings/" + sim_ctx.geometry + ".pgm");

			GpuGrid grid(sim_ctx, vk_context, flags);
			GpuPipelines pipelines(vk_context);
			GpuSimulationInfos gpu_sim_ctx(vk_context, sim_ctx);
			GpuSimulationContext cpu_sim_ctx;
			pipelines.init_descriptor_sets_and_commands(grid, gpu_sim_ctx, sim_ctx);

			matrix2dRM<float> output(grid._imax, grid._jmax, 0);

			std::cout << "Starting " << sim_ctx.problem << " simulation..." << std::endl << "0%" << std::flush;

			auto start = std::chrono::steady_clock::now();

			uint32_t command_index = 0;
			VkResult r;
			while (sim_ctx.t <= sim_ctx.t_end) {          /* main simulation loop*/
				VkUtil::execute_command_buffer(vk_context, pipelines.simulate_commands[command_index]);
				r = vkQueueWaitIdle(vk_context.queue);		VkUtil::check_vk_result(r);

				gpu_sim_ctx.get_simulation_context_data(&cpu_sim_ctx, sizeof(GpuSimulationContext));
				sim_ctx.t += cpu_sim_ctx.dt;
				sim_ctx.output_time += cpu_sim_ctx.dt;
				std::cout << "\r" << sim_ctx.t / sim_ctx.t_end * 100 << "%     " << std::flush;

				if (sim_ctx.output_time >= sim_ctx.dt_value) {
					write_gpu_vtkFile((output_filename + sim_ctx.problem).c_str(), sim_ctx, grid);
					sim_ctx.output_count++;
					sim_ctx.output_time -= sim_ctx.dt_value;
				}

				command_index ^= 1;
				grid.species_concentrations.swap_images();
				grid.temperature.swap_images();
			}

			std::cout << "\r100%       " << std::flush;

			std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now() - start;
			std::cout << std::endl << "Simulation done." << std::endl;
			std::cout << "Runtime: " << elapsed_seconds.count() << " seconds" << std::endl;
		}
		VkUtil::destroy_vk_context(vk_context);
	}
	// -------------------------------------------------
	//             Start of cpu simulation
	// -------------------------------------------------
	else {
		int threads;
		MPI_Comm_size(MPI_COMM_WORLD, &threads);

		if (sim_ctx.iproc * sim_ctx.jproc != threads) {
			if (sim_ctx.my_rank == 0)
				std::cout << "Abort: the number of threads specified in the .dat file does not fit the amount of threads started.\nAmount of threads started: " << threads << std::endl;
			MPI_Finalize();
			return -1;
		}

		std::cout << "Thread " << sim_ctx.my_rank << " at cell [" << sim_ctx.omg_i << "," << sim_ctx.omg_j << "] on i from " << sim_ctx.il << " to " << sim_ctx.ir << ", on j from " << sim_ctx.jb << " to " << sim_ctx.jt << std::endl;
		//std::cout << "Thread " << my_rank << " neighbours (left, top, right, bot): " << rank_l << "," << rank_t << "," << rank_r << "," << rank_b << std::endl;

		matrix2d<cell_flag> flags = read_pgm("simulation_settings/" + sim_ctx.geometry + ".pgm");

		Grid grid(sim_ctx, flags);

		matrix2d<double> u(grid._imax, grid._jmax, 0);        //instatiating a matrix of size (imax + 1)x(jmax + 1) with 0
		matrix2d<double> v(u), f(u), g(u), p(u), rs(u), t(u);

		MPI_Barrier(MPI_COMM_WORLD);

		if (sim_ctx.my_rank == 0) {
			std::cout << "Starting " << sim_ctx.problem << " simulation..." << std::endl << "0%" << std::flush;
		}


		auto start = std::chrono::steady_clock::now();

		while (sim_ctx.t <= sim_ctx.t_end) {          /* main simulation loop*/
			calculate_dt(sim_ctx, grid);

			boundaryvalues(sim_ctx, grid);

			calculate_temp(sim_ctx, grid, t);

			for (int i = 0; i < grid.species_concentrations.size(); ++i) {
				calculate_matrix(sim_ctx, grid, grid.species_concentrations[i]);
			}

			normalize_species_ratios(sim_ctx, grid);

			calculate_fg(sim_ctx, grid, f, g);

			calculate_rs(sim_ctx, f, g, rs);

			sor_loop(sim_ctx, grid, rs);

			calculate_uv(sim_ctx, grid, f, g);

			sim_ctx.t += sim_ctx.dt;
			sim_ctx.output_time += sim_ctx.dt;
			
			if (sim_ctx.my_rank == 0)
				std::cout << "\r" << sim_ctx.t / sim_ctx.t_end * 100 << "%     " << std::flush;

			if (sim_ctx.output_time >= sim_ctx.dt_value) {
				write_vtkFile((output_filename + sim_ctx.problem).c_str(), sim_ctx, grid);
				sim_ctx.output_count++;
				sim_ctx.output_time -= sim_ctx.dt_value;
			}
		}

		if (sim_ctx.my_rank == 0) {
			std::cout << "\r100%       " << std::flush;

			std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now() - start;
			std::cout << std::endl << "Simulation done." << std::endl;
			std::cout << "Runtime: " << elapsed_seconds.count() << " seconds" << std::endl;
		}

		// Finalize MPI: Wait until all ranks are here to safely exit
		MPI_Finalize();
	}
	return 0;
}