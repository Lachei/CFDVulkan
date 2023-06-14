#include "sor.hpp"
#include "boundary_val.hpp"
#include "par.hpp"
#include <cmath>
#include <mpi.h>

void sor(
        double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        Grid& grid,
        matrix<double> &RS,
        double *res
) {
    int i,j;
    double rloc;
    double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
    matrix<double> P;
    //grid.pressure(P);

    /* SOR iteration */
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax; j++) {
            P.at(i).at(j) = (1.0-omg)*P.at(i).at(j)
                            + coeff*(( P.at(i+1).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j));
        }
    }

    /* compute the residual */
    rloc = 0;
    for(i = 1; i <= imax; i++) {
        for(j = 1; j <= jmax; j++) {
            rloc += ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j))*
                    ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j));
        }
    }
    rloc = rloc/(imax*jmax);
    rloc = sqrt(rloc);
    /* set residual */
    *res = rloc;


    /* set boundary values */
    for(i = 1; i <= imax; i++) {
        P.at(i)[0] = P.at(i)[1];
        P.at(i).at(jmax+1) = P.at(i).at(jmax);
    }
    for(j = 1; j <= jmax; j++) {
        P[0].at(j) = P[1].at(j);
        P.at(imax+1).at(j) = P.at(imax).at(j);
    }

    //grid.set_pressure(P);

}

void sor_loop(SimulationContext& sim_ctx, Grid& grid, matrix2d<double>& RS)
{
	int i, j, counter;
	double rloc;
	double coeff = sim_ctx.omg / (2.0 * (1.0 / (sim_ctx.dx * sim_ctx.dx) + 1.0 / (sim_ctx.dy * sim_ctx.dy)));
	matrix2d<double>& P = grid.pressure;
	static std::vector<double> send, rec;
	MPI_Status mpi_stat;

	for (int it = 0; it < sim_ctx.itermax; ++it) {
		/* SOR iteration */
		for (i = 2; i < grid._imax - 1; i++) {
			for (j = 2; j < grid._jmax - 1; j++) {
                if(grid.cell_flagg(i,j)==cell_flag::FLUID || (int)grid.cell_flagg(i,j)>9){
				    P(i, j) = (1.0 - sim_ctx.omg) * P(i, j)
					    + coeff * ((P(i + 1, j) + P(i - 1, j)) / (sim_ctx.dx * sim_ctx.dx) + (P(i, j + 1) + P(i, j - 1)) / (sim_ctx.dy * sim_ctx.dy) - RS(i, j));
                }
			}
		}

		matrix_comm(P, sim_ctx.rank_l, sim_ctx.rank_r, sim_ctx.rank_b, sim_ctx.rank_t, send, rec, &mpi_stat, 0);

		/* compute the residual */
		rloc = 0;
        counter=0;
		for (i = 2; i < grid._imax - 1; i++) {
			for (j = 2; j < grid._jmax - 1; j++) {
                if(grid.cell_flagg(i,j)==cell_flag::FLUID || (int)grid.cell_flagg(i,j)>9){
				    rloc += ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (sim_ctx.dx * sim_ctx.dx) + (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (sim_ctx.dy * sim_ctx.dy) - RS(i, j)) *
				    	((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (sim_ctx.dx * sim_ctx.dx) + (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (sim_ctx.dy * sim_ctx.dy) - RS(i, j));
                    counter++;
                }
            }
		}
        /* normalize by number of fluid cells */
		/* summing up accross all processes*/
		double rcpy[2] = { rloc,double(counter) }, rec[2]{};
		MPI_Allreduce(rcpy, rec, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		rloc = rec[0] / rec[1];
		rloc = sqrt(rloc);

		/* set boundary values */
        for (i = 1; i < grid._imax; i++) {
            for (j = 1; j < grid._jmax; j++) {
				if(grid.cell_flagg(i,j) != cell_flag::FLUID && (int)grid.cell_flagg(i, j) <= 9)
					set_pressure_boundary(grid.neighbourhood(i,j),P,i,j);
            }
        }

		/* check residual */
		if (rloc < sim_ctx.eps) break;
	}
}

