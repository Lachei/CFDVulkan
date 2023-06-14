#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"
#include <cmath>
#include <stdio.h>


void setNoSlip(Neighbourhood n, Grid& grid, int i, int j){
    matrix2d<double>& u = grid.velocity[0],& v = grid.velocity[1],& p = grid.pressure;
    switch(n.val){
        case B_N:{
            v(i,j) = 0;
            u(i - 1,j) = -1 * u(i-1,j+1);
            u(i, j) = -1 * u(i, j + 1);
            p(i, j) = p(i, j + 1);
            break;
        }
        case B_S:{
            v(i,j-1) = 0;
            u(i-1,j) = -1 * u(i-1,j-1);
            u(i,j) = -1 * u(i,j-1);
            p(i,j) = p(i,j-1);
            break;
        }
        case B_W:{
            u(i-1,j) = 0;
            v(i,j-1) = -1* v(i-1,j-1);
            v(i,j) = -1 * v(i-1,j);
            p(i,j) = p(i-1,j);
            break;
        }
        case B_E:{
            u(i,j) = 0;
            v(i, j - 1) = -1* v(i+1,j-1);
            v(i,j) = -1 * v(i+1,j);
            p(i,j) = p(i+1,j);
            break;
        }
        case(B_N+B_E):{
            u(i, j) = 0;
            v(i, j) = 0;
            u(i - 1, j) = -1*u(i-1,j+1);
            v(i, j - 1) = -1*v(i+1,j-1);
            p(i,j) = .5*(p(i,j+1)+p(i+1,j));
            break;
        }
        case(B_S+B_E):{
            u(i, j) = 0;
            v(i,j-1) = 0;
            u(i - 1, j) = -1*u(i-1,j-1);
            v(i,j) = -1*v(i+1,j);
            p(i,j) = .5*(p(i,j-1)+p(i+1,j));
            break;
        }
        case(B_N+B_W):{
            v(i,j) = 0;
            u(i-1,j) = 0;
            u(i,j) = -1*u(i,j+1);
            v(i, j - 1) = -1*v(i-1,j-1);
            p(i,j)= .5*(p(i,j+1)+p(i-1,j));
            break;
        }
        case(B_S+B_W):{
            v(i,j-1) = 0;
            u(i-1,j) = 0;
            u(i,j) = -1*u(i,j-1);
            v(i,j)=-1*v(i-1,j);
            p(i,j)= .5*(p(i,j-1)+p(i-1,j));
            break;
        }
    }
}

void setFreeSlip(Neighbourhood n, Grid& grid, int i, int j){
    matrix2d<double>& u = grid.velocity[0],& v = grid.velocity[1],& p = grid.pressure;
    switch(n.val){
        case B_N:{
            v(i, j) = 0;
            u(i - 1, j) = u(i - 1, j + 1);
            u(i, j) = u(i, j + 1);
            p(i, j) = p(i, j + 1);
            break;
        }
        case B_S:{
            v(i, j - 1) = 0;
            u(i - 1, j) = u(i - 1, j - 1);
            u(i, j) = u(i, j - 1);
            p(i, j) = p(i, j - 1);
            break;
        }
        case B_W:{
            u(i - 1, j) = 0;
            v(i, j - 1) = v(i - 1, j - 1);
            v(i, j) = v(i - 1, j);
            p(i, j) = p(i - 1, j);
            break;
        }
        case B_E:{
            u(i, j) = 0;
            v(i, j - 1) = v(i + 1, j - 1);
            v(i, j) = v(i + 1, j);
            p(i, j) = p(i + 1, j);
            break;
        }
        case(B_N+B_E):{
            u(i, j) = 0;
            v(i, j) = 0;
            u(i - 1, j) = u(i - 1, j + 1);
            v(i, j - 1) = v(i + 1, j - 1);
            p(i, j) = .5 * (p(i, j + 1) + p(i + 1, j));
            break;
        }
        case(B_S+B_E):{
            u(i, j) = 0;
            v(i, j - 1) = 0;
            u(i - 1, j) = u(i - 1, j - 1);
            v(i, j) = v(i + 1, j);
            p(i, j) = .5 * (p(i, j - 1) + p(i + 1, j));
            break;
        }
        case(B_N+B_W):{
            v(i, j) = 0;
            u(i - 1, j) = 0;
            u(i, j) = u(i, j + 1);
            v(i, j - 1) = v(i - 1, j - 1);
            p(i, j) = .5 * (p(i, j + 1) + p(i - 1, j));
            break;
        }
        case(B_S+B_W):{
            v(i, j - 1) = 0;
            u(i - 1, j) = 0;
            u(i, j) = u(i, j - 1);
            v(i, j) = v(i - 1, j);
            p(i, j) = .5 * (p(i, j - 1) + p(i - 1, j));
            break;
        }
    }
}

void setInFlow(Neighbourhood n, Grid& grid, int i, int j, double u_in, double v_in){
    matrix2d<double>& u = grid.velocity[0],& v = grid.velocity[1],& p = grid.pressure;
    std::vector<matrix2d<double>>& s = grid.species_concentrations;
    for (int c = 0; c < s.size(); ++c) {
        s[c](i, j) = grid.species[c].initial_concentration;
    }
    switch(n.val){
        case B_N:{
            v(i, j) = v_in;
            u(i - 1, j) = 2 * u_in - 1 * u(i - 1, j + 1);
            u(i, j) = 2 * u_in - 1 * u(i, j + 1);
            p(i, j) = p(i, j + 1);
            break;
        }
        case B_S:{
            v(i, j - 1) = v_in;
            u(i - 1, j) = 2 * u_in - 1 * u(i - 1, j - 1);
            u(i, j) = 2 * u_in - 1 * u(i, j - 1);
            p(i, j) = p(i, j - 1);
            break;
        }
        case B_W:{
            u(i - 1, j) = u_in;
            v(i, j - 1) = 2 * v_in - 1 * v(i - 1, j - 1);
            v(i, j) = 2 * v_in - 1 * v(i - 1, j);
            p(i, j) = p(i - 1, j);
            break;
        }
        case B_E:{
            u(i, j) = u_in;
            v(i, j - 1) = 2 * v_in - 1 * v(i + 1, j - 1);
            v(i, j) = 2 * v_in - 1 * v(i + 1, j);
            p(i, j) = p(i + 1, j);
            break;
        }
        case(B_N+B_E):{
            u(i, j) = u_in;
            v(i, j) = v_in;
            u(i - 1, j) = 2 * u_in - 1 * u(i - 1, j + 1);
            v(i, j - 1) = 2 * v_in - 1 * v(i + 1, j - 1);
            p(i, j) = .5 * (p(i, j + 1) + p(i + 1, j));
            break;
        }
        case(B_S+B_E):{
            u(i, j) = u_in;
            v(i, j - 1) = v_in;
            u(i - 1, j) = 2 * u_in - 1 * u(i - 1, j - 1);
            v(i, j) = 2 * v_in - 1 * v(i + 1, j);
            p(i, j) = .5 * (p(i, j - 1) + p(i + 1, j));
            break;
        }
        case(B_N+B_W):{
            v(i, j) = v_in;
            u(i - 1, j) = u_in;
            u(i, j) = 2 * u_in - 1 * u(i, j + 1);
            v(i, j - 1) = 2 * v_in - 1 * v(i - 1, j - 1);
            p(i, j) = .5 * (p(i, j + 1) + p(i - 1, j));
            break;
        }
        case(B_S+B_W):{
            v(i, j - 1) = v_in;
            u(i - 1, j) = u_in;
            u(i, j) = 2 * u_in - 1 * u(i, j - 1);
            v(i, j) = 2 * v_in - 1 * v(i - 1, j);
            p(i, j) = .5 * (p(i, j - 1) + p(i - 1, j));
            break;
        }
    }
}

void setOutFlow(Neighbourhood n, Grid& grid, int i, int j){
    matrix2d<double>& u = grid.velocity[0],& v = grid.velocity[1],& p = grid.pressure;
    std::vector<matrix2d<double>>& s = grid.species_concentrations;
    switch(n.val){
        case B_N:{
            v(i, j) = v(i, j + 1);
            u(i - 1, j) = u(i - 1, j + 1);
            u(i, j) = u(i, j + 1);
            p(i, j) = p(i, j + 1);
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = s[c](i, j + 1);
            }
            break;
        }
        case B_S:{
            v(i, j - 1) = v(i, j - 2);
            u(i, j - 1) = u(i, j - 2);
            v(i,j) = v(i, j - 1);
            u(i, j) = u(i, j - 1);
            p(i, j) = p(i, j - 1);
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = s[c](i, j - 1);
            }
            break;
        }
        case B_W:{
            u(i - 1, j) = u(i - 2, j);
            v(i - 1, j) = v(i - 2, j);
            u(i, j) = u(i - 1, j);
            v(i, j) = v(i - 1, j);
            p(i, j) = p(i - 1, j);
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = s[c](i - 1, j);
            }
            break;
        }
        case B_E:{
            u(i, j) = u(i + 1, j);
            v(i, j - 1) =  v(i + 1, j - 1);
            v(i, j) = v(i + 1, j);
            p(i, j) = p(i + 1, j);
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = s[c](i + 1, j);
            }
            break;
        }
        case(B_N+B_E):{
            u(i, j) = u(i + 1, j);
            v(i, j) = v(i, j + 1);
            u(i - 1, j) = u(i - 1, j + 1);
            v(i, j - 1) = v(i + 1, j - 1);
            p(i, j) = .5 * (p(i, j + 1) + p(i + 1, j));
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = .5 * (s[c](i, j + 1) + s[c](i + 1, j));
            }
            break;
        }
        case(B_S+B_E):{
            u(i, j) = u(i - 1, j);
            v(i, j - 1) = v(i, j - 2);
            u(i - 1, j) = u(i - 1, j - 1);
            v(i, j) = v(i + 1, j);
            p(i, j) = .5 * (p(i, j - 1) + p(i + 1, j));
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = .5 * (s[c](i, j - 1) + s[c](i + 1, j));
            }
            break;
        }
        case(B_N+B_W):{
            v(i, j) = v(i, j + 1);
            u(i - 1, j) = (i - 2, j);
            u(i, j) = u(i, j + 1);
            v(i, j - 1) = v(i - 1, j - 1);
            p(i, j) = .5 * (p(i, j + 1) + p(i - 1, j));
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = .5 * (s[c](i, j + 1) + s[c](i - 1, j));
            }
            break;
        }
        case(B_S+B_W):{
            v(i, j - 1) = v(i, j - 2);
            u(i - 1, j) = u(i - 2, j);
            u(i, j) = u(i, j - 1);
            v(i, j) = v(i - 1, j);
            p(i, j) = .5 * (p(i, j - 1) + p(i - 1, j));
            for (int c = 0; c < s.size(); ++c) {
                s[c](i, j) = .5 * (s[c](i, j - 1) + s[c](i - 1, j));
            }
            break;
		}
    }
}

void setAdiabaticMultispecies(SimulationContext& sim_ctx, Grid& grid, int i, int j) {
    matrix2d<double>& u = grid.velocity[0], & v = grid.velocity[1], & p = grid.pressure, &t = grid.temperature;
    std::vector<matrix2d<double>>& s = grid.species_concentrations;
    switch (grid.neighbourhood(i,j).val) {
    case B_N: {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i,j)) - int(cell_flag::MULTISPECIES_0)].temp - t(i, j + 1);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - s[sp](i, j + 1);
        }
        break;
    }
    case B_S: {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].temp - t(i, j - 1);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - s[sp](i, j - 1);
        }
        break;
    }
    case B_W: {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].temp - t(i - 1, j);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - s[sp](i - 1, j);
        }
        break;
    }
    case B_E: {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].temp - t(i + 1, j);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - s[sp](i + 1, j);
        }
        break;
    }
    case(B_N + B_E): {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].temp - .5 * t(i, j + 1) - .5* t(i + 1, j);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - .5 * s[sp](i, j + 1) - .5 * s[sp](i + 1, j);
        }
        break;
    }
    case(B_S + B_E): {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].temp - .5 * t(i, j - 1) - .5 * t(i + 1, j);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - .5 * s[sp](i, j - 1) - .5 * s[sp](i + 1, j);
        }
        break;
    }
    case(B_N + B_W): {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].temp - .5 * t(i, j + 1) - .5 * t(i - 1, j);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - .5 * s[sp](i, j + 1) - .5 * s[sp](i - 1, j);
        }
        break;
    }
    case(B_S + B_W): {
        t(i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].temp - .5 * t(i, j - 1) - .5 * t(i - 1, j);
        for (int sp = 0; sp < s.size(); ++sp) {
            s[sp](i, j) = 2 * sim_ctx.multi_inserters[int(grid.cell_flagg(i, j)) - int(cell_flag::MULTISPECIES_0)].concentrations[sp] - .5 * s[sp](i, j - 1) - .5 * s[sp](i - 1, j);
        }
        break;
    }
    }
}

void setBoundaryCondition(int i, int j, SimulationContext& sim_ctx, Grid& grid){
    //if(i==0)printf("Setting cell[%d][%d]\n",i,j);
    switch(grid.cell_flagg(i,j)){
        case cell_flag::NO_SLIP: { setNoSlip(grid.neighbourhood(i,j), grid, i, j); break;}
        case cell_flag::INFLOW: { setInFlow(grid.neighbourhood(i, j), grid, i, j, sim_ctx.UIn, sim_ctx.VIn); break;}
        case cell_flag::FREE_SLIP: { setFreeSlip(grid.neighbourhood(i, j), grid, i, j); break;}
        case cell_flag::OUTFLOW: { setOutFlow(grid.neighbourhood(i, j), grid, i, j); break;}
        case cell_flag::MULTISPECIES_0: { setInFlow(grid.neighbourhood(i, j), grid, i, j, sim_ctx.multi_inserters[0].velocities[0], sim_ctx.multi_inserters[0].velocities[1]); break; }
        case cell_flag::MULTISPECIES_1: { setInFlow(grid.neighbourhood(i, j), grid, i, j, sim_ctx.multi_inserters[1].velocities[0], sim_ctx.multi_inserters[1].velocities[1]); break; }
        case cell_flag::MULTISPECIES_2: { setInFlow(grid.neighbourhood(i, j), grid, i, j, sim_ctx.multi_inserters[2].velocities[0], sim_ctx.multi_inserters[2].velocities[1]); break; }
        case cell_flag::MULTISPECIES_3: { setInFlow(grid.neighbourhood(i, j), grid, i, j, sim_ctx.multi_inserters[3].velocities[0], sim_ctx.multi_inserters[3].velocities[1]); break; }
        case cell_flag::MULTISPECIES_4: { setInFlow(grid.neighbourhood(i, j), grid, i, j, sim_ctx.multi_inserters[4].velocities[0], sim_ctx.multi_inserters[4].velocities[1]); break; }
        case cell_flag::FLUID:{
            return;
        }
        default: {
            return;
        }
    }
}


void setAdiabaticTemp(Grid& grid, int i, int j){
    switch (grid.neighbourhood(i,j).val)
    {
    case B_N: {
        grid.temperature(i,j) = grid.temperature(i,j+1);
        break;
    }
    case B_S: {
        grid.temperature(i,j) = grid.temperature(i,j-1);
        break;
    }
    case B_W: {
        grid.temperature(i,j) = grid.temperature(i-1,j);
        break;
    }
    case B_E: {
        grid.temperature(i,j) = grid.temperature(i+1,j);
        break;
    }
    case (B_N+B_E): {
        grid.temperature(i,j) =.5*(grid.temperature(i,j+1)+grid.temperature(i+1,j));
        break;
    }
    case (B_N+B_W): {
        grid.temperature(i,j) =.5*(grid.temperature(i,j+1)+grid.temperature(i-1,j));
        break;
    }
    case (B_S+B_E): {
        grid.temperature(i,j) =.5*(grid.temperature(i,j-1)+grid.temperature(i+1,j));
        break;
    }
    case (B_S+B_W): {
        grid.temperature(i,j) =.5*(grid.temperature(i,j-1)+grid.temperature(i-1,j));
        break;
    }
    
    default:
        break;
    }
}

void set_conversion_temp(int model, SimulationContext &sim_ctx, Grid &grid, int converter_index, int i, int j){
    switch(model){
        //Case 2 = FIRE
        case 2: {
            int o2_index=-1;
            for(int i=0; i<sim_ctx.species.size(); i++){
                if(sim_ctx.species[i].name.compare("O2")==0){
                    o2_index=i;
                    break;
                }
            }
            //printf("Setting temp\n");
            grid.temperature(i,j)= grid.species_concentrations[o2_index](i,j)>=0.16 ? sim_ctx.species_converters[converter_index].temp : grid.temperature(i,j); 
            return;
        }
        default: {
            grid.temperature(i,j)=sim_ctx.species_converters[converter_index].temp;
        }
    }
}

//assuming i/jmax to refer to size of inner grid
void boundaryvalues(SimulationContext& sim_ctx, Grid& grid) {

    //horizontal boundaries
    if (sim_ctx.jb == 0)
        for (int i = 1; i < grid._imax-1; i++) {
            if (grid.cell_flagg(i, 1) != cell_flag::FLUID && (int)grid.cell_flagg(i, 1) <= 9)setBoundaryCondition(i, 1, sim_ctx, grid);
            if (std::isnan(sim_ctx.WTI.N)) {
                grid.temperature(i, 1) = grid.temperature(i, 2);
            }
            else {
                grid.temperature(i,1) = 2 * sim_ctx.WTI.N - grid.temperature(i, 2);
            }
        }

     if (sim_ctx.jt == sim_ctx.jmax + 1)
         for (int i = 1; i < grid._imax-1; i++) {
             if (grid.cell_flagg(i, grid._jmax - 1) != cell_flag::FLUID && (int)grid.cell_flagg(i, grid._jmax - 1) <=9)setBoundaryCondition(i, grid._jmax - 1, sim_ctx, grid);
             if (std::isnan(sim_ctx.WTI.S)) {
                 grid.temperature(i, grid._jmax - 1) = grid.temperature(i, grid._jmax - 2);
             }
             else {
                 grid.temperature(i,grid._jmax - 1)  = 2 * sim_ctx.WTI.S - grid.temperature(i, grid._jmax - 2);
             }
         }

    //vertical boundaries
	if(sim_ctx.il == 0)
        for (int j = 1; j < grid._jmax-1; j++) {
            if (grid.cell_flagg(1, j) != cell_flag::FLUID && (int)grid.cell_flagg(1,j)<=9)setBoundaryCondition(1, j, sim_ctx, grid);
            if (std::isnan(sim_ctx.WTI.E)) {
                grid.temperature(1, j) = grid.temperature(2, j);
            }
            else {
                grid.temperature(1, j) = 2 * sim_ctx.WTI.E - grid.temperature(2, j);
            }
        }

    if (sim_ctx.ir == sim_ctx.imax + 1)
        for (int j = 1; j < grid._jmax-1; j++) {
			if (grid.cell_flagg(grid._imax - 1, j) != cell_flag::FLUID && (int)grid.cell_flagg(grid._imax-1,j)<=9)setBoundaryCondition(grid._imax - 1, j, sim_ctx, grid);
			if (std::isnan(sim_ctx.WTI.W)) {
                grid.temperature(grid._imax - 1, j) = grid.temperature(grid._imax - 2, j);
			}
			else {
                grid.temperature(grid._imax - 1, j) = 2 * sim_ctx.WTI.W - grid.temperature(grid._imax - 2, j);
			}
		}
    
	//Obstacles && converters
    for (int i = 2; i < grid._imax - 1; i++) {
        for (int j = 2; j < grid._jmax - 1; j++) {
            if (grid.cell_flagg(i, j) != cell_flag::FLUID && (int)grid.cell_flagg(i, j) <=9) {
                setBoundaryCondition(i, j, sim_ctx, grid);
                if (grid.cell_flagg(i, j) == cell_flag::FREE_SLIP || grid.cell_flagg(i,j) == cell_flag::NO_SLIP || grid.cell_flagg(i,j) == cell_flag::OUTFLOW)
                    setAdiabaticTemp(grid, i, j);
                else if (grid.cell_flagg(i,j) == cell_flag::INFLOW ) {
                    grid.temperature(i,j) = sim_ctx.TI;
                }
                else
                    setAdiabaticMultispecies(sim_ctx, grid, i, j);
            }else if((int)grid.cell_flagg(i,j)>9){
                if(!std::isnan(sim_ctx.species_converters[(int)grid.cell_flagg(i,j)-10].temp))set_conversion_temp(TEST_CASE,sim_ctx,grid,(int)grid.cell_flagg(i,j)-10,i , j);
            }
        }
    }
    //printf("Cell [0][17] is %f\n",grid.cell(0,17).velocity(velocity_type::U));
}



void set_fg_boundary(Grid &g, matrix2d<double> &F, matrix2d<double> &G, int i, int j){
    matrix2d<double>& u = g.velocity[0],& v = g.velocity[1];
    switch (g.neighbourhood(i,j).val)
    {
    case B_N: {
		G(i, j) = v(i, j);
        break;
    }
    case B_S: {
		G(i - 1,j) = v(i - 1, j);
        break;
    }
    case B_W: {
		F(i - 1,j) = u(i - 1, j);
        break;
    }
    case B_E: {
		F(i, j) = u(i,j);
        break;
    }
    case (B_N+B_E): {
		G(i, j) = v(i, j);
		F(i, j) = u(i, j);
        break;
    }
    case (B_N+B_W): {
		F(i - 1,j) = u(i - 1, j);
		G(i, j) = v(i, j);
        break;
    }
    case (B_S+B_E): {
		G(i,j - 1) = v(i, j - 1);
		F(i, j) = u(i, j);
        break;
    }
    case (B_S+B_W): {
		G(i,j - 1) = v(i, j - 1);
		F(i - 1,j) = u(i - 1, j);
        break;
    }
    
    default:
        break;
    }
}

void set_pressure_boundary(Neighbourhood n,matrix2d<double> &P, int i, int j){

    double zero=.0; 
    switch(n.val){
        case B_N:{
            P(i, j)=P(i, j + 1);
            break;
        }
        case B_S:{
            P(i, j)=P(i,j-1);
            break;
        }
        case B_W:{
            P(i, j)=P(i-1,j);
            break;
        }
        case B_E:{
            P(i, j)=P(i+1,j);
            break;
        }
        case(B_N+B_E):{
            P(i, j)= .5*(P(i, j + 1)+P(i+1,j));
            break;
        }
        case(B_S+B_E):{
            P(i, j)= .5*(P(i,j-1)+P(i+1, j));
            break;
        }
        case(B_N+B_W):{
            P(i, j)= .5*(P(i, j + 1)+P(i-1,j));
            break;
        }
        case(B_S+B_W):{
            P(i, j)= .5*(P(i,j-1)+P(i-1,j));
            break;
        }
    }
}
