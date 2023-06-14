#include "grid.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

Grid::Grid(SimulationContext& sim_ctx, matrix2d<cell_flag>& cell_flags)
    :_imax(sim_ctx.ir - sim_ctx.il + 2), _jmax(sim_ctx.jt - sim_ctx.jb + 2)
{

    // intstantiating the grid cells
    velocity[int(velocity_type::U)] = matrix2d<double>(_imax, _jmax, sim_ctx.UI);
    velocity[int(velocity_type::V)] = matrix2d<double>(_imax, _jmax, sim_ctx.VI);
    pressure = matrix2d<double>(_imax, _jmax, sim_ctx.PI);
    temperature = matrix2d<double>(_imax, _jmax, sim_ctx.TI);
    density = matrix2d<double>(_imax, _jmax, 0);
    neighbourhood = matrix2d<Neighbourhood>(_imax, _jmax, {});
    cell_flagg = matrix2d<cell_flag>(_imax, _jmax, cell_flag::FLUID);
    species_concentrations.resize(sim_ctx.species.size());
    for (int i = 0; i < species_concentrations.size(); ++i) {
        species_concentrations[i] = matrix2d<double>(_imax, _jmax, sim_ctx.species[i].initial_concentration);
    }
    species = sim_ctx.species;
    //filling in the flags
    for (int it = 0; it < _imax; ++it) {
        for (int jt = 0; jt < _jmax; ++jt) {
            int i = std::max(sim_ctx.il + it - 1, 0), j = std::max(sim_ctx.jb + jt - 1, 0), id = sim_ctx.imax + 1, jd = sim_ctx.jmax + 1;
            cell_flagg(it, jt) = cell_flags(float(i) / (id), float(j) / (jd));
            if (cell_flagg(it, jt) != cell_flag::FLUID && (int)cell_flagg(it, jt) <= 9) {
                pressure(it, jt) = 0;
                velocity[int(velocity_type::U)](it, jt) = 0;
                velocity[int(velocity_type::V)](it, jt) = 0;
                if (j < jd && (cell_flags(float(i) / (id), float(j + 1) / (jd)) == cell_flag::FLUID || (int)cell_flags(float(i) / (id), float(j + 1) / (jd)) > 9)) neighbourhood(it, jt).val = neighbourhood(it, jt).val | B_N;
                if (j > 0 && (cell_flags(float(i) / (id), float(j - 1) / (jd)) == cell_flag::FLUID || (int)cell_flags(float(i) / (id), float(j - 1) / (jd)) > 9)) neighbourhood(it, jt).val = neighbourhood(it, jt).val | B_S;
                if (i < id && (cell_flags(float(i + 1) / (id), float(j) / (jd)) == cell_flag::FLUID || (int)cell_flags(float(i + 1) / (id), float(j) / (jd)) > 9)) neighbourhood(it, jt).val = neighbourhood(it, jt).val | B_E;
                if (i > 0 && (cell_flags(float(i - 1) / (id), float(j) / (jd)) == cell_flag::FLUID || (int)cell_flags(float(i - 1) / (id), float(j) / (jd)) > 9)) neighbourhood(it, jt).val = neighbourhood(it, jt).val | B_W;
            }
        }
    }
};

void Grid::print_velocity(velocity_type type) {
    if (!(type == velocity_type::U || type == velocity_type::V)) {
        std::cerr << "Wrong velocity type" << std::endl;
    }
    else {
        for (int y = _jmax - 1; y >= 1; y--) {
            for (int x = 1; x < _imax; x++) {
                // Accessing velocity
                std::cout << velocity[int(type)](x, y) << " ";
            }
            // Print new line
            std::cout << std::endl;
        }
    }
}


void Grid::print_pressure() {
    for (int y = _jmax - 1; y >= 0; y--) {
        for (int x = 0; x < _imax; x++) {
            // Accessing pressure
            std::cout << pressure(x, y) << " ";
        }
        // Print new line
        std::cout << std::endl;
    }
}

void Grid::print_cell_flags()
{
    for (int y = _jmax - 1; y >= 1; y--) {
        for (int x = 1; x < _imax; x++) {
            // Accessing pressure
            cell_flag c = cell_flagg(x, y);
            switch (c) {
            case cell_flag::FLUID: std::cout << "F"; break;
            case cell_flag::FREE_SLIP: std::cout << "S"; break;
            case cell_flag::INFLOW: std::cout << "I"; break;
            case cell_flag::NO_SLIP: std::cout << "N"; break;
            case cell_flag::OUTFLOW: std::cout << "O"; break;
            case cell_flag::MULTISPECIES_0: std::cout << "M"; break;
            case cell_flag::MULTISPECIES_1: std::cout << "M"; break;
            case cell_flag::MULTISPECIES_2: std::cout << "M"; break;
            case cell_flag::MULTISPECIES_3: std::cout << "M"; break;
            case cell_flag::MULTISPECIES_4: std::cout << "M"; break;
            case cell_flag::SPECIES_CONVERTER_0: std::cout << "T"; break;
            case cell_flag::SPECIES_CONVERTER_1: std::cout << "T"; break;
            case cell_flag::SPECIES_CONVERTER_2: std::cout << "T"; break;
            case cell_flag::SPECIES_CONVERTER_3: std::cout << "T"; break;
            case cell_flag::SPECIES_CONVERTER_4: std::cout << "T"; break;
            }
        }
        // Print new line
        std::cout << std::endl;
    }
}

void Grid::print_neighbourhood()
{
    for (int y = _jmax - 1; y >= 1; y--) {
        for (int x = 1; x < _imax; x++) {
            // Accessing pressure
            switch (neighbourhood(x, y).val) {
            case B_N: std::cout << "NN "; break;
            case B_E: std::cout << "EE "; break;
            case B_S: std::cout << "SS "; break;
            case B_W: std::cout << "WW "; break;
            case B_N + B_E: std::cout << "NE "; break;
            case B_N + B_W: std::cout << "NW "; break;
            case B_S + B_E: std::cout << "SE "; break;
            case B_S + B_W: std::cout << "SW "; break;
            default: std::cout << "   "; break;
            }
        }
        // Print new line
        std::cout << std::endl;
    }
}


GpuGrid::GpuGrid(SimulationContext& sim_ctx, VkUtil::VkContext& vk_context, matrix2dRM<cell_flag>& flags)
    :_imax(sim_ctx.imax + 2), _jmax(sim_ctx.jmax + 2),
    // intstantiating the grid cells
    velocity(vk_context, _imax, _jmax, VK_FORMAT_R32G32_SFLOAT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    pressure(vk_context, _imax, _jmax, VK_FORMAT_R32_SFLOAT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    temperature(vk_context, _imax, _jmax, VK_FORMAT_R32_SFLOAT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    density(vk_context, _imax, _jmax, VK_FORMAT_R32_SFLOAT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    cell_flags(vk_context, _imax, _jmax, VK_FORMAT_R8_UINT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    neighbourhood(vk_context, _imax, _jmax, VK_FORMAT_R8_UINT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    species_concentrations(vk_context, _imax, _jmax, sim_ctx.species.size(), VK_FORMAT_R32_SFLOAT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    fg(vk_context, _imax, _jmax, VK_FORMAT_R32G32_SFLOAT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0),
    rs(vk_context, _imax, _jmax, VK_FORMAT_R32_SFLOAT, VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT, 0)
{
    // uploading data
    matrix2dRM<uint8_t> neigh(_imax, _jmax, {});
    matrix2dRM<cell_flag> cell_f(_imax, _jmax, cell_flag::FLUID);
    matrix2dRM<float> tmp(_imax, _jmax, sim_ctx.PI);
    matrix2dRM<std::pair<float, float>> vel(_imax, _jmax, { sim_ctx.UI, sim_ctx.VI });
	for (int i = 0; i < _imax; ++i) {
		for (int j = 0; j < _jmax; ++j) {
			cell_f(i, j) = flags(float(i) / (_imax - 1), float(j) / (_jmax - 1));
			if (cell_f(i, j) != cell_flag::FLUID && (int)cell_f(i, j) <= 9) {
                if (j < _jmax - 1 && (flags(float(i) / (_imax - 1), float(j + 1) / (_jmax - 1)) == cell_flag::FLUID || (int)flags(float(i) / (_imax - 1), float(j + 1) / (_jmax - 1)) > 9)) neigh(i, j) = neigh(i,j) | B_N;
                if(j>0 && (flags(float(i) / (_imax - 1), float(j-1) / (_jmax - 1))==cell_flag::FLUID || (int) flags(float(i) / (_imax - 1), float(j-1) / (_jmax - 1)) > 9)) neigh(i, j) = neigh(i, j) | B_S;
                if(i< _imax - 1 && (flags(float(i+1) / (_imax - 1), float(j) / (_jmax - 1))==cell_flag::FLUID||(int) flags(float(i+1) / (_imax - 1), float(j) / (_jmax - 1)) > 9)) neigh(i, j) = neigh(i, j) | B_W;
                if(i>0 && (flags(float(i-1) / (_imax - 1), float(j) / (_jmax - 1))==cell_flag::FLUID || (int) flags(float(i-1) / (_imax - 1), float(j) / (_jmax - 1)) > 9)) neigh(i, j) = neigh(i, j) | B_E;
			}
		}
	}
    cell_flagg = cell_f;
    neighbours = neigh;
    
    velocity.set_image(vel.data(), vel.size() * sizeof(std::pair<float, float>));
    pressure.set_image(tmp.data(), tmp.size() * sizeof(float));
    tmp = matrix2dRM<float>(_imax, _jmax, 0);
    density.set_image(tmp.data(), tmp.size() * sizeof(float));
    neighbourhood.set_image(neigh.data(), neigh.size() * sizeof(uint8_t));
    cell_flags.set_image(cell_f.data(), cell_f.size() * sizeof(cell_flag));
    tmp = matrix2dRM<float>(_imax, _jmax, sim_ctx.TI);
    temperature.set_cur_image(tmp.data(), tmp.size() * sizeof(float));
    temperature.swap_images();              // the data should be in the back image, as the values are read out from the back image and stored in the front image    
    temperature.set_cur_image(tmp.data(), tmp.size() * sizeof(float));
    for (int i = 0; i < sim_ctx.species.size(); ++i) {
        tmp = matrix2dRM<float>(_imax, _jmax, sim_ctx.species[i].initial_concentration);
        species_concentrations.set_cur_image(i, tmp.data(), tmp.size() * sizeof(float));
    }
    species_concentrations.swap_images();   // the data should be in the back image, as the values are read out from the back image and stored in the front image 
    for (int i = 0; i < sim_ctx.species.size(); ++i) {
        tmp = matrix2dRM<float>(_imax, _jmax, sim_ctx.species[i].initial_concentration);
        species_concentrations.set_cur_image(i, tmp.data(), tmp.size() * sizeof(float));
    }
};