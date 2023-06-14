#pragma once
#include "datastructures.hpp"
#include <vector>
#include <mpi.h>
#include <iostream>

void init_parallel(SimulationContext& sim_ctx);

void matrix_comm(matrix2d<double>& M, int rank_l, int rank_r, int rank_b, int rank_t, std::vector<double>& bufSend, std::vector<double>& bufRecv, MPI_Status* status, int chunk);

void matrix_pair_comm(matrix2d<double>& M, matrix2d<double>& N, int rank_l, int rank_r, int rank_b, int rank_t, std::vector<double>& bufSend, std::vector<double>& bufRecv, MPI_Status* status, int chunk);

void uv_comm(matrix2d<double>& U, matrix2d<double>& V, int rank_l, int rank_r, int rank_b, int rank_t, std::vector<double>& bufSend, std::vector<double>& bufRecv, MPI_Status* status, int chunk);