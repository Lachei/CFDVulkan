#include "par.hpp"

void init_parallel(SimulationContext& sim_ctx)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &sim_ctx.my_rank);

	int imax = sim_ctx.imax + 2;
	int jmax = sim_ctx.jmax + 2;
	
	sim_ctx.omg_i = sim_ctx.my_rank % sim_ctx.iproc;
	sim_ctx.omg_j = sim_ctx.my_rank / sim_ctx.iproc;
	

	sim_ctx.il = (double(imax) / sim_ctx.iproc) * (sim_ctx.omg_i);
	sim_ctx.ir = (double(imax) / sim_ctx.iproc) * (sim_ctx.omg_i + 1) - 1;
	sim_ctx.jb = (double(jmax) / sim_ctx.jproc) * (sim_ctx.omg_j);
	sim_ctx.jt = (double(jmax) / sim_ctx.jproc) * (sim_ctx.omg_j + 1) - 1;

	sim_ctx.rank_l = sim_ctx.my_rank - 1;
	sim_ctx.rank_r = sim_ctx.my_rank + 1;
	sim_ctx.rank_b = sim_ctx.my_rank - sim_ctx.iproc;
	sim_ctx.rank_t = sim_ctx.my_rank + sim_ctx.iproc;	

	if (sim_ctx.omg_i == 0) sim_ctx.rank_l = MPI_PROC_NULL;
	if (sim_ctx.omg_i == sim_ctx.iproc - 1) sim_ctx.rank_r = MPI_PROC_NULL;
	if (sim_ctx.omg_j == 0) sim_ctx.rank_b = MPI_PROC_NULL;
	if (sim_ctx.omg_j == sim_ctx.jproc - 1) sim_ctx.rank_t = MPI_PROC_NULL;
}


void matrix_comm(matrix2d<double>& M, int rank_l, int rank_r, int rank_b, int rank_t, std::vector<double>& bufSend, std::vector<double>& bufRecv, MPI_Status* status, int chunk)
{
	// send/reseive left/right
	bufSend.resize(M.height());
	bufRecv.resize(M.height());

	//send left - receive right
	int c = 0;
	for (int j = 0; j < M.height(); ++j) {
		bufSend[c++] = M(2,j);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_l, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_r != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int j = 0; j < M.height(); ++j) {
			M(M.width() - 1, j) = bufRecv[c++];
		}
	}

	//send right - receive left
	c = 0;
	for (int j = 0; j < M.height(); ++j) {
		bufSend[c++] = M(M.width() - 2, j);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_r, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_l != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int j = 0; j < M.height(); ++j) {
			M(1,j) = bufRecv[c++];
		}
	}

	// send/reseive top/bottom
	bufSend.resize(M.width());
	bufRecv.resize(M.width());

	//send top - receive bottom
	c = 0;
	for (int i = 0; i < M.width(); ++i) {
		bufSend[c++] = M(i,M.height() - 2);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_t, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_b != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int i = 0; i < M.width(); ++i) {
			M(i,1) = bufRecv[c++];
		}
	}

	//send bottom - receive top
	c = 0;
	for (int i = 0; i < M.width(); ++i) {
		bufSend[c++] = M(i,2);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_b, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_t != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int i = 0; i < M.width(); ++i) {
			M(i,M.height() - 1) = bufRecv[c++];
		}
	}
}

void matrix_pair_comm(matrix2d<double>& M, matrix2d<double>& N, int rank_l, int rank_r, int rank_b, int rank_t, std::vector<double>& bufSend, std::vector<double>& bufRecv, MPI_Status* status, int chunk)
{
	// send/reseive left/right
	//bufSend.resize(2 * (jt - jb + 1));
	//bufRecv.resize(2 * (jt - jb + 1));
    //
	////send left - receive right
	//int c = 0;
	//for (int j = jb; j < jt + 1; ++j) {
	//	bufSend[c++] = M[il][j];
	//	bufSend[c++] = N[il][j];
	//}
	//if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_l, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
	//	std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
	//	exit(-1);
	//}
	//c = 0;
	//if (rank_r != MPI_PROC_NULL) {         /* only copy if ther was something received*/
	//	for (int j = jb; j < jt + 1; ++j) {
	//		M[ir + 1][j] = bufRecv[c++];
	//		N[ir + 1][j] = bufRecv[c++];
	//	}
	//}
    //
	////send right - receive left
	//c = 0;
	//for (int j = jb; j < jt + 1; ++j) {
	//	bufSend[c++] = M[ir][j];
	//	bufSend[c++] = N[ir][j];
	//}
	//if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_r, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
	//	std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
	//	exit(-1);
	//}
	//c = 0;
	//if (rank_l != MPI_PROC_NULL) {         /* only copy if ther was something received*/
	//	for (int j = jb; j < jt + 1; ++j) {
	//		M[il - 1][j] = bufRecv[c++];
	//		N[il - 1][j] = bufRecv[c++];
	//	}
	//}
    //
	//// send/reseive top/bottom
	//bufSend.resize(2 * (ir - il + 1));
	//bufRecv.resize(2 * (ir - il + 1));
    //
	////send top - receive bottom
	//c = 0;
	//for (int i = il; i < ir + 1; ++i) {
	//	bufSend[c++] = M[i][jt];
	//	bufSend[c++] = N[i][jt];
	//}
	//if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_t, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
	//	std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
	//	exit(-1);
	//}
	//c = 0;
	//if (rank_b != MPI_PROC_NULL) {         /* only copy if ther was something received*/
	//	for (int i = il; i < ir + 1; ++i) {
	//		M[i][jb - 1] = bufRecv[c++];
	//		N[i][jb - 1] = bufRecv[c++];
	//	}
	//}
    //
	////send bottom - receive top
	//c = 0;
	//for (int i = il; i < ir + 1; ++i) {
	//	bufSend[c++] = M[i][jb];
	//	bufSend[c++] = N[i][jb];
	//}
	//if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_b, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
	//	std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
	//	exit(-1);
	//}
	//c = 0;
	//if (rank_t != MPI_PROC_NULL) {         /* only copy if ther was something received*/
	//	for (int i = il; i < ir + 1; ++i) {
	//		M[i][jt + 1] = bufRecv[c++];
	//		N[i][jt + 1] = bufRecv[c++];
	//	}
	//}
}

void uv_comm(matrix2d<double>& U, matrix2d<double>& V, int rank_l, int rank_r, int rank_b, int rank_t, std::vector<double>& bufSend, std::vector<double>& bufRecv, MPI_Status* status, int chunk)
{
	// send/reseive left/right
	bufSend.resize(2 * U.height());
	bufRecv.resize(2 * U.height());

	//send left - receive right
	int c = 0;
	for (int j = 0; j < U.height(); ++j) {
		bufSend[c++] = U(2,j);
		bufSend[c++] = V(2,j);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_l, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_r != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int j = 0; j < U.height(); ++j) {
			U(U.width() - 1,j) = bufRecv[c++];
			V(U.width() - 1,j) = bufRecv[c++];
		}
	}

	//send right - receive left
	bufSend.resize(4 * U.height());
	bufRecv.resize(4 * U.height());
	c = 0;
	for (int j = 0; j < U.height(); ++j) {
        bufSend[c++] = U(U.width() - 2, j);
		bufSend[c++] = U(U.width() - 3, j);
		bufSend[c++] = V(V.width() - 2, j);
		bufSend[c++] = V(V.width() - 3, j);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_r, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_l != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int j = 0; j < U.height(); ++j) {
			U(1,j) = bufRecv[c++];
			U(0,j) = bufRecv[c++];
			V(1,j) = bufRecv[c++];
			V(0,j) = bufRecv[c++];
		}
	}

	// send/reseive top/bottom
	bufSend.resize(4 * U.width());
	bufRecv.resize(4 * U.width());

	//send top - receive bottom
	c = 0;
	for (int i = 0; i < U.width(); ++i) {
        bufSend[c++] = U(i, U.height() - 2);
		bufSend[c++] = U(i, U.height() - 3);
		bufSend[c++] = V(i, V.height() - 2);
		bufSend[c++] = V(i, V.height() - 3);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_t, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_b != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int i = 0; i < U.width(); ++i) {
			U(i,1) = bufRecv[c++];
			U(i,0) = bufRecv[c++];
			V(i,1) = bufRecv[c++];
			V(i,0) = bufRecv[c++];
		}
	}

	//send bottom - receive top
	bufSend.resize(2 * U.width());
	bufRecv.resize(2 * U.width());
	c = 0;
	for (int i = 0; i < U.width(); ++i) {
		bufSend[c++] = U(i,2);
		bufSend[c++] = V(i,2);
	}
	if (MPI_Sendrecv(bufSend.data(), bufSend.size(), MPI_DOUBLE, rank_b, chunk, bufRecv.data(), bufRecv.size(), MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status) != MPI_SUCCESS) {
		std::cout << "We do have a problem with sending/receiving.\nMPI status: " << status->MPI_ERROR << std::endl;
		exit(-1);
	}
	c = 0;
	if (rank_t != MPI_PROC_NULL) {         /* only copy if ther was something received*/
		for (int i = 0; i < U.width(); ++i) {
			U(i,U.height() - 1) = bufRecv[c++];
			V(i,U.height() - 1) = bufRecv[c++];
		}
	}
}
