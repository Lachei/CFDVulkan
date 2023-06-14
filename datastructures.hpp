#include "vector"
#include <cassert>
#include "enums.hpp"
#include <string>

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#define GAS_CONSTANT 8.314

extern int TEST_CASE;

template<typename T>
using matrix = std::vector<std::vector<T>>;

template<typename T>
class matrix2d {
public:
    matrix2d() :_storage(), _width(0), _height(0) {};
    matrix2d(int w, int h, T clear_value) :_storage(w* h, clear_value), _width(w), _height(h) {};

    int width() { return _width; };
    int height() { return _height; };

    T& operator()(int i, int j) {
        assert(i >= 0 && i < _width);
        assert(j >= 0 && j < _height);
        return _storage[i * _height + j];
    };
    // addressing the matrix via normalized uv coordinates
    T& operator()(float u, float v) {
        assert(u >= 0 && u <= 1 && v >= 0 && v <= 1);
        return _storage[(int(u * _width) - int(u)) * _height + int(v * _height) - int(v)];
    };
private:
    std::vector<T> _storage;
    int _width, _height;
};


template<typename T>
class matrix2dRM {
public:
    matrix2dRM() :_storage(), _width(0), _height(0) {};
    matrix2dRM(int w, int h, T clear_value) :_storage(w* h, clear_value), _width(w), _height(h) {};

    int width() { return _width; };
    int height() { return _height; };
    T* data() { return _storage.data(); };
    auto begin() { return _storage.begin(); };
    auto end() { return _storage.end(); };
    uint32_t size() { return _storage.size(); };

	T& operator()(int i, int j) {
		assert(i >= 0 && i < _width);
		assert(j >= 0 && j < _height);
		return _storage[i + j * _width];
	};
	// addressing the matrix via normalized uv coordinates
	T& operator()(float u, float v) {
		assert(u >= 0 && u <= 1 && v >= 0 && v <= 1);
		return _storage[(int(u * _width) - int(u)) + (int(v * _height) - int(v)) * _width];
	};
private:
	std::vector<T> _storage;
	int _width, _height;
};

struct Species{
    std::string name;
    double molar_mass;
    double initial_concentration;
};

struct MultiInserter{
    double temp;
    std::vector<double> concentrations; 
    double velocities[2];
};

struct SpeciesCoverter {
    double temp;
    std::vector<double> reduction_per_qm_sec;
    std::vector<std::vector<double>> conversion_ratio;
};


struct SimulationContext {
    double Re;                /* reynolds number   */
    double UI;                /* velocity x-direction */
    double VI;                /* velocity y-direction */
    double PI;                /* pressure */
    double TI;                /* temperature */
    double UIn;               /* inflow velocity u*/
    double VIn;               /* inflow velocity v*/
    WallTmps WTI;             /* wall temperatures*/
    double GX;                /* gravitation x-direction */
    double GY;                /* gravitation y-direction */
    double t_end;             /* end time */
    double xlength;           /* length of the domain x-dir.*/
    double ylength;           /* length of the domain y-dir.*/
    double dt;                /* time step */
    double dx;                /* length of a cell x-dir. */
    double dy;                /* length of a cell y-dir. */
    int imax;                 /* number of cells x-direction*/
    int jmax;                 /* number of cells y-direction*/
    double alpha;             /* uppwind differencing factor*/
    double omg;               /* relaxation factor */
    double tau;               /* safety factor for time step*/
    int itermax;              /* max. number of iterations for pressure per time step */
    double eps;               /* accuracy bound for pressure*/
    double dt_value;          /* time for output */
    std::string problem;      /* problem string*/
    typedef ::simulation simulation;    /* parsed problem (represents the problem string as simulation enum)*/
    std::string geometry;     /* geometry file name*/
    double t = 0;             /* current time*/
    double output_time = 0;   /* timer for output*/
    int output_count = 0;     /* output count*/

    double Pr;                /* Prandtl number*/
    double beta;              /* beta for thermal*/

    int gpu_index;            /* If the gpu index is larger than -1 the user specifies the gpu which should be used for simulation*/
    int my_rank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i, omg_j, iproc, jproc;  /* variables for mpi information*/
    std::vector<Species> species;
    std::vector<MultiInserter> multi_inserters;
    std::vector<SpeciesCoverter> species_converters;
};

struct GpuSimulationContext {
    float Re, UI, VI, PI, TI, UIn, VIn, WTN, WTE, WTS, WTW, GX, GY, t_end, xlenght, ylength;
    float dt;           //old dt
    float dx, dy;
    int imax, jmax;
    uint32_t alpha;         //alpha has to be set to 0 after each iteration of simulation
    float omg, tau;
    int itermax;
    float eps, dt_value, Pr, beta;
    float sor_eps, cur_sor_eps, prev_sor_eps;
    uint32_t divider, lock, sor_counter;
    uint32_t new_dt;        // new dt is interpreted as uint to allow atomic_min. This has to be set to t_end at the end of each iteration of simulation
    int patch_x, patch_y, patch_z, patch_x_origin, patch_y_origin, patch_z_origin;
    int model, o2index;
    int amt_species;
};

#endif //DATA_STRUCTURES_H
