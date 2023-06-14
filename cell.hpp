#ifndef CFDLAB_CELL_HPP
#define CFDLAB_CELL_HPP
#include <array>
#include "enums.hpp"

class Cell {
public:
    // Constructors
    Cell();
    Cell(double& PI, double& UI, double& VI);
	Cell(double& PI, double& UI, double& VI, double& TI);
    Cell(double& PI, double& UI, double& VI, cell_flag flag);

    // Get + Set pressure
    double& pressure();
    void set_pressure(double& value);

    // Get + Set temperature
    double& temperature();
    void set_temperature(double& value);

    // Get + Set density
    double& density();
    void set_density(double& value);

    // Get + Set velocity
    double& velocity(velocity_type type);
    void set_velocity(double& value, velocity_type type);

    // Get + Set vorder
    bool& border(border_position position);
    void set_border(border_position position);

	// Get + Set flag
	cell_flag& flag();
	void set_flag(cell_flag flag);

    Neighbourhood& get_neighbourhood();
    void add_neighbour(unsigned int n);

private:
    // one pressure value per call
    double _pressure = 0;

    //one temperature value per cell(-center)
    double _temperature=0;

    //one density value per cell for boussinesq approximation
    double _density=0;

    // Fixed size velocity
    std::array<double, 2> _velocity = {0};
    
    // Fixed number of borders
    std::array<bool, 4> _border = {false};

	// cell flag
	cell_flag _flag = cell_flag::FLUID;

	Neighbourhood _neighbourhood{};
};

#endif //CFDLAB_CELL_HPP
