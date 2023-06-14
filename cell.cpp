#include "cell.hpp"

Cell::Cell() {};

Cell::Cell(double& PI, double& UI, double& VI):
_pressure(PI) {
    set_velocity(UI, velocity_type::U);
    set_velocity(VI, velocity_type::V);
}
Cell::Cell(double& PI, double& UI, double& VI, double& TI) :
_pressure(PI), _temperature(TI){
	set_velocity(UI, velocity_type::U);
	set_velocity(VI, velocity_type::V);

}
Cell::Cell(double& PI, double& UI, double& VI, cell_flag flag):
_pressure(PI), _flag(flag) {
	set_velocity(UI, velocity_type::U);
	set_velocity(VI, velocity_type::V);
}
;

// Pressure Get and Set
double& Cell::pressure() {
    return _pressure;
}

void Cell::set_pressure(double& value) {
    _pressure = value;
}

// Velocity Get and Set
double& Cell::velocity(velocity_type type) {
    return _velocity[static_cast<int>(type)];
};

void Cell::set_velocity(double& value, velocity_type type) {
    _velocity[static_cast<int>(type)] = value;
};

// borders Get and Set
bool& Cell::border(border_position position) {
    return _border[static_cast<int>(position)];
};

// Set border
void Cell::set_border(border_position position) {
    _border[static_cast<int>(position)] = true;
}
cell_flag& Cell::flag()
{
	return _flag;
}
void Cell::set_flag(cell_flag flag)
{
	_flag = flag;
}
;

// Get + Set temperature
double& Cell::temperature(){
    return _temperature;
}

void Cell::set_temperature(double& value){
    _temperature=value;
};

double& Cell::density(){
    return _density;
}

void Cell::set_density(double& value){
    _density=value;
};

Neighbourhood& Cell::get_neighbourhood(){
    return _neighbourhood;
}
void Cell::add_neighbour(unsigned int n){
    _neighbourhood.val= _neighbourhood.val | n;
}
