//
// Created by Moritz Gnisia on 05.04.20.
//

#ifndef CFDLAB_ENUMS_H
#define CFDLAB_ENUMS_H

#define B_N 1
#define B_S 2
#define B_W 4
#define B_E 8

enum class velocity_type {
    U,
    V
};

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT
};

enum class matrix_selection {
    ROW,
    COLUMN
};

enum class read_type {
    INT,
    DOUBLE
};

enum class simulation {
	MOVING_LID,
	PLANAR_SHEAR_FLOW,
	KARMAN_VORTEX_STREET,
	STEP,
	NATURAL_CONVECTION,
	FLUID_TRAP,
	RAYLEIGH_BENARD_CONVECTION,
	NO_SIM
};

enum class cell_flag : unsigned char {
	NO_SLIP,
	FREE_SLIP,
	OUTFLOW,
	INFLOW,
	FLUID,
    MULTISPECIES_0,
    MULTISPECIES_1,
    MULTISPECIES_2,
    MULTISPECIES_3,
    MULTISPECIES_4,
    SPECIES_CONVERTER_0,
    SPECIES_CONVERTER_1,
    SPECIES_CONVERTER_2,
    SPECIES_CONVERTER_3,
    SPECIES_CONVERTER_4
};

 /// B_E, B_W, B_S, B_N
struct Neighbourhood {
    unsigned int val : 4; 
};

struct WallTmps{
    double N=0;
    double S=0;
    double W=0;
    double E=0;
};

#endif //CFDLAB_ENUMS_H
