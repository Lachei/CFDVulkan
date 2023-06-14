#define MAXSPECIES 20
#define GAS_CONSTANT 8.314
#define SOR_INT_TO_FLOAT(x) (x) * 1e-10f
#define SOR_FLOAT_TO_INT(x) uint((x) * 1e10f)

#define LOCAL_SIZE_X 32
#define LOCAL_SIZE_Y 32

struct SimCtx{
    float Re, UI, VI, PI, TI, UIn, VIn, WTN, WTE, WTS, WTW, GX, GY, t_end, xlenght, ylength;
    float dt;           //old dt
    float dx, dy;
    int imax, jmax;
    uint alpha;         //alpha has to be set to 0 after each iteration of simulation
    float omg, tau;  
    int itermax;
    float eps, dt_value, Pr, beta;
    float sor_eps, cur_sor_eps, prev_sor_eps;
    uint divider, lock, sor_counter;
    uint new_dt;        // new dt is interpreted as uint to allow atomic_min. This has to be set to t_end at the end of each iteration of simulation
    int patch_x, patch_y, patch_z, patch_x_origin, patch_y_origin, patch_z_origin;
    int model, o2index;          // used for fire simulations
    int amt_species;
};

const int NO_SLIP             = 0;
const int FREE_SLIP           = 1;
const int OUTFLOW             = 2;
const int INFLOW              = 3;
const int FLUID               = 4;
const int MULTISPECIES_0      = 5;
const int MULTISPECIES_1      = 6;
const int MULTISPECIES_2      = 7;
const int MULTISPECIES_3      = 8;
const int MULTISPECIES_4      = 9;
const int SPECIES_CONVERTER_0 = 10;
const int SPECIES_CONVERTER_1 = 11;
const int SPECIES_CONVERTER_2 = 12;
const int SPECIES_CONVERTER_3 = 13;
const int SPECIES_CONVERTER_4 = 14;

const int B_N = 1;
const int B_S = 2;
const int B_E = 4;
const int B_W = 8;

const ivec2 top          = ivec2(0,1);
const ivec2 right        = ivec2(1,0);
const ivec2 bot          = ivec2(0,-1);
const ivec2 left         = ivec2(-1,0);
const ivec2 top_right    = ivec2(1,1);
const ivec2 top_left     = ivec2(-1,1);
const ivec2 bot_right    = ivec2(1,-1);
const ivec2 bot_left     = ivec2(-1,-1);