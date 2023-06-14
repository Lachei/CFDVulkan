#include "helper.hpp"
#include "init.hpp"
#include "iostream"
#include <fstream>
#include <algorithm>

int TEST_CASE=0;

std::ifstream& read_double(std::ifstream& f, double& d) {
	std::string tmp;
	f >> tmp;
	std::for_each(tmp.begin(), tmp.end(), [](unsigned char c) {return std::tolower(c); });
	if (tmp == "nan") {
		d = std::nan("");
		return f;
	}
	if (tmp == "-nan") {
		d = -std::nan("");
		return f;
	}
	if (tmp == "inf") {
		d = std::numeric_limits<double>::infinity();
		return f;
	}
	if (tmp == "-inf") {
		d = -std::numeric_limits<double>::infinity();
		return f;
	}
	try {
		d = std::stod(tmp);
		return f;
	}
	catch (...) {
		f.setstate(f.failbit);
		return f;
	}
}
void set_test_case(std::string& problem){
    if(problem == "death_valley")TEST_CASE=1;
    else if(problem == "indoor_fire")TEST_CASE=2;
}

int read_parameters( std::string&       szFileName,       /* name of the file */
                     SimulationContext& sim_ctx)          /* simulation context to be filled*/
{
	bool default_species=true;
	sim_ctx.gpu_index = -1;		/* by default cpu computing is used*/
	std::ifstream file(szFileName);
	if (!file.is_open()) return -1;
	std::string var;
	while (!file.eof() && file.good()) {
		file >> var;
		if (var[0] == '#') {     /* comment*/
			file.ignore(1024, '\n');
		}
		else {
			if (var == "xlength")              file >> sim_ctx.xlength;
			if (var == "ylength")              file >> sim_ctx.ylength;
			if (var == "Re")                   file >> sim_ctx.Re;
			if (var == "t_end")                file >> sim_ctx.t_end;
			if (var == "dt")                   file >> sim_ctx.dt;
			if (var == "omg")                  file >> sim_ctx.omg;
			if (var == "eps")                  file >> sim_ctx.eps;
			if (var == "tau")                  file >> sim_ctx.tau;
			if (var == "alpha")                file >> sim_ctx.alpha;
			if (var == "beta")                 file >> sim_ctx.beta;
			if (var == "Pr")                   file >> sim_ctx.Pr;
			if (var == "dt_value")             file >> sim_ctx.dt_value;
			if (var == "UI")                   file >> sim_ctx.UI;
			if (var == "VI")                   file >> sim_ctx.VI;
			if (var == "GX")                   file >> sim_ctx.GX;
			if (var == "GY")                   file >> sim_ctx.GY;
			if (var == "PI")                   file >> sim_ctx.PI;
			if (var == "TI")                   file >> sim_ctx.TI;
			if (var == "itermax")              file >> sim_ctx.itermax;
			if (var == "imax")                 file >> sim_ctx.imax;
			if (var == "jmax")                 file >> sim_ctx.jmax;
			if (var == "iproc")                file >> sim_ctx.iproc;
			if (var == "jproc")                file >> sim_ctx.jproc;
			if (var == "problem")              file >> sim_ctx.problem;
			if (var == "geometry")             file >> sim_ctx.geometry;
			if (var == "WTI")                  { read_double(file, sim_ctx.WTI.S); read_double(file, sim_ctx.WTI.W); read_double(file, sim_ctx.WTI.N); read_double(file, sim_ctx.WTI.E); }
			if (var == "UIn")                  file >> sim_ctx.UIn;
			if (var == "VIn")                  file >> sim_ctx.VIn;
            if (var == "Species")              { default_species=false; Species s{}; file >> s.name; file >> s.molar_mass; file >> s.initial_concentration; sim_ctx.species.push_back(s); }
            if (var == "MultiSpeciesInserter") { MultiInserter m{}; m.concentrations.resize(sim_ctx.species.size()); 
                                                    for (int i = 0; i < sim_ctx.species.size(); ++i) file >> m.concentrations[i]; 
                                                    file >> m.temp; file >> m.velocities[0]; file >> m.velocities[1];
                                                    sim_ctx.multi_inserters.push_back(m); }
			if (var == "SpeciesConverter")		{SpeciesCoverter sc{}; sc.conversion_ratio.resize(sim_ctx.species.size()); sc.reduction_per_qm_sec.resize(sim_ctx.species.size());
													for (int i = 0; i < sim_ctx.species.size(); ++i)file >> sc.reduction_per_qm_sec[i];
													for (int i = 0; i < sim_ctx.species.size(); ++i){
														sc.conversion_ratio[i].resize(sim_ctx.species.size());
														for (int j = 0; j < sim_ctx.species.size(); ++j)file >> sc.conversion_ratio[i][j];
													}
													read_double(file, sc.temp);
													sim_ctx.species_converters.push_back(sc); }
			if (var == "GpuIndex")				file >> sim_ctx.gpu_index;
		}
	}

	set_test_case(sim_ctx.problem);
	if(default_species){Species s{"Default_Species",18.01528,1.0}; sim_ctx.species.push_back(s);}
	
	sim_ctx.dx = sim_ctx.xlength / (double)(sim_ctx.imax);
	sim_ctx.dy = sim_ctx.ylength / (double)(sim_ctx.jmax);

	if (!file.good() && !file.eof()) return -1;

	return 1;
}
