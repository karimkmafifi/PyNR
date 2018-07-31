#ifndef INPUT_H
#define INPUT_H

#include "common.h"

struct Input
{
public:
	std::string receptor_path;
	std::string flex_path;
	std::string ligand_path;
	std::string unit_length;
	std::string center[3];
	std::string size[3];
	//std::string spacing;
	std::string out_dir;
	std::string out_name;
	std::string log_name;
	std::string cpu;
	std::string seed;
	std::string exhaustiveness;
	std::string num_modes;
	std::string energy_range;
	bool score_only;
	//bool default_grid;
};

struct InputFinal
{
public:
	std::string receptor_path;
	std::string flex_path;
	std::string ligand_path;
    float center[3];
    float size[3];
    //float spacing;

	std::string out_path;
	std::string log_path;

	int cpu;
	int seed;
	int exhaustiveness;
	int verbosity;
	int num_modes;

    float energy_range;

	bool flexb;
	bool logb;
	bool score_only;
	//bool default_grid;
};

#endif
