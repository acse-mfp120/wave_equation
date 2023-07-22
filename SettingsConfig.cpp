// local includes
#include "SettingsConfig.h"
#include "utils.h"
// system includes
#include <sstream>
#include <algorithm>
#include <fstream>

int Configs::ParseConfigFile()
{
	std::string line;
	std::vector<std::string> cleaned_file{};
	std::ifstream in;
	// open the file
	in.open("config.txt");
	if (in.good()) {
		while (std::getline(in, line)) {
			// if the line is a comment, ignore
			if (line[0] == '/') continue;
			// else add to the cleaned 'file'
			cleaned_file.push_back(line);
		}
		// now that we only have the relevant lines, store each of the settings
		configs.boundary_ = decode_key[split_by_sep(cleaned_file[0], ':')[1]];
		configs.x_boundary_ = decode_key[split_by_sep(cleaned_file[1], ':')[1]];
		configs.y_boundary_ = decode_key[split_by_sep(cleaned_file[2], ':')[1]];
		std::string temp;
		temp = split_by_sep(cleaned_file[3], ':')[1];
		configs.x_max_ = atoi(temp.c_str());
		temp = split_by_sep(cleaned_file[4], ':')[1];
		configs.y_max_ = atoi(temp.c_str());
		temp = split_by_sep(cleaned_file[5], ':')[1];
		configs.dx_ = atof(temp.c_str());
		temp = split_by_sep(cleaned_file[6], ':')[1];
		configs.dy_ = atof(temp.c_str());
		temp = split_by_sep(cleaned_file[7], ':')[1];
		configs.t_res_ = atof(temp.c_str());
		temp = split_by_sep(cleaned_file[8], ':')[1];
		configs.t_save_period_ = atof(temp.c_str());
		temp = split_by_sep(cleaned_file[9], ':')[1];
		configs.c_ = atof(temp.c_str());
		temp = split_by_sep(cleaned_file[10], ':')[1];
		configs.t_max_ = atof(temp.c_str());

		// now do initial conditions/disturbances
		int counter{};
		for (size_t i = 11; i < cleaned_file.size(); i++)
		{
			temp = split_by_sep(cleaned_file[i], ':')[1];
			// skip over the first and last chars (the brackets)
			temp = temp.substr(1, temp.size() - 2);
			std::vector<std::string> disturb = split_by_sep(temp, ',');
			InitialDisturbances dist;
			dist.x = atof(disturb[0].c_str());
			dist.y = atof(disturb[1].c_str());
			dist.rad = atof(disturb[2].c_str());
			init_disturbs[counter] = dist;
			counter++;
		}
		num_init_disturbances_ = counter;
	}
	else {
		return 1;
	}
	return 0;
}

void Configs::BroadcastConfigs()
{
	// constructs and committs the MPI data type
	build_MPI_Config_Types();
	// broadcast the data
	MPI_Bcast(&configs, 1, Configs::MPI_Sim_Configs, 0, MPI_COMM_WORLD);
	// free the derived type
	MPI_Type_free(&MPI_Sim_Configs);
	// broadcast the disturbances
	// first let everyone know how many disturbances there are
	MPI_Bcast(&num_init_disturbances_, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// now send them in the pre-allocated array
	MPI_Bcast(init_disturbs, 100, Configs::MPI_Init_Disturbs, 0, MPI_COMM_WORLD);
	MPI_Type_free(&MPI_Init_Disturbs);
}

void Configs::build_MPI_Config_Types()
{
	// define the list of types in the struct
	MPI_Datatype typelist[11] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT,
								MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
								MPI_DOUBLE, MPI_DOUBLE };
	// create the arrays to store the num_items, addresses and offsetts
	int block_lengths[11] = { 1,1,1,1,1,1,1,1,1,1,1 };
	MPI_Aint displacements[11];
	MPI_Aint addresses[11], add_start;

	Configerations temp;
	// Get the addresses
	// boundary settings
	MPI_Get_address(&temp.boundary_, &addresses[0]);
	MPI_Get_address(&temp.x_boundary_, &addresses[1]);
	MPI_Get_address(&temp.y_boundary_, &addresses[2]);
	// grid size and resolution
	MPI_Get_address(&temp.x_max_, &addresses[3]);
	MPI_Get_address(&temp.y_max_, &addresses[4]);
	MPI_Get_address(&temp.dx_, &addresses[5]);
	MPI_Get_address(&temp.dy_, &addresses[6]);
	// output save period
	MPI_Get_address(&temp.t_save_period_, &addresses[7]);
	// wave speed
	MPI_Get_address(&temp.c_, &addresses[8]);
	// max simulation length
	MPI_Get_address(&temp.t_max_, &addresses[9]);
	// time resolution
	MPI_Get_address(&temp.t_res_, &addresses[10]);

	// get the start address
	MPI_Get_address(&temp, &add_start);
	// calculate the offsets and populate the displacement array
	for (int i = 0; i < 11; i++)
	{
		displacements[i] = addresses[i] - add_start;
	}
	// Create the data type and commit
	MPI_Type_create_struct(11, block_lengths, displacements, typelist, &MPI_Sim_Configs);
	MPI_Type_commit(&MPI_Sim_Configs);

	// now create the initial disturbance type
	// Note that we can have multiple initial disturbances!
	MPI_Datatype typelist2[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
	int lengths[3] = { 1,1,1 };
	MPI_Aint displacements2[3], addresses2[3], add_start2;
	
	InitialDisturbances temp2;
	MPI_Get_address(&temp2.x, &addresses2[0]);
	MPI_Get_address(&temp2.y, &addresses2[1]);
	MPI_Get_address(&temp2.rad, &addresses2[2]);
	MPI_Get_address(&temp2, &add_start2);
	for (int i = 0; i < 3; i++)
	{
		displacements2[i] = addresses2[i] - add_start2;
	}
	// Create the data type and commit
	MPI_Type_create_struct(3, lengths, displacements2, typelist2, &MPI_Init_Disturbs);
	MPI_Type_commit(&MPI_Init_Disturbs);
}

std::vector<std::string> Configs::split_by_sep(std::string input, char sep)
{
	std::stringstream ss{ input };
	std::vector<std::string> vstrTokens{};
	std::string token;
	while (std::getline(ss, token, sep)) {
		if (!token.empty()) {
			// remove any whitespace before we add it to the vector
			auto start = std::remove_if(std::begin(token), std::end(token), [=](unsigned char c) {return c == ' '; });
			token.erase(start, std::end(token));
			// add cleaned token
			vstrTokens.push_back(token);
		}
	}
	return vstrTokens;
}
