#ifndef __SETTINGS_CONFIG_H__
#define __SETTINGS_CONFIG_H__
/*
* Author: Miguel Fidel Pereira
* Course: ACSE-6
* Assessment 2: MPI Assignment
* Date Created: 2021-04-05
*/
// Note: the comments are like this because I was going to use doxygen to
// produce nice looking documentation, but later scrapped that idea.
/**
* @file SettingsConfig.h
*
* @brief This file contains definitions required to parse and manage the settings throughout the program.
*/
// local includes
// system includes
#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

class Configs
{
// public functions and variables included in the MPI datatype
public:
	/// Configerations is a struct that will be used to create a custom
	/// MPI datatype which will then be used to pass the configs along
	/// to all of the other procs.
	// initialised with sensible variables to give users an
	// idea of the range of values they might use in the config.txt file
	struct Configerations
	{
		/// boundary code:
		/// 0 = Dirichlet
		/// 1 = Neumann
		/// 2 = Robin
		int boundary_{ 0 };
		/// Fixed (0) or Periodic (1) boundary for x and y
		int x_boundary_{ 0 };
		int y_boundary_{ 0 };
		/// Grid Size
		int x_max_{ 20 };
		int y_max_{ 20 };
		/// Grid Resolution
		double dx_{ 1 };
		double dy_{ 1 };
		/// Time Resolution of the simulation
		double t_res_{ 0.5 };
		/// Output Period
		double t_save_period_{ 2 };
		/// Speed of the wave
		double c_{ 1 };
		/// Max Runtime
		double t_max_{ 4 };
	};

	/// Custom struct used to create MPI datatype for
	/// communicating the initial disturbances to all procs.
	struct InitialDisturbances
	{
		double x{};
		double y{};
		double rad{};
	};

	/// Internal instance of the config struct
	/// Gets populated by parse_config_file()
	Configerations configs{};
	/// Internal instance to hold vector of the initial disturbances object
	InitialDisturbances init_disturbs[100];
	int num_init_disturbances_{};

// public functions and variables not included in the MPI datatype
public:
	/// Parses the config file found in the project's root directory.
	/// @returns exit code indicating success (0) or failure (1).
	int ParseConfigFile();

	/// Broadcast the Config object to all other processes
	void BroadcastConfigs();

private:
	/// Config Map used to assign numeric codes to each string definition
	/// or key used in the Config.txt file
	std::unordered_map<std::string, int> decode_key{
		{"dir", 0},
		{"neu", 1},
		{"robin", 2},
		{"fixed", 0},
		{"periodic", 1}
	};

	// Note: not making these static since the whole concept
	// of a Settings/Config class is that you have one per program
	// lifetime. I.e. one settings file which defines the entire
	// simulation, i.e. program execution.
	MPI_Datatype MPI_Sim_Configs;
	MPI_Datatype MPI_Init_Disturbs;

	/// Build and commit the derived MPI datatype
	void build_MPI_Config_Types();

	/// Splits a string into tokens using the given separation character
	/// Really just a util function for parsing the config files, but
	/// doesn't belong in utils namespace because it's only related to
	/// what goes on in class Configs
	std::vector<std::string> split_by_sep(std::string input, char sep);
};

#endif
