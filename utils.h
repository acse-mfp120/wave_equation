#ifndef __utils__
#define __utils__
/*
* Author: Miguel Fidel Pereira
* Course: ACSE-6
* Assessment 2: MPI Assignment
* Date Created: 2021-04-05
*/

// File comment for doxygen.
/**
* @file utils.h
* 
* @brief Utils contains useful routines required by the main executable, but which are not related to solving the wave equation.
*/

// local includes
#include "SettingsConfig.h"

// system includes
#include <chrono>
#include <vector>
#include <string>

namespace utils {
	/// Create output folder to store the printed output files
	bool create_data_folder();

	/// Prints grid to file
	/// @param proc			- Indicates the number of the process calling the function
	/// @param count		- ID that gets incremented after each output file written by proc
	/// @param x_dim		- Size of the x-domain incl. ghost cells
	/// @param y_dim		- Size of the y-domain incl. ghost cells
	/// @param neighbours   - Vector declaring ID of cells L,R,U,D neighbours
	/// @param file_size	- Array indicating the number of xblocks (cols) and yblocks (rows) to be printed
	/// @param grid			- Pointer to the contents of the grid stored as flattened array
	/// @returns exit code indicating success (0) or failure (1).
	int print_grid_to_file(int proc, int count, int x_dim, int y_dim, int* neighbours, int* file_size, const double* grid);

	/// Print cartesian proc arrangment details to a special config file (not where the other settings
	/// are kept) so that the post processing script can know the global arrangment of the output files
	/// @param num_procs	- The total number of procs
	/// @param grid_layout	- Two member array describing the number of procs arranged along each dimension
	int print_cart_layout_to_file(int num_procs, int* dims);

	/// Does halo exchange and border updates according to given border conditions
	/// @param grid			- Pointer to main grid stored by proc
	/// @param neighbours	- Array holding the name of the neighbours of the proc. Order is: Left, Right, Up, Down.
	/// @param xblocks		- The total number of columns assigned to this proc (includes ghost cells etc.)
	/// @param yblocks		- Same as above, but for rows.
	/// @param boundary_cond- Code to indicate the boundary condition. 0 = Dirichlet, 1 = Neumann.
	/// @param grid_placement - Position of proc within cartesian grid.
	/// @param comm			- The MPI communicator for the Cartesian grid arrangement that was used to determine the neighbours param.
	/// @param id			- The ID of the calling proc.
	void halo_and_borders(double* grid, int* neighbours, unsigned int xblocks, unsigned int yblocks, int boundary_cond, int *grid_placement, MPI_Comm &comm, unsigned int id);

	/// Do an update of the proc's grid
	/// @param new_grid			- The 'new' grid about to be populated with the latest values
	/// @param grid				- The current grid
	/// @param old_grid			- The previous grid
	/// @param xblocks			- The number of columns in the entire grid (ghost cells included)
	/// @param yblocks			- The number of rows in the entire grid (ghost cells included)
	/// @param config			- The configuration object
	void update_grid(double* new_grid, double* grid, double* old_grid, int xblocks, int yblocks, Configs::Configerations config);
}



#endif // !__utils__
