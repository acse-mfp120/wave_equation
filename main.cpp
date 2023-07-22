// local includes
#include "utils.h"
#include "SettingsConfig.h"

// system includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <mpi.h>
#define _USE_MATH_DEFINES
#include <math.h>

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	double start_time, end_time;
	start_time = MPI_Wtime();

	int id, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// create the config object
	Configs config;

	if (id == 0)
	{
		// have proc 0 read the file and populate its config object
		config.ParseConfigFile();
	}

	// now broadcast the config object to everyone
	config.BroadcastConfigs();
	double t{ 0.0 };
	double t_out{ config.configs.t_save_period_ };
	utils::create_data_folder();
	// Use MPI to pick an arrangement for the procs
	// dims will store the suggested number of procs to be arranged along each dimension.
	// Note: initialising our dims array to 0 ensures that MPI knows it has free reign to
	// choose an optimal way to divide up the domain. The answers are generally
	// good, but for prime numbers of processors, it becomes clear that for critical computations,
	// one might want to use an irregular decomposition instead of a cartesian one, or simply
	// restrict users to a number of procs and grid size.
	int dims[2] = { 0, 0 };
	if (MPI_Dims_create(p, 2, dims) != MPI_SUCCESS)
	{
		// bail out if this fails, since if it does, the rest of the logic breaks
		// conceptually
		throw std::runtime_error("Could not create valid processor grid arrangement");
	};
	// Ensure that we don't have too many procs compared to blocks per proc
	// imax and jmax are the total number of horiz and verticle blocks
	int temp_x_extents_bl = config.configs.x_max_ / config.configs.dx_;
	int temp_y_extents_bl = config.configs.y_max_ / config.configs.dy_;
	// initial blocks per proc calc. These figures will be augmented based for ghost cells
	int xblocks_per_proc{ temp_x_extents_bl / dims[0] };
	int yblocks_per_proc{ temp_y_extents_bl / dims[1] };
	// for now, we're just going to notifiy the user via a print statement that some of the
	// settings they chose will have to be altered in order to run the program effectively.
	// In the future, it would be nicer if these outputs could be written to a file, and not
	// std::out.
	if (id == 0)
	{
		// add the factor of 1000 in order to ensure that we can turn dx and dy into ints without rounding errors
		if (1000 * config.configs.x_max_ % (int)(1000 * config.configs.dx_) != 0 ||
			1000 * config.configs.y_max_ % (int)(1000 * config.configs.dy_) != 0)
		{
			// what we want is for the total grid size to be multiple of the dx and dy resolution
			std::cout << "WARNING: The blocks won't add up to your requested grid size!" << std::endl;
			/// TODO: turn this warning into a file that gets printed at the end if in debug mode.
		}
		if (xblocks_per_proc == 0 || yblocks_per_proc == 0)
		{
			std::cout << "Too many procs on the dancefloor" << std::endl;
			throw std::runtime_error("Too many procs for the given grid size.");
		}
		if (temp_x_extents_bl % dims[0] != 0 || temp_y_extents_bl % dims[1] != 0)
		{
			std::cout << "WARNING: The number of blocks assigned to each proc does not divide evenly. Grid will be smaller than you intended" << std::endl;
		}
	}
	// Work backwards and obtain the grid size that each proc thinks they're dealing with
	int actual_global_grid_bl[2] = { xblocks_per_proc * dims[0], yblocks_per_proc * dims[1] }; // global grid size in blocks
	int actual_blocks_per_proc_ng[2] = { xblocks_per_proc, yblocks_per_proc }; // actual blocks per grid excl. ghost cells
	double actual_res_units[2] = { (double)config.configs.x_max_ / actual_global_grid_bl[0] , (double)config.configs.y_max_ / actual_global_grid_bl[1] };

	// Now create new cartesian proc arrangement and a new communicator for this group
	// What this means is that we're relying on MPI to take the grid layout that it
	// suggested, and then create a custom communicator for that proc arrangement
	// which will include helpful meta information for domain decomp applications.
	// I.e. this communicator group will allow us to determine the neighbours of
	// each proc much easier than if we manually implemented it, especially for
	// periodic x or y dimensions.

	// specify if x and y boundaries are periodic (periodic = true/1)
	int periods[2] = { config.configs.x_boundary_, config.configs.y_boundary_ };
	int reorder = false; // allows MPI to assign arb ranks
	// now create the group and the communicator for the group
	MPI_Comm MPI_COMM_CART;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &MPI_COMM_CART);
	// and get the updated rank info and co-ordinates of this new proc arrangement
	int rank;
	int grid_placement[2];
	MPI_Comm_rank(MPI_COMM_CART, &rank);
	MPI_Cart_coords(MPI_COMM_CART, rank, 2, grid_placement);
	// Note: the grid placement will be used later for the initial disturb placement

	// create enum to easily retrieve the desired neighbour from the neighbour arr
	enum directions { L, R, U, D };
	std::vector<int> neighbours(4, 0);
	// note: direction 0 means moving along the x axis (left then right)
	MPI_Cart_shift(MPI_COMM_CART, 0, 1, &neighbours[L], &neighbours[R]);
	MPI_Cart_shift(MPI_COMM_CART, 1, 1, &neighbours[D], &neighbours[U]);

	// the offsets represent the global (x,y) location of the bottom
	// left position of the proc's grid
	double x_offset_units{ (double)grid_placement[0] * actual_blocks_per_proc_ng[0] * actual_res_units[0] };
	double y_offset_units{ (double)grid_placement[1] * actual_blocks_per_proc_ng[1] * actual_res_units[1] };
	// now figure out how many extra cols and rows we need for ghosts.
	// Note: could use smarter method to avoid repeated if statements
	// but this would require clever logic in order to avoid short
	// circuiting errors due to the logical && and || operators
	if (neighbours[L] != MPI_PROC_NULL)
	{
		// neighbour, thus require extra col for ghost cells
		xblocks_per_proc += 1;
	}
	if (neighbours[R] != MPI_PROC_NULL)
	{
		// same situation as above
		xblocks_per_proc += 1;
	}
	if (neighbours[U] != MPI_PROC_NULL)
	{
		// upper neighbour, thus require extra row for ghost cells
		yblocks_per_proc += 1;
	}
	if (neighbours[D] != MPI_PROC_NULL)
	{
		// same as above
		yblocks_per_proc += 1;
	}

	// now create the correctly sized grids for each proc
	int size_for_proc{ xblocks_per_proc * yblocks_per_proc };
	/// TODO: replace this with dynamically alloc raw pointers after debug complete!
	std::vector<double> grid(size_for_proc, 0);
	std::vector<double> old_grid(size_for_proc, 0);
	std::vector<double> new_grid(size_for_proc, 0);
	// dynamicall allocate raw array and init to 0
	//double* grid{ new double[size_for_proc] {} };
	//double* old_grid{ new double[size_for_proc] {} };
	//double* new_grid{ new double[size_for_proc] {} };

	// Store the particulars of the cart layout to file
	// for the post-processing step
	if (id == 0)
	{
		// post-processing script will make use of this
		// file
		utils::print_cart_layout_to_file(p, dims);
	}

	/// Now, we want to setup the initial state of the grid
	/// This could likely be speed up by means of placing a bounding box
	/// around each of the disturbances and then only searching within
	/// that proc instead of the entire proc grid limits.
	/// Implementing that however proved quite complicated, hence this
	/// simplified implementation for now.
	for (int k = 0; k < config.num_init_disturbances_; k++)
	{
		for (int i = 1; i < yblocks_per_proc - 1; i++)
		{
			for (int j = 1; j < xblocks_per_proc - 1; j++)
			{
				// check if the point falls within the circle. 
				// Recall that we're taking the co-ords of each 'block'/'pixel' as being the value of 
				// the point in the bottom left corner of the block/pixel, where that value is assigned 
				// according to a cartesian co-rdinate system with its origin at the bottom left of the
				// global grid
				int adj_x = (neighbours[L] != MPI_PROC_NULL) ? 1 : 0;
				int adj_y = (neighbours[D] != MPI_PROC_NULL) ? 1 : 0;

				double x = x_offset_units + (double(j) - adj_x) * actual_res_units[0];
				double y = y_offset_units + (double(i) - adj_y) * actual_res_units[1];

				double dist = sqrt(pow((x - config.init_disturbs[k].x), 2) + 
							       pow((y - config.init_disturbs[k].y), 2));

				if (dist <= config.init_disturbs[k].rad)
				{
					grid[i * xblocks_per_proc + j] = 5.0 * (cos(dist / config.init_disturbs[k].rad * M_PI) + 1.0);
					old_grid[i * xblocks_per_proc + j] = 5.0 * (cos(dist / config.init_disturbs[k].rad * M_PI) + 1.0);
				}
			}
		}
	}
	// have to wait for everyone to get to this point before we can start exchanging ghost cells
	MPI_Barrier(MPI_COMM_CART);

	// now do halo exchange to finish getting the initial grid setup
	utils::halo_and_borders(grid.data(), neighbours.data(), xblocks_per_proc, yblocks_per_proc, config.configs.boundary_, grid_placement, MPI_COMM_CART, rank);
	utils::halo_and_borders(old_grid.data(), neighbours.data(), xblocks_per_proc, yblocks_per_proc, config.configs.boundary_, grid_placement, MPI_COMM_CART, rank);
	utils::print_grid_to_file(id, 0, xblocks_per_proc, yblocks_per_proc, neighbours.data(), actual_blocks_per_proc_ng, grid.data());

	// now we can start with the main algorithm since everything is finally setup
	int count{ 1 };
	while (t < config.configs.t_max_)
	{
		utils::update_grid(new_grid.data(), grid.data(), old_grid.data(), xblocks_per_proc, yblocks_per_proc, config.configs);
		utils::halo_and_borders(new_grid.data(), neighbours.data(), xblocks_per_proc, yblocks_per_proc, config.configs.boundary_, grid_placement, MPI_COMM_CART, rank);
		// incrementing t
		t += config.configs.t_res_;
		// swop the grids
		old_grid.swap(new_grid);
		old_grid.swap(grid);
		// print the grid out at regular intervals
		if (t_out <= t)
		{
			utils::print_grid_to_file(id, count, xblocks_per_proc, yblocks_per_proc, neighbours.data(), actual_blocks_per_proc_ng, grid.data());
			count++;
			t_out += config.configs.t_save_period_;
		}
	}

	// not using raw pointers probably introduces a slow down for very large grid sizes
	//delete[] grid;
	//delete[] old_grid;
	//delete[] new_grid;

	end_time = MPI_Wtime();
	if (id == 0) {
		printf("Elapsed time is %fs\n", end_time - start_time);
	}
	MPI_Finalize();
}