// local includes
#include "utils.h"
// system includes
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#define _USE_MATH_DEFINES
#include <math.h>
//#include <experimental/filesystem>

#ifdef _WIN32
#include <Windows.h>
#endif // _WIN32

bool utils::create_data_folder()
{
#ifdef _WIN32
    std::string s{ "data" };
    std::wstring folder = std::wstring(s.begin(), s.end());
    return CreateDirectory(folder.c_str(), NULL);
#else
    // note: this code is problematic for anyone on gcc -v < 8.xxx
    // Hence commenting this out to avoid any trouble when examiners
    // run this code. In order to account for this however, you
    // need to have a data folder already created.
    std::filesystem::current_path(std::filesystem::temp_directory_path());
    return std::filesystem::create_directories("data");
#endif // _WIN32
}

// The format for the filename is as follows: 'proc_xdim_ydim_timestamp'
int utils::print_grid_to_file(int proc, int count, int xblocks, int yblocks, int *neighbours, int* im_size, const double* grid)
{
    // if neighbour below, then starting point = 0 + xblocks_per_proc
    // if neighbour to the left, then starting point is offset by 1 to account for the ghost cell, i.e. +1
    // if neighbour on the right, then stop one before xblock_per_proc
    // if neighbour above, then stop one before yblocks_per_proc

    int start_row = (neighbours[3] != MPI_PROC_NULL) ? 1 : 0;
    int end_row = (neighbours[2] != MPI_PROC_NULL) ? yblocks - 1 : yblocks;
    int start_col = (neighbours[0] != MPI_PROC_NULL) ? 1 : 0;
    int end_col = (neighbours[1] != MPI_PROC_NULL) ? xblocks - 1 : xblocks;

    // make an object of the input file stream class
    std::ofstream out;
    // create the filename
    std::string filename = std::to_string(proc) + "_" + std::to_string(im_size[0]) + "_" + std::to_string(im_size[1]) + "_" + std::to_string(count) + ".dat";
    // explicitely create the file new each time, i.e. create new file and overwrite.
#ifdef _WIN32
    // yes, this does work in Windows.
    filename = "data//" + filename;
#else
    // and this too.
    filename = "data/" + filename;
#endif // _WIN32

    out.open(filename, std::ofstream::out);
    if (out.good()) {
        for (int i = start_row; i < end_row; i++)
        {
            for (int j = start_col; j < end_col; j++)
            {
                out << grid[i * xblocks + j];
                if (j < end_col - 1) out << " ";
            }
            out << "\n";
        }
    }
    else {
        return 1;
    }
    // no need to explicitely close file since ofstream object is RAII
    return 0;
}

int utils::print_cart_layout_to_file(int num_procs, int* dims)
{
    // make an object of the input file stream class
    std::ofstream out;
    // create the filename
    std::string filename{ "cart_layout.txt" };
    // explicitely create the file new each time, i.e. create new file and overwrite.
    out.open(filename, std::ofstream::out);
    if (out.good()) {
        out << "// Output Arrangment for Cartesian Grid." << std::endl;
        out << num_procs << std::endl;
        out << dims[0] << "," << dims[1] << std::endl;
    }
    else {
        return 1;
    }
    // no need to explicitely close file since ofstream object is RAII
    return 0;
}

void utils::halo_and_borders(double* grid, int* neighbours, unsigned int xblocks, unsigned int yblocks, int boundary_cond, int *grid_placement, MPI_Comm &comm, unsigned int id)
{
    MPI_Request* request = new MPI_Request[4 * 2];
    unsigned int counter{};

    /// Create the vectors containing all of the required starting points
    /// Order of the arrays is same as neighbours vector: L, R, U, D
    // Start and Receive indeces
    std::vector<unsigned int> send_starting_pts{ xblocks + 1, xblocks - 2 + xblocks, xblocks * (yblocks - 2) + 1, xblocks + 1 };
    std::vector<unsigned int> rec_starting_pts{ xblocks, 2 * xblocks - 1, xblocks * (yblocks - 1) + 1, 1 };
    // Column and Border offsets for border start position
    unsigned int col_border_offset = (grid_placement[1] == 0) ? 0 : 1;
    unsigned int row_border_offset = (grid_placement[0] == 0) ? 0 : 1;
    // Border start positions (index)
    std::vector<unsigned int> border_starting_pts{ col_border_offset * xblocks, 
                                                   xblocks - 1 + col_border_offset * xblocks,
                                                   xblocks * (yblocks - 1) + row_border_offset,
                                                   row_border_offset };
    // Number of border blocks to update each time
    std::vector<unsigned int> num_border_blocks{ yblocks - 2, yblocks - 2, xblocks - 2, xblocks - 2 };
    // Offset to account for placement of proc within grid for border calculation
    std::vector<int> border_block_offsets{ 1, -1, -1 * (int)xblocks, (int)xblocks };
    // Strides and the number of values being sent for both ghost cells and borders
    std::vector<unsigned int> strides{ xblocks, xblocks, 1, 1 };
    std::vector<unsigned int> num_vals{ yblocks - 2, yblocks - 2, xblocks - 2, xblocks - 2 };

    for (int i = 0; i < 4; i++)
    {
        if (neighbours[i] != MPI_PROC_NULL)
        {
            MPI_Datatype MPI_type;
            MPI_Type_vector(num_vals[i], 1, strides[i], MPI_DOUBLE, &MPI_type);
            MPI_Type_commit(&MPI_type);
            // now send and receive into the data type
            MPI_Isend(&grid[send_starting_pts[i]], 1, MPI_type, neighbours[i], 0, comm, &request[counter * 2]);
            MPI_Irecv(&grid[rec_starting_pts[i]], 1, MPI_type, neighbours[i], 0, comm, &request[counter * 2 + 1]);
            counter++;
        }
        else
        {
            // select the value to update the border with depending on the boundary condition
            double val = (boundary_cond == 0) ? 0 : grid[border_starting_pts[i] + border_block_offsets[i]];
            // update border values
            for (int k = 0; k < num_border_blocks[i]; k++)
            {
                grid[(int)border_starting_pts[i] + (int)strides[i] * k] = val;
            }
        }

    }
    // wait for everything to finish
    MPI_Waitall(counter * 2, request, MPI_STATUSES_IGNORE);
    delete[] request;
}

void utils::update_grid(double* new_grid, double* grid, double* old_grid, int xblocks, int yblocks, Configs::Configerations config)
{
    // extract the settings
    double dt = config.t_res_;
    double c = config.c_;
    double dx = config.dx_;
    double dy = config.dy_;
    // do the updates
    for (int i = 1; i < xblocks - 1; i++)
    {
	    for (int j = 1; j < yblocks - 1; j++)
	    {
		    // centre
		    double centre = grid[i * xblocks + j];
		    // left
		    double left = grid[i * xblocks + j - 1];
		    // right
		    double right = grid[i * xblocks + j + 1];
		    // top
		    double top = grid[(i + 1) * xblocks + j];
		    // bottom
		    double bottom = grid[(i -1) * xblocks + j];
            // update the cell
		    new_grid[i * xblocks + j] = pow(dt * c, 2.0) * ((right - 2.0 * centre + left) / pow(dx, 2.0) + 
											    (top - 2.0 * centre + bottom) / pow(dy, 2.0)) + 
											    2.0 * centre - old_grid[i * xblocks + j];
	    }
    }
}
