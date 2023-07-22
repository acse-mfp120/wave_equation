/*
* Author: Miguel Fidel Pereira
* Course: ACSE-6
* Assessment 2: MPI Assignment
* Date Created: 2021-04-05
*/

// File comment for doxygen.
/**
* @file tests.cpp
*
* @brief Contains unit tests which support and verify code.
*/

// local includes
#include "utils.h"
// system includes
#include <assert.h>
#include <vector>
#include <fstream>
#include <string>

using namespace std::chrono_literals;

void test_file_printing()
{
	// grid of size 4x3
	std::vector<double> grid{ 0,1,1,  0,1,0,  0,0.4,5,  1,1,0};
	// now print to file
	utils::print_grid_to_file(0, std::chrono::microseconds(5000), 4, 3, grid.data());
	// now read the resulting file and check that it prints to the correct location
	std::ifstream in;
	in.open("0_4_3_5000.txt");
	if (in.good()) {
		std::string line;
		std::vector<std::string> ans{ {"0 1 1"}, {"0 1 0"}, {"0 0.4 5"}, {"1 1 0"} };
		int counter{ 0 };
		while (std::getline(in, line)) {
			// this code will execute every time a newline is reached!
			assert(line == ans[counter]);
			counter++;
		};
	}
	else {
		// throw error
		assert(false);
	}
}


int tests()
{
	// test the file-writing
	test_file_printing();

	return 0;

}