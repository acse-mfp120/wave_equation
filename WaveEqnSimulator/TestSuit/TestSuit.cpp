// necessary unit test includes
#include "pch.h"
#include "CppUnitTest.h"
// local includes
#include "../../utils.cpp"
// system includes
#include <string>
#include <fstream>
#include <vector>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace TestSuit
{
	TEST_CLASS(TestSuit)
	{
	public:
		
		TEST_METHOD(Test_File_Printing)
		{
			// grid of size 4x3
			std::vector<double> grid{ 0,1,1,  0,1,0,  0,0.4,5,  1,1,0 };
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
					Assert::AreEqual(ans[counter], line);
					counter++;
				};
			}
			else {
				// throw error
				Assert::Fail();
			}
		}
		//TEST_METHOD(Test_Config_Parsing)
		//{
		//	Assert::AreEqual(0, utils::parse_config_file());
		//}
	};
}
