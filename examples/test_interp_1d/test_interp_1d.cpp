#include <dcubic>

#include <iostream>

#include <random>
#include <ctime>

using namespace std;
using namespace dcu;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main() {

	// Seed
	srand( time( NULL ) );

	// Dim
	Dimension1D dim(0.0,1.0,30);

	// Grid
	Grid grid({&dim});

	// Fill randomly
	std::vector<int> v({0});
	for (int i=0; i<30; i++) {
		v[0] = i;
		grid.get_grid_point_ref(v).set_ordinate(fRand(0.0,5.0));
	};

	// Write
	grid.write_to_file("test_interp_1d.txt");

	// Get val
	double x = grid.get_val({0.71});
	std::cout << "val @ 0.71 = " << x << std::endl;

	return 0;
};