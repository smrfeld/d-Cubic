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
	IdxSet v(1);
	for (v[0]=1; v[0]<=30; v[0]++) {
		grid.get_grid_point_inside(v)->set_ordinate(fRand(0.0,5.0));
	};

	// Write
	grid.write_to_file("test_interp_1d.txt");

	// Get val
	double* x = new double[1];
	x[0] = 0.71;
	double y = grid.get_val(x);
	std::cout << "val @ 0.71 = " << y << std::endl;

	// Clean up
	delete[] x;
	x = nullptr;

	return 0;
};