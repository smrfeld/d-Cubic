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

	// Dims
	Dimension1D dim1(0.0,1.0,22);
	Dimension1D dim2(0.0,1.0,8);

	// Grid
	Grid grid({dim1,dim2});

	// Fill randomly
	IdxSet v(2);
	for (v[0]=1; v[0]<=22; v[0]++) {
		for (v[1]=1; v[1]<=8; v[1]++) {
			grid.get_grid_point_inside(v)->set_ordinate(fRand(0.0,5.0));
		};
	};

	// Write
	grid.write_to_file("test_interp_2d.txt");

	// Get val
	double* x = new double[2];
	x[0] = 0.71;
	x[1] = 0.33;
	double y = grid.get_val(x);
	std::cout << "Val @ 0.71, 0.33 = " << y << std::endl;

	// Clean up
	delete[] x;
	x = nullptr;

	return 0;
};