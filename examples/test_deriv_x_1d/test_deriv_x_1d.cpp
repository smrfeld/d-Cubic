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

	// Make dim
	Dimension1D dim(0.0,1.0,30);

	// Make grid
	Grid grid({dim});

	// Fill randomly
	IdxSet v(1);
	for (v[0]=1; v[0]<=30; v[0]++) {
		grid.get_grid_point_inside(v)->set_ordinate(fRand(-4.0,-2.0));
	};

	// Write
	grid.write_to_file("test_deriv_x_1d.txt");

	// Point to evaluate at
	double* x = new double[1];
	x[0] = 0.46;

	// Derivs
	double x_deriv = grid.get_deriv_wrt_x(x,0);
	std::cout << "deriv @ " << x[0] << " wrt x = " << x_deriv << std::endl;

	// Clean up
	delete[] x;
	x = nullptr;

	return 0;
};