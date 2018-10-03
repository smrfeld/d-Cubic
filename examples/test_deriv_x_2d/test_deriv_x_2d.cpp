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

	// Make dims
	Dimension1D dim_1(0.0,1.0,30);
	Dimension1D dim_2(0.0,1.0,14);

	// Make grid
	Grid grid({&dim_1,&dim_2});

	// Fill randomly
	IdxSet v(2);
	for (v[0]=1; v[0]<=30; v[0]++) {
		for (v[1]=1; v[1]<=14; v[1]++) {
			grid.get_grid_point_inside(v)->set_ordinate(fRand(-4.0,-2.0));
		};
	};

	// Write
	grid.write_to_file("test_deriv_x_2d.txt");

	// Point to evaluate at
	double* x = new double[2];
	x[0] = 0.46;
	x[1] = 0.33;

	// Derivs
	double x_deriv = grid.get_deriv_wrt_x(x,0);
	std::cout << "deriv @ " << x[0] << "," << x[1] << " wrt x = " << x_deriv << std::endl;
	double y_deriv = grid.get_deriv_wrt_x(x,1);
	std::cout << "deriv @ " << x[0] << "," << x[1] << " wrt y = " << y_deriv << std::endl;

	// Clean up
	delete[] x;
	x = nullptr;

	return 0;
};