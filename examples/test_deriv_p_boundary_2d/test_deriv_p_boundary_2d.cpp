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
	Dimension1D dim_1(0.0,1.0,11);
	Dimension1D dim_2(0.0,1.0,21);

	// Make grid
	Grid grid({&dim_1,&dim_2});

	// Fill randomly
	IdxSet v(2);
	for (v[0]=1; v[0]<=11; v[0]++) {
		for (v[1]=1; v[1]<=21; v[1]++) {
			grid.get_grid_point_inside(v)->set_ordinate(fRand(0.0,5.0));
		};
	};

	// Write
	grid.write_to_file("test_deriv_p_boundary_2d.txt");

	// Point to evaluate at
	double* x = new double[2];
	x[0] = 0.023;
	x[1] = 0.993;

	// NOTE: derivative wrt p0x or px3 will throw an error, since it is out of the grid
	IdxSet local_idxs(2);
	double x_deriv;
	for (local_idxs[0]=1; local_idxs[0]<=3; local_idxs[0]++) {
		for (local_idxs[1]=0; local_idxs[1]<=2; local_idxs[1]++) {
			x_deriv = grid.get_deriv_wrt_pt_value(x,local_idxs);
			std::cout << "deriv @ " << x[0] << "," << x[1] << " wrt p" << local_idxs[0] << local_idxs[1] << " = " << x_deriv << std::endl;
		};
	};

	// Clean up
	delete[] x;
	x = nullptr;

	return 0;
};