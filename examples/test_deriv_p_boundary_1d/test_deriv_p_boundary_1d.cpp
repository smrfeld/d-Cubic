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
	Dimension1D dim(0.0,1.0,11);

	// Make grid
	Grid grid({&dim});

	// Fill randomly
	std::vector<int> v({0});
	for (int i=0; i<11; i++) {
		v[0] = i;
		grid.get_grid_point_ref(v).set_ordinate(fRand(0.0,5.0));
	};

	// Write
	grid.write_to_file("test_deriv_p_boundary_1d.txt");

	// Point to evaluate at
	std::vector<double> abcissa({0.023});

	// NOTE: derivative wrt p0 will throw an error, since it is out of the grid
	IdxSet4 local_idxs({0});
	double x_deriv;
	for (int i=1; i<=3; i++) {
		local_idxs[0] = i;
		x_deriv = grid.get_deriv_wrt_pt_value(abcissa,local_idxs);
		std::cout << "deriv @ " << abcissa[0] << " wrt p" << local_idxs[0] << " = " << x_deriv << std::endl;
	};

	return 0;
};