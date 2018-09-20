#include <dCubic>

#include <iostream>

#include <random>

using namespace std;
using namespace dcu;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main() {

	// Make dims
	Dimension1D dim_1(0.0,1.0,30);
	Dimension1D dim_2(0.0,1.0,17);

	// Make grid
	Grid grid({dim_1,dim_2});

	// Fill randomly
	std::vector<int> v({0,0});
	for (int i=0; i<30; i++) {
		for (int j=0; j<17; j++) {
			v[0] = i;
			v[1] = j;
 			grid.set_grid_point_ordinate(grid.get_grid_point(v), fRand(0.0,5.0));
		};
	};

	// Write to file
	grid.write_to_file("test_deriv_p_2d.txt");

	// Point to evaluate at
	std::vector<double> abcissa({0.71,0.23});

	// Derivs
	std::vector<int> local_idxs({0,0});
	double x_deriv;
	for (int i=0; i<=3; i++) {
		for (int j=0; j<=3; j++) {
			local_idxs[0] = i;
			local_idxs[1] = j;
			x_deriv = grid.get_deriv_wrt_pt_value(abcissa,local_idxs);
			std::cout << "deriv @ " << abcissa[0] << "," << abcissa[1] << " wrt p" << local_idxs[0] << local_idxs[1] << " = " << x_deriv << std::endl;
		};
	};

	return 0;
};