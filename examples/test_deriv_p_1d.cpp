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

	Dimension1D dim(0.0,1.0,30);

	Grid grid({dim});

	std::vector<int> v({0});
	for (int i=0; i<30; i++) {
		v[0] = i;
		grid.set_grid_point_ordinate(grid.get_grid_point(v), fRand(0.0,5.0));
	};

	grid.write_to_file("test_deriv_p_1d.txt");

	std::vector<double> abcissa({0.71});

	// Deriv WRT p0
	std::vector<int> local_idxs({0});
	double x_deriv = grid.get_deriv_wrt_pt_value(abcissa,local_idxs);
	std::cout << "deriv @ " << abcissa[0] << " wrt p" << local_idxs[0] << " = " << x_deriv << std::endl;

	// Deriv WRT p1
	local_idxs[0] = 1;
	x_deriv = grid.get_deriv_wrt_pt_value(abcissa,local_idxs);
	std::cout << "deriv @ " << abcissa[0] << " wrt p" << local_idxs[0] << " = " << x_deriv << std::endl;

	// Deriv WRT p2
	local_idxs[0] = 2;
	x_deriv = grid.get_deriv_wrt_pt_value(abcissa,local_idxs);
	std::cout << "deriv @ " << abcissa[0] << " wrt p" << local_idxs[0] << " = " << x_deriv << std::endl;

	// Deriv WRT 3
	local_idxs[0] = 3;
	x_deriv = grid.get_deriv_wrt_pt_value(abcissa,local_idxs);
	std::cout << "deriv @ " << abcissa[0] << " wrt p" << local_idxs[0] << " = " << x_deriv << std::endl;

	return 0;
};