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
	Grid grid({&dim});

	// Fill randomly
	std::vector<int> v({0});
	for (int i=0; i<30; i++) {
		v[0] = i;
		grid.get_grid_point_ref(v).set_ordinate(fRand(-4.0,-2.0));
	};

	// Write
	grid.write_to_file("test_deriv_x_1d.txt");

	// Point to evaluate at
	std::vector<double> abcissa({0.46});

	// Derivs
	double x_deriv = grid.get_deriv_wrt_x(abcissa,0);
	std::cout << "deriv @ " << abcissa[0] << " wrt x = " << x_deriv << std::endl;

	return 0;
};