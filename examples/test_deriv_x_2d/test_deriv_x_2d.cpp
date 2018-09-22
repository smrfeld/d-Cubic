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
	std::vector<int> v({0,0});
	for (int i=0; i<30; i++) {
		for (int j=0; j<14; j++) {
			v[0] = i;
			v[1] = j;
			grid.get_grid_point_ref(v).set_ordinate(fRand(-4.0,-2.0));
		};
	};

	// Write
	grid.write_to_file("test_deriv_x_2d.txt");

	// Point to evaluate at
	std::vector<double> abcissa({0.46,0.33});

	// Derivs
	double x_deriv = grid.get_deriv_wrt_x(abcissa,0);
	std::cout << "deriv @ " << abcissa[0] << "," << abcissa[1] << " wrt x = " << x_deriv << std::endl;
	double y_deriv = grid.get_deriv_wrt_x(abcissa,1);
	std::cout << "deriv @ " << abcissa[0] << "," << abcissa[1] << " wrt y = " << y_deriv << std::endl;

	return 0;
};