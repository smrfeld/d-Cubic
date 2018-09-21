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
	std::vector<int> v({0,0});
	for (int i=0; i<22; i++) {
		for (int j=0; j<8; j++) {
			v[0] = i;
			v[1] = j;
			grid.get_grid_point_ref(v).set_ordinate(fRand(0.0,5.0));
		};
	};

	// Write
	grid.write_to_file("test_interp_2d.txt");

	// Get val
	double x = grid.get_val({0.71,0.33});
	std::cout << "Val @ 0.71, 0.33 = " << x << std::endl;

	return 0;
};