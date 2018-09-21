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

	Dimension1D dim1(0.0,1.0,22);
	Dimension1D dim2(0.0,1.0,8);

	Grid grid({dim1,dim2});

	std::vector<int> v({0,0});
	for (int i=0; i<22; i++) {
		for (int j=0; j<8; j++) {
			v[0] = i;
			v[1] = j;
			grid.set_grid_point_ordinate(grid.get_grid_point(v), fRand(0.0,5.0));
		};
	};

	grid.write_to_file("test_interp_2d.txt");

	double x = grid.get_val({0.71,0.33});
	std::cout << "Val @ 0.71, 0.33 = " << x << std::endl;

	return 0;
};