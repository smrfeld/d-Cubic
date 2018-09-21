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

	grid.write_to_file("test_interp.txt");

	double x = grid.get_val({0.71});
	std::cout << "val @ 0.71 = " << x << std::endl;

	return 0;
};