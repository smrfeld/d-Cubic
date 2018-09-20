#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forwards
	class Dimension1D;
	class GridPt;
	class GridPtOut;
	class GridPtKey;
	class IdxSet;
	struct Nbr4;
	struct Nbr2;

	/****************************************
	Grid
	****************************************/

	class Grid {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		Grid(std::vector<Dimension1D> dims);
		Grid(const Grid& other);
		Grid(Grid&& other);
		Grid& operator=(const Grid &other);
		Grid& operator=(Grid &&other);
		~Grid();

		/********************
		Get dims
		********************/

		int get_no_dims() const;
		const std::vector<Dimension1D>& get_dims() const;

		/********************
		Get grid points
		********************/

		const std::map<GridPtKey, GridPt*>& get_grid_points() const;
		const GridPt* get_grid_point(std::vector<int> grid_idxs) const;
		const GridPt* get_grid_point(IdxSet idx_set) const;
		const GridPt* get_grid_point(GridPtKey key) const;

		const std::map<GridPtKey, GridPtOut*>& get_grid_points_outside() const;
		const GridPtOut* get_grid_point_outside(std::vector<int> grid_idxs) const;
		const GridPtOut* get_grid_point_outside(IdxSet idx_set) const;
		const GridPtOut* get_grid_point_outside(GridPtKey key) const;

		/********************
		Set grid point values
		********************/

		void set_grid_point_ordinate(const GridPt* grid_pt, double val);

		/********************
		Get grid points surrounding a point
		********************/

		// Second arg = fractions between 0,1 in each dim
		std::pair<Nbr2,std::vector<double>> get_surrounding_2_grid_pts(std::vector<double> abscissas) const;
		std::pair<Nbr4,std::vector<double>> get_surrounding_4_grid_pts(std::vector<double> abscissas) const;

		/********************
		Get a point by interpolating
		********************/

		double get_val(std::vector<double> abscissas) const;

		/********************
		Get derivative
		********************/

		double get_deriv_wrt_pt_value(std::vector<double> abscissas, std::vector<int> grid_idxs);
		double get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet idx_set);
		double get_deriv_wrt_pt_value(std::vector<double> abscissas, GridPtKey grid_pt_key);

		double get_deriv_wrt_x(std::vector<double> abscissas, int k);

		/********************
		1D funcs
		********************/

		double interpolate_1d(double x, double p0, double p1, double p2, double p3) const;
		double interpolate_1d_by_ref(const double &x, const double &p0, const double &p1, const double &p2, const double &p3) const;

		/********************
		Read/write grid
		********************/

		void read_from_file(std::string fname);
		void write_to_file(std::string fname) const;
	};

};