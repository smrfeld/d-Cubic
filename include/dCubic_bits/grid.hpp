#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef UNORDERED_MAP_H
#define UNORDERED_MAP_H
#include <unordered_map>
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
	enum class GridPtType: unsigned int;
	class IdxSet;
	enum class LocInDim: unsigned int;

	/****************************************
	Hash for an unordered map
	****************************************/

	// Hash
	struct hash_gpk {
	    size_t operator() ( const GridPtKey &grid_pt_key ) const;
	};

	/****************************************
	Neighborhood of points surrounding a point, 2 in each dim
	****************************************/

	struct Nbr2 {
		std::unordered_map<GridPtKey, const GridPt*, hash_gpk> in;
	};

	/****************************************
	Neighborhood of points surrounding a point, 4 in each dim
	****************************************/

	struct Nbr4 {
		std::unordered_map<GridPtKey, GridPtType, hash_gpk> types;
		std::unordered_map<GridPtKey, const GridPt*, hash_gpk> in;
		std::unordered_map<GridPtKey, const GridPtOut*, hash_gpk> out;
	};

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

		const std::unordered_map<GridPtKey, GridPt*, hash_gpk>& get_grid_points() const;
		const GridPt* get_grid_point(std::vector<int> grid_idxs) const;
		const GridPt* get_grid_point(IdxSet idx_set) const;
		const GridPt* get_grid_point(GridPtKey key) const;

		const std::unordered_map<GridPtKey, GridPtOut*, hash_gpk>& get_grid_points_outside() const;
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
		Get derivative wrt a point value p
		********************/

		// Here: grid_idxs/idx_set/grid_pt_key are 0,1,2,3 each
		// i.e. refer to p! not global idxs
		double get_deriv_wrt_pt_value(std::vector<double> abscissas, std::vector<int> local_grid_idxs);
		double get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet local_idx_set);
		double get_deriv_wrt_pt_value(std::vector<double> abscissas, GridPtKey local_grid_pt_key);

		/********************
		Get derivative wrt x
		********************/

		double get_deriv_wrt_x(std::vector<double> abscissas, int k);

		/********************
		1D funcs
		********************/

		double interpolate_1d(double x_frac, double p0, double p1, double p2, double p3) const;
		double interpolate_1d_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const;

		double get_deriv_wrt_p_1d(double x_frac, int p, LocInDim loc);
		double get_deriv_wrt_p_1d_by_ref(const double &x_frac, int p, LocInDim loc) const;
		
		double get_deriv_wrt_x_1d(double x_frac, double p0, double p1, double p2, double p3);
		double get_deriv_wrt_x_1d_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const;

		/********************
		Read/write grid
		********************/

		void read_from_file(std::string fname);
		void write_to_file(std::string fname) const;
	};

};