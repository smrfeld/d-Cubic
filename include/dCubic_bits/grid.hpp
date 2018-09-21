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

#ifndef IDX_SET_H
#define IDX_SET_H
#include "idx_set.hpp"
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
		IdxSet idxs_i;
		std::unordered_map<GridPtKey, const GridPt*, hash_gpk> in;
		std::vector<double> frac_abscissas;
	};

	/****************************************
	Neighborhood of points surrounding a point, 4 in each dim
	****************************************/

	struct Nbr4 {
		IdxSet idxs_i;
		std::unordered_map<GridPtKey, GridPtType, hash_gpk> types;
		std::unordered_map<GridPtKey, const GridPt*, hash_gpk> in;
		std::unordered_map<GridPtKey, const GridPtOut*, hash_gpk> out;
		std::vector<double> frac_abscissas;
	};

	/****************************************
	Ordinates associated with a Nbr2
	****************************************/

	struct P2 {
		std::unordered_map<GridPtKey, double, hash_gpk> p;

		P2(Nbr2 nbr2);
	};

	/****************************************
	Ordinates associated with a Nbr4
	****************************************/

	struct P4 {
		std::unordered_map<GridPtKey, double, hash_gpk> p;

		P4(Nbr4 nbr4);
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
		Nbr2 get_surrounding_2_grid_pts(std::vector<double> abscissas) const;
		Nbr4 get_surrounding_4_grid_pts(std::vector<double> abscissas) const;

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
		1D
		********************/

		double f1d_interpolate(double x_frac, double p0, double p1, double p2, double p3) const;
		double f1d_interpolate_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const;

		// p_idx = 0,1,2,3, depending on loc
		double f1d_deriv_pt_value(double x_frac, int p_idx) const;
		double f1d_deriv_pt_value_by_ref(const double &x_frac, int p_idx) const;

		double f1d_deriv_x(double x_frac, double p0, double p1, double p2, double p3) const;
		double f1d_deriv_x_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const;

	    /********************
		Apply M or P mappings
		********************/

	    IdxSet apply_m_mapping_by_ref(const IdxSet &idxs) const;
	    IdxSet apply_p_mapping_by_ref(const IdxSet &idxs) const;

		/********************
		Read/write grid
		********************/

		void read_from_file(std::string fname);
		void write_to_file(std::string fname) const;
	};

};