#include <string>
#include <vector>
#include <map>

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
	Neighborhood of points surrounding a point, 2 in each dim
	****************************************/

	struct Nbr2 {
		IdxSet idxs_i;
		std::map<IdxSet2, const GridPt*> in;
		std::vector<double> frac_abscissas;

		Nbr2(IdxSet idxs_i);
	};

	/****************************************
	Neighborhood of points surrounding a point, 4 in each dim
	****************************************/

	struct Nbr4 {
		IdxSet idxs_i;
		std::map<IdxSet4, GridPtType> types;
		std::map<IdxSet4, const GridPt*> in;
		std::map<IdxSet4, const GridPtOut*> out;
		std::vector<double> frac_abscissas;

		Nbr4(IdxSet idxs_i);
	};

	/****************************************
	Ordinates associated with a Nbr2
	****************************************/

	struct P2 {
		std::map<IdxSet2, double> p;

		P2(Nbr2 nbr2);
	};

	/****************************************
	Ordinates associated with a Nbr4
	****************************************/

	struct P4 {
		std::map<IdxSet4, double> p;

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
		Print
		********************/

		void print_grid_pts_inside() const;

		/********************
		Get dims
		********************/

		int get_no_dims() const;
		const std::vector<Dimension1D>& get_dims() const;

		/********************
		Get grid points
		********************/

		/*****
		Const ptrs
		*****/

		const GridPt* get_grid_point(std::vector<int> grid_idxs) const;
		const GridPt* get_grid_point(IdxSet idx_set) const;
		const GridPt* get_grid_point(GridPtKey key) const;

		const GridPtOut* get_grid_point_outside(std::vector<int> grid_idxs) const;
		const GridPtOut* get_grid_point_outside(IdxSet idx_set) const;
		const GridPtOut* get_grid_point_outside(GridPtKey key) const;

		/*****
		Refs
		*****/

		GridPt& get_grid_point_ref(std::vector<int> grid_idxs);
		GridPt& get_grid_point_ref(IdxSet idx_set);
		GridPt& get_grid_point_ref(GridPtKey key);

		GridPtOut& get_grid_point_outside_ref(std::vector<int> grid_idxs);
		GridPtOut& get_grid_point_outside_ref(IdxSet idx_set);
		GridPtOut& get_grid_point_outside_ref(GridPtKey key);

		/********************
		Get grid points surrounding a point
		********************/

		Nbr2 get_surrounding_2_grid_pts(std::vector<double> abscissas) const;
		Nbr2 get_surrounding_2_grid_pts_by_ref(const std::vector<double>& abscissas) const;

		Nbr4 get_surrounding_4_grid_pts(std::vector<double> abscissas) const;
		Nbr4 get_surrounding_4_grid_pts_by_ref(const std::vector<double>& abscissas) const;

		/********************
		Get a point by interpolating
		********************/

		double get_val(std::vector<double> abscissas) const;
		double get_val_by_ref(const std::vector<double>& abscissas) const;

		/********************
		Get derivative wrt a point value p
		********************/

		double get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet4 idxs_k);
		double get_deriv_wrt_pt_value_by_ref(const std::vector<double>& abscissas, const IdxSet4& idxs_k);

		/********************
		Get derivative wrt x
		********************/

		double get_deriv_wrt_x(std::vector<double> abscissas, int k);
		double get_deriv_wrt_x_by_ref(const std::vector<double>& abscissas, int k);

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
		Copy ordinates from another grid
		********************/

	    // Warning: grid dimensions are not checked!
	    void copy_ordinates(const Grid* other);

		/********************
		Read/write grid
		********************/

		void read_from_file(std::string fname);
		void write_to_file(std::string fname) const;
	};

};