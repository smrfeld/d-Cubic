#include <vector>

#ifndef IDX_SET_H
#define IDX_SET_H
#include "idx_set.hpp"
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forwards
	class GridPt;

	/****************************************
	Neighborhood of points surrounding a point, 2 in each dim
	****************************************/

	class Nbr2 {

	private:

		// No dims
		int _no_dims;

		// The point
		IdxSet _idxs_i;

		// The grid pts
		int _no_grid_pts;
		GridPt** _grid_pts;

		// Frac abscissa
		double* _frac_abscissas;

		// Constructor helpers
		void _clean_up();
		void _copy(const Nbr2& other);
		void _move(Nbr2& other);
		void _shared_constructor();

	public:

		// Constructor
		Nbr2(IdxSet idxs_i);
		Nbr2(IdxSet idxs_i, double* frac_abscissas);
		Nbr2(IdxSet idxs_i, std::vector<double> frac_abscissas);
		Nbr2(const Nbr2& other);
		Nbr2(Nbr2&& other);
		Nbr2& operator=(const Nbr2 &other);
		Nbr2& operator=(Nbr2 &&other);
		~Nbr2();

		/********************
		Access
		********************/

		// Grid pts
		void set_grid_point(IdxSet idxs, GridPt *grid_pt);
		GridPt* get_grid_point(IdxSet idxs) const;
		GridPt* get_grid_point(int i) const;

		// Idx sets
		IdxSet convert_idx_set(int i) const;
		int convert_idx_set(IdxSet idxs) const;

		// No grid pts
		int get_no_grid_pts() const;
		
		// Frac abscissas
		double get_frac_abscissa(int dim) const;

		// Center idx
		IdxSet get_idxs_i() const;
		int get_idx_i(int dim) const;
	};



















































	/****************************************
	Neighborhood of points surrounding a point, 4 in each dim
	****************************************/
	
	class Nbr4 {

	private:

		// No dims
		int _no_dims;

		// The point
		IdxSet _idxs_i;

		// The grid pts
		int _no_grid_pts;
		GridPt** _grid_pts;

		// Frac abscissa
		double* _frac_abscissas;

		// Constructor helpers
		void _clean_up();
		void _copy(const Nbr4& other);
		void _move(Nbr4& other);
		void _shared_constructor();

	public:

		// Constructor
		Nbr4(IdxSet idxs_i);
		Nbr4(IdxSet idxs_i, double* frac_abscissas);
		Nbr4(IdxSet idxs_i, std::vector<double> frac_abscissas);
		Nbr4(const Nbr4& other);
		Nbr4(Nbr4&& other);
		Nbr4& operator=(const Nbr4 &other);
		Nbr4& operator=(Nbr4 &&other);
		~Nbr4();

		/********************
		Access
		********************/

		// Grid pts
		void set_grid_point(IdxSet idxs, GridPt *grid_pt);
		GridPt* get_grid_point(IdxSet idxs) const;
		GridPt* get_grid_point(int i) const;

		// No grid pts
		int get_no_grid_pts() const;

		// Idx sets
		IdxSet convert_idx_set(int i) const;
		int convert_idx_set(IdxSet idxs) const;

		// Frac abscissas
		double get_frac_abscissa(int dim) const;

		// Center idx
		IdxSet get_idxs_i() const;
		int get_idx_i(int dim) const;

		// Check if all pts are interior
		bool check_are_all_pts_inside() const;
	};
};