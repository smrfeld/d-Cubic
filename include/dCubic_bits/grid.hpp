#include <string>
#include <vector>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forwards
	class Dimension1D;
	class GridPt;
	class GridPtIn;
	class GridPtOut;
	class Nbr4;
	class Nbr2;
	class IdxSet;
	enum class Loc: unsigned int;

	/****************************************
	Grid
	****************************************/

	class Grid {

	private:

		// No dims
		int _no_dims;

		// Dims
		Dimension1D** _dims;

		// No pts in each dim
		int* _no_pts_in_dim;

		// Grid points
		GridPtIn** _grid_pts_in;
		GridPtOut** _grid_pts_out;

		// Total number grid pts
		int _no_grid_pts;

		// Iterating funcs
		void _iterate_make_grid_pt_inside(int dim, IdxSet idxs, double* abscissas);
		void _iterate_make_grid_pt_outside(int dim, IdxSet idxs, double* abscissas, Loc* locs);		
		void _iterate_get_surrounding_2_grid_pts(int dim, IdxSet idxs_local, const IdxSet &idxs_lower, Nbr2 &nbr2) const;
		void _iterate_get_surrounding_4_grid_pts(int dim, IdxSet idxs_local, const IdxSet &idxs_p0, Nbr4 &nbr4) const;
		double _iterate_interpolate(int delta, int d, Nbr4 &nbr4, IdxSet idxs_j) const;
		double _iterate_deriv_pt_value(int delta, int d, Nbr4 &nbr4, IdxSet idxs_j, IdxSet idxs_k) const;
		double _iterate_deriv_x(int delta, int k, int d, Nbr4 &nbr4, IdxSet idxs_j) const;

		// Constructor helpers
		void _clean_up();
		void _copy(const Grid& other);
		void _move(Grid& other);
		void _shared_constructor();

	public:

		/********************
		Constructor
		********************/

		// Note: ownership of Dimension1D is not transferred!
		Grid(int no_dims, Dimension1D** dims);
		Grid(std::vector<Dimension1D*> dims);
		Grid(const Grid& other);
		Grid(Grid&& other);
		Grid& operator=(const Grid &other);
		Grid& operator=(Grid &&other);
		~Grid();

		/********************
		Print
		********************/

		void print() const;

		/********************
		Get dims
		********************/

		int get_no_dims() const;
		Dimension1D* get_dim(int dim) const;

		/********************
		Get/set grid points
		********************/

		GridPt* get_grid_point(IdxSet idxs) const;
		GridPtIn* get_grid_point_inside(IdxSet idxs) const;
		GridPtOut* get_grid_point_outside(IdxSet idxs) const;

		// MOVE a grid point in
		// This means Grid class assumes ownership!
		void move_grid_point_inside(IdxSet idxs, GridPtIn* grid_pt);
		void move_grid_point_outside(IdxSet idxs, GridPtOut* grid_pt);

		/********************
		Get grid points surrounding a point
		********************/

		Nbr2 get_surrounding_2_grid_pts(double* abscissas) const;

		Nbr4 get_surrounding_4_grid_pts(double* abscissas) const;

		/********************
		Get a point by interpolating
		********************/

		double get_val(double* abscissas) const;

		/********************
		Get derivative wrt a point value p
		********************/

		double get_deriv_wrt_pt_value(double* abscissas, IdxSet idxs_k);

		/********************
		Get derivative wrt x
		********************/

		double get_deriv_wrt_x(double* abscissas, int k);

		/********************
		1D functions
		********************/

		double f1d_interpolate(double x_frac, double p0, double p1, double p2, double p3) const;

		// p_idx = 0,1,2,3, depending on loc
		double f1d_deriv_pt_value(double x_frac, int p_idx) const;

		double f1d_deriv_x(double x_frac, double p0, double p1, double p2, double p3) const;

	    /********************
		Apply M or P mappings
		********************/

	    IdxSet apply_m_mapping(IdxSet idxs) const;
	    IdxSet apply_p_mapping(IdxSet idxs) const;

		/********************
		Read/write grid
		********************/

		void read_from_file(std::string fname);
		void write_to_file(std::string fname) const;
	};

};