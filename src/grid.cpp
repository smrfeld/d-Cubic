#include "../include/dCubic_bits/grid.hpp"

// Other headers
#include "../include/dCubic_bits/dimension_1d.hpp"
#include "../include/dCubic_bits/grid_pt.hpp"
#include "../include/dCubic_bits/grid_pt_out.hpp"
#include "../include/dCubic_bits/grid_pt_key.hpp"
#include "../include/dCubic_bits/idx_set.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include "cmath"

#define DIAG_PROJ 0

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Grid implementation
	****************************************/

	class Grid::Impl {
	private:

		// Dimension1Dension of the grid
		int _dim_grid;

		// No pts in each dim of the grid
		std::vector<Dimension1D> _dims;		

		// Size of the grid
		int _no_pts_grid;

		// Grid points
		std::map<GridPtKey, GridPt*> _grid_pts;

		// Outside grid points
		std::map<GridPtKey, GridPtOut*> _grid_pts_out;

		/********************
		Make grid
		********************/

		void _make_grid_pts();

		void _iterate_make_grid_pt(IdxSet &grid_pt_idxs, int dim);

		/********************
		Make outside grid points
		********************/

		void _make_outside_grid_pts();

		void _iterate_make_grid_pt_outside(IdxSet &grid_pt_idxs, IdxSet &idxs_of_dims_outside, int dim);
		void _iterate_which_dim_are_outside(IdxSet &idxs_of_dims_outside, int idx, int no_dim_outside);

		/********************
		Get surrounding
		********************/

		void _iterate_get_surrounding_2_grid_pts(IdxSet &idxs_local, IdxSet &idxs_lower, IdxSet &idxs_upper, Nbr2 &map, int dim) const;

		void _iterate_get_surrounding_4_grid_pts(IdxSet &idxs_local, IdxSet &idxs_0, IdxSet &idxs_1, IdxSet &idxs_2, IdxSet &idxs_3, Nbr4 &nbr4, int dim) const;

		/********************
		Form the coeffs
		********************/

		// Start with coeff_p = 1.0
		// Dimension1D to interpolate in starts at _grid_dim-1
		// void _iterate_form_coeffs(std::map<GridPtKey, double> &coeffs_store, std::vector<double> &frac_abscissas, int dim_to_iterate_interpolate_in, IdxSet &idxs_p, double coeff_p);

		/********************
		Interpolate
		********************/

		double _iterate_interpolate(int delta, int d, std::vector<double> &frac_abscissas, IdxSet &idxs_j, Nbr4 &nbrs_p) const;

		/********************
		Constructor helpers
		********************/

		void _clean_up();
		void _copy(const Impl& other);
		void _move(Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(std::vector<Dimension1D> dims);
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl &other);
		Impl& operator=(Impl &&other);
		~Impl();

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





















	/****************************************
	Implementation
	****************************************/

	Grid::Impl::Impl(std::vector<Dimension1D> dims) {
		// Store
		_dim_grid = dims.size();
		_dims = dims;
		_no_pts_grid = 1;
		for (auto const &dim: _dims) {
			_no_pts_grid *= dim.get_no_pts();
		};

		// Make grid
		_make_grid_pts();
		_make_outside_grid_pts();
	};
	Grid::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	Grid::Impl::Impl(Impl&& other) {
		_move(other);
	};
    Grid::Impl& Grid::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Grid::Impl& Grid::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	Grid::Impl::~Impl()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void Grid::Impl::_clean_up()
	{
		// Nothing...
	};
	void Grid::Impl::_copy(const Impl& other)
	{
		_dim_grid = other._dim_grid;
		_dims = other._dims;
		_no_pts_grid = other._no_pts_grid;
		_grid_pts = other._grid_pts;
		_grid_pts_out = other._grid_pts_out;
	};
	void Grid::Impl::_move(Impl& other)
	{
		_copy(other);

		// Reset other
		other._dim_grid = 0;
		other._dims.clear();
		other._no_pts_grid = 0;
		other._grid_pts.clear();
		other._grid_pts_out.clear();
	};

	/********************
	Get dims
	********************/

	int Grid::Impl::get_no_dims() const {
		return _dims.size();
	};
	const std::vector<Dimension1D>& Grid::Impl::get_dims() const {
		return _dims;
	};

	/********************
	Make grid
	********************/

	void Grid::Impl::_make_grid_pts() {
		// Clear old
		_grid_pts.clear();

		// Iterate
		IdxSet grid_pt_idxs(_dim_grid);
		_iterate_make_grid_pt(grid_pt_idxs,0);
	};

	void Grid::Impl::_iterate_make_grid_pt(IdxSet &grid_pt_idxs, int dim) {
		if (dim != _dim_grid) {
			// Deeper!
			for (grid_pt_idxs[dim]=0; grid_pt_idxs[dim]<_dims[dim].get_no_pts(); grid_pt_idxs[dim]++) {
				_iterate_make_grid_pt(grid_pt_idxs, dim+1);
			};
		} else {
			// Do something

			// Make the abscissa
			std::vector<double> abscissas;
			for (auto dim2=0; dim2<_dim_grid; dim2++) {
				abscissas.push_back(_dims[dim2].get_pt_at_idx(grid_pt_idxs[dim2]));
			};

			// Make the grid pt
			_grid_pts[GridPtKey(grid_pt_idxs,_dims)] = new GridPt(grid_pt_idxs,abscissas);
		};
	};

	/********************
	Make outside grid points
	********************/

	void Grid::Impl::_make_outside_grid_pts() {
		// Different kinds of outside points:
		// Only 1 dim is outside
		// 2 dims are outside
		// etc

		// Iterate over how many dimensions are outside for each grid point
		// Min = 1
		// Max = _dim_grid
		for (int no_outside=1; no_outside <= _dim_grid; no_outside++) {
			// The idxs of the grid point that are outside
			IdxSet idxs_of_dims_outside(no_outside);
			// Iterate
			_iterate_which_dim_are_outside(idxs_of_dims_outside, 0, no_outside);
		};	
	};

	void Grid::Impl::_iterate_make_grid_pt_outside(IdxSet &grid_pt_idxs, IdxSet &idxs_of_dims_outside, int dim) {

		if (dim != _dim_grid) {
			// Further down!

			// Check if this dim is outside
			if (!idxs_of_dims_outside.find(dim)) {
				// Inside
				// Loop all interior pts
				for (grid_pt_idxs[dim]=0; grid_pt_idxs[dim]<_dims[dim].get_no_pts(); grid_pt_idxs[dim]++) {
					_iterate_make_grid_pt_outside(grid_pt_idxs,idxs_of_dims_outside,dim+1);
				};
			} else {
				// Outside
				// Loop just -1 and _dims[dim]->get_no_pts()
				grid_pt_idxs[dim] = -1;
				_iterate_make_grid_pt_outside(grid_pt_idxs,idxs_of_dims_outside,dim+1);
				grid_pt_idxs[dim] = _dims[dim].get_no_pts();
				_iterate_make_grid_pt_outside(grid_pt_idxs,idxs_of_dims_outside,dim+1);
			};

		} else {
			// Do something

			// Print the outside grid pt idxs
			// std::cout << grid_pt_idxs << std::endl;

			// Make the abscissa
			std::vector<double> abscissas;
			for (auto dim2=0; dim2<_dim_grid; dim2++) {
				abscissas.push_back(_dims[dim2].get_pt_at_idx(grid_pt_idxs[dim2]));
			};

			// Find the two pts
			IdxSet p1_idxs(_dim_grid), p2_idxs(_dim_grid);
			for (auto dim2=0; dim2<_dim_grid; dim2++) {
				if (grid_pt_idxs[dim2] == -1) {
					p1_idxs[dim2] = grid_pt_idxs[dim2] + 1;
					p2_idxs[dim2] = grid_pt_idxs[dim2] + 2;
				} else if (grid_pt_idxs[dim2] == _dims[dim2].get_no_pts()) {
					p1_idxs[dim2] = grid_pt_idxs[dim2] - 1;
					p2_idxs[dim2] = grid_pt_idxs[dim2] - 2;
				} else {
					p1_idxs[dim2] = grid_pt_idxs[dim2];
					p2_idxs[dim2] = grid_pt_idxs[dim2];
				};
			};
			const GridPt *p1 = get_grid_point(p1_idxs);
			const GridPt *p2 = get_grid_point(p2_idxs);

			// Make the outside grid point
			_grid_pts_out[GridPtKey(grid_pt_idxs,_dims)] = new GridPtOut(grid_pt_idxs,abscissas,p1,p2);
			// std::cout << "Made outside pt: " << grid_pt_idxs << " = " << _grid_pts_out[GridPtKey(grid_pt_idxs,GridPtType::OUTSIDE,_dims)]->print_abscissa() << std::endl;
		};
	};

	void Grid::Impl::_iterate_which_dim_are_outside(IdxSet &idxs_of_dims_outside, int idx, int no_dim_outside) {
		if (idx != no_dim_outside) {
			// Further down!

			int idxStart;
			if (idx == 0) {
				idxStart = 0;
			} else {
				idxStart = idxs_of_dims_outside[idx-1]+1; // no duplicates in for loop
				// This makes for loops like:
				// for (i=0; ...) 
				//    for (j = i+1; ....)
				//       ...
			};

			for (idxs_of_dims_outside[idx]=idxStart; idxs_of_dims_outside[idx]<_dim_grid; idxs_of_dims_outside[idx]++) {
				_iterate_which_dim_are_outside(idxs_of_dims_outside,idx+1,no_dim_outside);
			};

		} else {
			// Do something
			
			// Print which dims are outside
			/*
			std::cout << ">>> Dimension1Ds outside: ";
			for (auto j_outside=0; j_outside<idxs_of_dims_outside.size(); j_outside++) {
				std::cout << idxs_of_dims_outside[j_outside] << " ";
			};
			std::cout << "" << std::endl;
			*/

			// Go over all possible grid points satisfying having this many outside points
			IdxSet grid_pt_idxs(_dim_grid);
			_iterate_make_grid_pt_outside(grid_pt_idxs,idxs_of_dims_outside,0);
		};
	};

	/********************
	Get grid point
	********************/

	const std::map<GridPtKey, GridPt*>& Grid::Impl::get_grid_points() const {
		return _grid_pts;
	};
	const GridPt* Grid::Impl::get_grid_point(std::vector<int> grid_idxs) const {
		return get_grid_point(IdxSet(grid_idxs));
	};
	const GridPt* Grid::Impl::get_grid_point(IdxSet idx_set) const {
		return get_grid_point(GridPtKey(idx_set,_dims));
	};
	const GridPt* Grid::Impl::get_grid_point(GridPtKey key) const {
		auto it = _grid_pts.find(key);
		if (it == _grid_pts.end()) {
			return nullptr;
		};
		return it->second;
	};

	const std::map<GridPtKey, GridPtOut*>& Grid::Impl::get_grid_points_outside() const {
		return _grid_pts_out;
	};
	const GridPtOut* Grid::Impl::get_grid_point_outside(std::vector<int> grid_idxs) const {
		return get_grid_point_outside(IdxSet(grid_idxs));
	};
	const GridPtOut* Grid::Impl::get_grid_point_outside(IdxSet idx_set) const {
		return get_grid_point_outside(GridPtKey(idx_set,_dims));
	};
	const GridPtOut* Grid::Impl::get_grid_point_outside(GridPtKey key) const {
		auto it = _grid_pts_out.find(key);
		if (it == _grid_pts_out.end()) {
			return nullptr;
		};
		return it->second;	
	};

	/********************
	Set grid point values
	********************/

	void Grid::Impl::set_grid_point_ordinate(const GridPt* grid_pt, double val) {
		auto it = _grid_pts.find(GridPtKey(grid_pt->get_idxs(),_dims));
		if (it == _grid_pts.end()) {
			std::cerr << ">>> Error: Grid::Impl::set_grid_point_ordinate <<< could not find grid pt: " << grid_pt->print_abscissa() << std::endl;
			exit(EXIT_FAILURE);
		};
		it->second->set_ordinate(val);
	};

	/********************
	Get grid points surrounding a point
	********************/

	std::pair<Nbr2,std::vector<double>> Grid::Impl::get_surrounding_2_grid_pts(std::vector<double> abscissas) const {
		// Check size
		if (abscissas.size() != _dim_grid) {
			std::cerr << ">>> Error:Grid::Impl::get_surrounding_2_grid_pts <<< Abscissa size should equal grid size." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Frac abscissas
		std::vector<double> frac_abscissas;

		// Get bounding idxs
		IdxSet idxs_lower(_dim_grid), idxs_upper(_dim_grid);
		std::pair<bool,std::pair<int,int>> bounds;
		for (auto dim=0; dim<_dim_grid; dim++) {
			bounds = _dims[dim].get_surrounding_idxs(abscissas[dim]);
			if (!bounds.first) {
				// Outside grid
				std::cerr << ">>> Error:Grid::Impl::get_surrounding_2_grid_pts <<< Abscissa in dim: " << dim << " value: " << abscissas[dim] << " is outside the grid: " << _dims[dim].get_start_pt() << " to: " << _dims[dim].get_end_pt() << std::endl;
				exit(EXIT_FAILURE);
			};

			idxs_lower[dim] = bounds.second.first;
			idxs_upper[dim] = bounds.second.second;

			// Frac
			frac_abscissas.push_back((abscissas[dim] - _dims[dim].get_pt_at_idx(idxs_lower[dim])) / (_dims[dim].get_pt_at_idx(idxs_upper[dim]) - _dims[dim].get_pt_at_idx(idxs_lower[dim])));
		};

		// Iterate to fill out the map
		IdxSet idxs_local(_dim_grid);
		Nbr2 ret;
		_iterate_get_surrounding_2_grid_pts(idxs_local,idxs_lower,idxs_upper,ret,0);

		return std::make_pair(ret,frac_abscissas);
	};

	void Grid::Impl::_iterate_get_surrounding_2_grid_pts(IdxSet &idxs_local, IdxSet &idxs_lower, IdxSet &idxs_upper, Nbr2 &map, int dim) const {
		if (dim != _dim_grid) {
			// Deeper!
			// Can be lower (=0) or higher (=+1) in this dim
			idxs_local[dim] = 0;
			_iterate_get_surrounding_2_grid_pts(idxs_local,idxs_lower,idxs_upper,map,dim+1);
			idxs_local[dim] = 1;
			_iterate_get_surrounding_2_grid_pts(idxs_local,idxs_lower,idxs_upper,map,dim+1);

		} else {
			// Do something

			// Get grid point idxs
			IdxSet idxs_grid_pt(_dim_grid);
			for (auto dim2=0; dim2<_dim_grid; dim2++) {
				if (idxs_local[dim2] == 0) {
					idxs_grid_pt[dim2] = idxs_lower[dim2];
				} else if (idxs_local[dim2] == 1) {
					idxs_grid_pt[dim2] = idxs_upper[dim2];
				};
			};

			// Add to map
			map.in[GridPtKey(idxs_local,2)] = get_grid_point(idxs_grid_pt);
		};
	};

	std::pair<Nbr4,std::vector<double>> Grid::Impl::get_surrounding_4_grid_pts(std::vector<double> abscissas) const {
		// Check size
		if (abscissas.size() != _dim_grid) {
			std::cerr << ">>> Error:Grid::Impl::get_surrounding_4_grid_pts <<< Abscissa size should equal grid size." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Frac abscissas
		std::vector<double> frac_abscissas;

		// Get bounding idxs
		IdxSet idxs_0(_dim_grid), idxs_1(_dim_grid), idxs_2(_dim_grid), idxs_3(_dim_grid);
		std::pair<bool,std::pair<int,int>> bounds;
		for (auto dim=0; dim<_dim_grid; dim++) {
			bounds = _dims[dim].get_surrounding_idxs(abscissas[dim]);
			if (!bounds.first) {
				// Outside grid
				std::cerr << ">>> Error:Grid::Impl::get_surrounding_4_grid_pts <<< Abscissa in dim: " << dim << " value: " << abscissas[dim] << " is outside the grid: " << _dims[dim].get_start_pt() << " to: " << _dims[dim].get_end_pt() << std::endl;
				exit(EXIT_FAILURE);
			};

			idxs_1[dim] = bounds.second.first;
			idxs_2[dim] = bounds.second.second;
			idxs_0[dim] = idxs_1[dim]-1;
			idxs_3[dim] = idxs_2[dim]+1;

			// Frac
			frac_abscissas.push_back((abscissas[dim] - _dims[dim].get_pt_at_idx(idxs_1[dim])) / (_dims[dim].get_pt_at_idx(idxs_2[dim]) - _dims[dim].get_pt_at_idx(idxs_1[dim])));
		};

		// Iterate to fill out the map
		IdxSet idxs_local(_dim_grid);
		Nbr4 ret;
		_iterate_get_surrounding_4_grid_pts(idxs_local,idxs_0,idxs_1,idxs_2,idxs_3,ret,0);

		return std::make_pair(ret,frac_abscissas);
	};	

	void Grid::Impl::_iterate_get_surrounding_4_grid_pts(IdxSet &idxs_local, IdxSet &idxs_0, IdxSet &idxs_1, IdxSet &idxs_2, IdxSet &idxs_3, Nbr4 &nbr4, int dim) const {
		if (dim != _dim_grid) {
			// Deeper!
			// Can be lower (=0,1) or higher (=2,3) in this dim
			idxs_local[dim] = 0;
			_iterate_get_surrounding_4_grid_pts(idxs_local,idxs_0,idxs_1,idxs_2,idxs_3,nbr4,dim+1);
			idxs_local[dim] = 1;
			_iterate_get_surrounding_4_grid_pts(idxs_local,idxs_0,idxs_1,idxs_2,idxs_3,nbr4,dim+1);
			idxs_local[dim] = 2;
			_iterate_get_surrounding_4_grid_pts(idxs_local,idxs_0,idxs_1,idxs_2,idxs_3,nbr4,dim+1);
			idxs_local[dim] = 3;
			_iterate_get_surrounding_4_grid_pts(idxs_local,idxs_0,idxs_1,idxs_2,idxs_3,nbr4,dim+1);

		} else {
			// Do something

			// Get grid point idxs
			IdxSet idxs_grid_pt(_dim_grid);
			for (auto dim2=0; dim2<_dim_grid; dim2++) {
				if (idxs_local[dim2] == 0) {
					idxs_grid_pt[dim2] = idxs_0[dim2];
				} else if (idxs_local[dim2] == 1) {
					idxs_grid_pt[dim2] = idxs_1[dim2];
				} else if (idxs_local[dim2] == 2) {
					idxs_grid_pt[dim2] = idxs_2[dim2];
				} else if (idxs_local[dim2] == 3) {
					idxs_grid_pt[dim2] = idxs_3[dim2];
				};
			};

			// Check: is it inside or out?
			bool inside=true;
			for (auto dim2=0; dim2<_dim_grid; dim2++) {
				if (idxs_grid_pt[dim2] < 0 || idxs_grid_pt[dim2] > _dims[dim2].get_no_pts()-1) {
					// Out
					inside = false;
					break;
				};
			};

			// Add to nbr4
			GridPtKey key(idxs_local,4);
			if (inside) {
				nbr4.types[key] = GridPtType::INSIDE;
				nbr4.in[key] = get_grid_point(idxs_grid_pt);
			} else {
				nbr4.types[key] = GridPtType::OUTSIDE;
				nbr4.out[key] = get_grid_point_outside(idxs_grid_pt);
			};
		};
	};

	/********************
	Get a point by interpolating
	********************/

	double Grid::Impl::_iterate_interpolate(int delta, int d, std::vector<double> &frac_abscissas, IdxSet &idxs_j, Nbr4 &nbrs_p) const {

		if (delta == d) {
			// Arrived; regular 1D interpolation
			double p0,p1,p2,p3;

			idxs_j[0] = 0;
			GridPtKey key = GridPtKey(idxs_j,_dims);
			GridPtType type = nbrs_p.types[key];
			if (type == GridPtType::INSIDE) {
				p0 = nbrs_p.in[key]->get_ordinate();
			} else {
				p0 = nbrs_p.out[key]->get_ordinate();
			};

			idxs_j[0] = 1;
			key = GridPtKey(idxs_j,_dims);
			type = nbrs_p.types[key];
			if (type == GridPtType::INSIDE) {
				p1 = nbrs_p.in[key]->get_ordinate();
			} else {
				p1 = nbrs_p.out[key]->get_ordinate();
			};

			idxs_j[0] = 2;
			key = GridPtKey(idxs_j,_dims);
			type = nbrs_p.types[key];
			if (type == GridPtType::INSIDE) {
				p2 = nbrs_p.in[key]->get_ordinate();
			} else {
				p2 = nbrs_p.out[key]->get_ordinate();
			};

			idxs_j[0] = 3;
			key = GridPtKey(idxs_j,_dims);
			type = nbrs_p.types[key];
			if (type == GridPtType::INSIDE) {
				p3 = nbrs_p.in[key]->get_ordinate();
			} else {
				p3 = nbrs_p.out[key]->get_ordinate();
			};

			return interpolate_1d_by_ref(frac_abscissas[d-delta+1-1],p0,p1,p2,p3);
		} else {
			double p0,p1,p2,p3;

			// Deeper
			idxs_j[d-delta+1-1] = 0;
			p0 = _iterate_interpolate(delta, d, frac_abscissas, idxs_j, nbrs_p);

			idxs_j[d-delta+1-1] = 1;
			p1 = _iterate_interpolate(delta, d, frac_abscissas, idxs_j, nbrs_p);

			idxs_j[d-delta+1-1] = 2;
			p2 = _iterate_interpolate(delta, d, frac_abscissas, idxs_j, nbrs_p);			

			idxs_j[d-delta+1-1] = 3;
			p3 = _iterate_interpolate(delta, d, frac_abscissas, idxs_j, nbrs_p);

			return interpolate_1d_by_ref(frac_abscissas[d-delta+1-1],p0,p1,p2,p3);
		};
	};

	double Grid::Impl::get_val(std::vector<double> abscissas) const {
		// Delta
		int delta = 1;

		// Stopping
		int d = get_no_dims();

		// Idxs
		IdxSet idxs_j(d);

		// nbrs p and frac
		auto pr = get_surrounding_4_grid_pts(abscissas);
		Nbr4 nbrs_p = pr.first;
		std::vector<double> frac_abscissas = pr.second;

		// Iterate
		return _iterate_interpolate(delta,d,frac_abscissas,idxs_j,nbrs_p);
	};

	/********************
	Get derivative
	********************/

	double Grid::Impl::get_deriv_wrt_pt_value(std::vector<double> abscissas, std::vector<int> grid_idxs) {
		return 0.0;
	};
	double Grid::Impl::get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet idx_set) {
		return 0.0;
	};
	double Grid::Impl::get_deriv_wrt_pt_value(std::vector<double> abscissas, GridPtKey grid_pt_key) {
		return 0.0;
	};

	double Grid::Impl::get_deriv_wrt_x(std::vector<double> abscissas, int k) {
		return 0.0;
	};

	/********************
	1D funcs
	********************/

	double Grid::Impl::interpolate_1d(double x, double p0, double p1, double p2, double p3) const {
		return (-0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3)*pow(x,3) + (p0 - 2.5*p1 + 2.0*p2 - 0.5*p3)*pow(x,2) + (-0.5*p0 + 0.5*p2)*x + p1;
	};
	double Grid::Impl::interpolate_1d_by_ref(const double &x, const double &p0, const double &p1, const double &p2, const double &p3) const {
		return (-0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3)*pow(x,3) + (p0 - 2.5*p1 + 2.0*p2 - 0.5*p3)*pow(x,2) + (-0.5*p0 + 0.5*p2)*x + p1;
	};

	/********************
	Read/write grid
	********************/

	void Grid::Impl::read_from_file(std::string fname) {
		std::ofstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Grid::Impl::read_from_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all grid pts
		std::vector<double> abscissas;
		double ordinate;
		int i_count = 0;
		for (auto &pr: _grid_pts) {

			// Get
			abscissas = pr.second->get_abscissas();
			ordinate = pr.second->get_ordinate();

			// Write
			for (auto dim=0; dim<_dim_grid; dim++) {
				f << abscissas[dim] << " ";
			};
			f << ordinate;

			if (i_count != _grid_pts.size()-1) {
				f << "\n";
			};
			i_count++;
		};

		// Close
		f.close();
	};
	void Grid::Impl::write_to_file(std::string fname) const {
		std::ofstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Grid::Impl::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all grid pts
		std::vector<double> abscissas;
		double ordinate;
		int i_count = 0;
		for (auto &pr: _grid_pts) {

			// Get
			abscissas = pr.second->get_abscissas();
			ordinate = pr.second->get_ordinate();

			// Write
			for (auto dim=0; dim<_dim_grid; dim++) {
				f << abscissas[dim] << " ";
			};
			f << ordinate;

			if (i_count != _grid_pts.size()-1) {
				f << "\n";
			};
			i_count++;
		};

		// Close
		f.close();
	};

	/********************
	Project
	********************/

	/*
	void Grid::Impl::project() {

		if (DIAG_PROJ) {
			std::cout << ">>> Grid::Impl::project <<<" << std::endl;
		};

		// Entries in the vec/matrix
		std::map<GridPtKey, std::map<GridPtKey,double> > a_matr;
		std::map<GridPtKey, double> b_vec;

		// Global coeff map
		std::map<GridPtKey, double> global_coeffs;

		// Declarations
		std::vector<double> abscissas,abscissas_next,frac_abscissas;
		double path;
		std::map<GridPtKey, double> local_coeffs;
		std::shared_ptr<DataPt> data_pt, data_pt_next;
		Nbr4 nbr4;
		double coeff;
		IdxSet idxs_local(_dim_grid);
		IdxSet idxs_global(_dim_grid);
		IdxSet idxs_dep_1(_dim_grid),idxs_dep_2(_dim_grid);
		std::vector<GridPtKey> keys;
		double xk,xl;

		// Go through all data points except the last
		for (auto i_pt=0; i_pt<_data_pts.size()-1; i_pt++) {
			if (DIAG_PROJ) {
				std::cout << ">>> Grid::Impl::project <<< Doing pt: " << i_pt << " / " << _data_pts.size()-2 << std::endl;
			};

			data_pt = _data_pts[i_pt];
			data_pt_next = _data_pts[i_pt+1];

			// Get abscissas
			abscissas = data_pt->get_abscissas();
			abscissas_next = data_pt_next->get_abscissas();

			// Calculate path element
			path = 0.0;
			for (auto dim=0; dim<_dim_grid; dim++) {
				path += pow(abscissas_next[dim]-abscissas[dim],2);
			};
			path = sqrt(path);

			if (DIAG_PROJ) {
				std::cout << ">>> Grid::Impl::project <<<     Calculated path length: " << path << std::endl;
			};

			// Get frac abscissas
			frac_abscissas = data_pt->get_frac_abscissas();

			// Iterate such that local_coeffs holds the local map of coeffs for the surrounding nbr4
			local_coeffs.clear();
			_iterate_form_coeffs(local_coeffs,frac_abscissas,_dim_grid-1,idxs_local,1.0);

			if (DIAG_PROJ) {
				std::cout << ">>> Grid::Impl::project <<<     Local coeffs:" << std::endl;
				for (auto &pr: local_coeffs) {
					std::cout << "                             " << pr.first << " " << pr.second << std::endl;
				};
			};

			// Now grab the nbr4
			nbr4 = data_pt->get_nbr4();

			// Go through each point; add INSIDE to the global map
			// Skip OUTSIDE
			global_coeffs.clear();
			for (auto &pr: nbr4.types) {
				if (pr.second == GridPtType::INSIDE) {
					// Get global idxs
					idxs_global = nbr4.in[pr.first]->get_idxs();

					// Get coeff
					coeff = local_coeffs[pr.first];

					// std::cout << "Getting coeff of local idx: " << pr.first << " : " << coeff << std::endl;

					// Add to global map
					global_coeffs[GridPtKey(idxs_global,_dims)] = coeff;
				};
			};

			// Corrections for OUTSIDE pts
			for (auto &pr: nbr4.types) {
				if (pr.second == GridPtType::OUTSIDE) {
					// Need to correct the coeffs

					// Get the coeff from the local map
					coeff = local_coeffs[pr.first];

					// Get the idxs of the two dep pts
					idxs_dep_1 = nbr4.out[pr.first]->get_dep_p1()->get_idxs();
					idxs_dep_2 = nbr4.out[pr.first]->get_dep_p2()->get_idxs();

					// Adjust coeffs
					global_coeffs[GridPtKey(idxs_dep_1,_dims)] += 2.0 * coeff;
					global_coeffs[GridPtKey(idxs_dep_2,_dims)] -= coeff;
				};
			};

			if (DIAG_PROJ) {
				std::cout << ">>> Grid::Impl::project <<<     Got coeffs for all neighborhood pts" << std::endl;
			};

			// All keys in global_coeffs
			keys.clear();
			for (auto &pr: global_coeffs) {
				keys.push_back(pr.first);	
			};

			if (DIAG_PROJ) {
				std::cout << ">>> Grid::Impl::project <<<     Adding to matrix,vector...." << std::endl;
			};

			// Now add to the matrix/vector
			for (auto k=0; k<keys.size(); k++) {
				GridPtKey idx_k = keys[k];
				xk = global_coeffs[idx_k];

				// Add to vec
				auto it = b_vec.find(idx_k);
				if (it == b_vec.end()) {
					b_vec[idx_k] = xk * data_pt->get_ordinate() * path;
				} else {
					b_vec[idx_k] += xk * data_pt->get_ordinate() * path;
				};

				for (auto l=k; l<keys.size(); l++) {
					GridPtKey idx_l = keys[l];
					xl = global_coeffs[idx_l];

					// Add to matr
					auto it1 = a_matr.find(idx_k);
					if (it1 == a_matr.end()) {
						a_matr[idx_k][idx_l] = xk * xl * path;
						if (k!=l) {
							a_matr[idx_l][idx_k] = xk * xl * path;
						};
					} else {
						auto it2 = it1->second.find(idx_k);
						if (it2 == it1->second.end()) {
							a_matr[idx_k][idx_l] = xk * xl * path;
							if (k!=l) {
								a_matr[idx_l][idx_k] = xk * xl * path;
							};
						} else {
							a_matr[idx_k][idx_l] += xk * xl * path;
							if (k!=l) {
								a_matr[idx_l][idx_k] += xk * xl * path;
							};
						};
					};
				};
			};
		};

		if (DIAG_PROJ) {
			std::cout << ">>> Grid::Impl::project <<< Done looping through data." << std::endl;
			std::cout << ">>> Grid::Impl::project <<< The A matrix size is: " << a_matr.size() << std::endl;
			std::cout << ">>> Grid::Impl::project <<< The b vector size is: " << b_vec.size() << std::endl;
			std::cout << ">>> Grid::Impl::project <<< Now constructing the armadillo objects!" << std::endl;
		};

		// Get size
		int n = a_matr.size();
		if (a_matr.size() != b_vec.size()) {
			std::cerr << ">>> Error: Grid::Impl::project <<< b vec size = " << b_vec.size() << " but A matrix size = " << a_matr.size() << std::endl;
			exit(EXIT_FAILURE);
		};

		// arma::sp_mat A(n,n);
		arma::mat A(n,n);
		arma::vec b(n),x(n);

		// Grab all the keys
		std::vector<GridPtKey> keys2;
		for (auto &pr: b_vec) {
			keys2.push_back(pr.first);
		};

		// Make the matrix,vector
		for (auto i=0; i<keys2.size(); i++) {

			// vec
			b(i) = b_vec[keys2[i]];

			for (auto j=0; j<keys2.size(); j++) {

				// matr
				A(i,j) = a_matr[keys2[i]][keys2[j]];

			};
		};

		if (DIAG_PROJ) {
			A.print("A:");
			b.print("b:");
			std::cout << ">>> Grid::Impl::project <<< Solving..." << std::endl;
		};

		// Solve
		// x = spsolve(A, b, "lapack");
		x = solve(A,b);

		if (DIAG_PROJ) {
			std::cout << ">>> Grid::Impl::project <<< Done!" << std::endl;
			x.print("x:");
			std::cout << ">>> Grid::Impl::project <<< Next: setting in grid" << std::endl;
		};

		// Set solution in the grid
		for (auto &pr: _grid_pts) {
			auto it = std::find(keys2.begin(), keys2.end(), pr.first);
			if (it == keys2.end()) {
				// Not found -> 0
				pr.second->set_ordinate(0.0);
			} else {
				// Found; set from x
				int i = it - keys2.begin();
				pr.second->set_ordinate(x(i));
			};
		};

		if (DIAG_PROJ) {
			std::cout << ">>> Grid::Impl::project <<< Finished!" << std::endl;
		};

	};

	void Grid::Impl::_iterate_form_coeffs(std::map<GridPtKey, double> &coeffs_store, std::vector<double> &frac_abscissas, int dim_to_iterate_interpolate_in, IdxSet &idxs_p, double coeff_p) {

		if (dim_to_iterate_interpolate_in != -1) {
			// Deeper!

			// Grab the pt in this dim
			double x = frac_abscissas[dim_to_iterate_interpolate_in];

			// p0
			idxs_p[dim_to_iterate_interpolate_in] = 0;		
			_iterate_form_coeffs( coeffs_store, frac_abscissas, dim_to_iterate_interpolate_in-1, idxs_p, coeff_p*(-0.5*x+pow(x,2)-0.5*pow(x,3)) );

			// p1
			idxs_p[dim_to_iterate_interpolate_in] = 1;		
			_iterate_form_coeffs( coeffs_store, frac_abscissas, dim_to_iterate_interpolate_in-1, idxs_p, coeff_p*(1.0-2.5*pow(x,2)+1.5*pow(x,3)) );

			// p2
			idxs_p[dim_to_iterate_interpolate_in] = 2;		
			_iterate_form_coeffs( coeffs_store, frac_abscissas, dim_to_iterate_interpolate_in-1, idxs_p, coeff_p*(0.5*x + 2.0 * pow(x,2) -1.5*pow(x,3)) );

			// p3
			idxs_p[dim_to_iterate_interpolate_in] = 3;		
			_iterate_form_coeffs( coeffs_store, frac_abscissas, dim_to_iterate_interpolate_in-1, idxs_p, coeff_p*(-0.5*pow(x,2)+0.5*pow(x,3)) );

		} else {
			// Do something

			// Print the idx and coeff
			// std::cout << idxs_p << " coeff: " << coeff_p << std::endl;

			// Add to the map
			coeffs_store[GridPtKey(idxs_p,4)] = coeff_p;
		};
	};
	*/






































	/****************************************
	Impl forwards
	****************************************/

	/********************
	Constructor
	********************/

	Grid::Grid(std::vector<Dimension1D> dims) : _impl(new Impl(dims)) {};
	Grid::Grid(const Grid& other) : _impl(new Impl(*other._impl)) {};
	Grid::Grid(Grid&& other) : _impl(std::move(other._impl)) {};
	Grid& Grid::operator=(const Grid &other) {
        _impl.reset(new Impl(*other._impl));
        return *this; 
	};
	Grid& Grid::operator=(Grid &&other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	Grid::~Grid() = default;

	/********************
	Get dims
	********************/

	int Grid::get_no_dims() const {
		return _impl->get_no_dims();
	};
	const std::vector<Dimension1D>& Grid::get_dims() const {
		return _impl->get_dims();
	};

	/********************
	Get grid pts
	********************/

	const std::map<GridPtKey, GridPt*>& Grid::get_grid_points() const {
		return _impl->get_grid_points();
	};
	const GridPt* Grid::get_grid_point(std::vector<int> grid_idxs) const {
		return _impl->get_grid_point(grid_idxs);
	};
	const GridPt* Grid::get_grid_point(IdxSet idx_set) const {
		return _impl->get_grid_point(idx_set);
	};
	const GridPt* Grid::get_grid_point(GridPtKey key) const {
		return _impl->get_grid_point(key);
	};

	const std::map<GridPtKey, GridPtOut*>& Grid::get_grid_points_outside() const {
		return _impl->get_grid_points_outside();
	};
	const GridPtOut* Grid::get_grid_point_outside(std::vector<int> grid_idxs) const {
		return _impl->get_grid_point_outside(grid_idxs);
	};
	const GridPtOut* Grid::get_grid_point_outside(IdxSet idx_set) const {
		return _impl->get_grid_point_outside(idx_set);
	};
	const GridPtOut* Grid::get_grid_point_outside(GridPtKey key) const {
		return _impl->get_grid_point_outside(key);
	};

	/********************
	Set grid point values
	********************/

	void Grid::set_grid_point_ordinate(const GridPt* grid_pt, double val) {
		_impl->set_grid_point_ordinate(grid_pt,val);
	};

	/********************
	Get grid points surrounding a point
	********************/

	std::pair<Nbr2,std::vector<double>> Grid::get_surrounding_2_grid_pts(std::vector<double> abscissas) const {
		return _impl->get_surrounding_2_grid_pts(abscissas);
	};
	std::pair<Nbr4,std::vector<double>> Grid::get_surrounding_4_grid_pts(std::vector<double> abscissas) const {
		return _impl->get_surrounding_4_grid_pts(abscissas);
	};

	/********************
	Get a point by interpolating
	********************/

	double Grid::get_val(std::vector<double> abscissas) const {
		return _impl->get_val(abscissas);
	};

	/********************
	Get derivative
	********************/

	double Grid::get_deriv_wrt_pt_value(std::vector<double> abscissas, std::vector<int> grid_idxs) {
		return _impl->get_deriv_wrt_pt_value(abscissas,grid_idxs);
	};
	double Grid::get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet idx_set) {
		return _impl->get_deriv_wrt_pt_value(abscissas,idx_set);
	};
	double Grid::get_deriv_wrt_pt_value(std::vector<double> abscissas, GridPtKey grid_pt_key) {
		return _impl->get_deriv_wrt_pt_value(abscissas,grid_pt_key);
	};

	double Grid::get_deriv_wrt_x(std::vector<double> abscissas, int k) {
		return _impl->get_deriv_wrt_x(abscissas,k);
	};

	/********************
	1D funcs
	********************/

	double Grid::interpolate_1d(double x, double p0, double p1, double p2, double p3) const {
		return interpolate_1d(x,p0,p1,p2,p3);
	};
	double Grid::interpolate_1d_by_ref(const double &x, const double &p0, const double &p1, const double &p2, const double &p3) const {
		return interpolate_1d_by_ref(x,p0,p1,p2,p3);
	};


	/********************
	Read/write grid
	********************/

	void Grid::read_from_file(std::string fname) {
		_impl->read_from_file(fname);
	};
	void Grid::write_to_file(std::string fname) const {
		_impl->write_to_file(fname);
	};
};