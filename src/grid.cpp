#include "../include/dcubic_bits/grid.hpp"

// Other headers
#include "../include/dcubic_bits/dimension_1d.hpp"
#include "../include/dcubic_bits/grid_pt.hpp"
#include "../include/dcubic_bits/grid_pt_out.hpp"
#include "../include/dcubic_bits/grid_pt_key.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include "cmath"
#include <unordered_map>

#define DIAG_PROJ 0

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Hash for an unordered map
	****************************************/

	// Hash
	struct hash_gpk {
	    size_t operator() (const GridPtKey &grid_pt_key ) const {
	        return std::hash<int>()(grid_pt_key.get_linear());
	    };
	};

	/****************************************
	Neighborhood of points surrounding a point, 2 in each dim
	****************************************/

    Nbr2::Nbr2(IdxSet idxs_i) : idxs_i(idxs_i) {};

	/****************************************
	Neighborhood of points surrounding a point, 4 in each dim
	****************************************/

    Nbr4::Nbr4(IdxSet idxs_i) : idxs_i(idxs_i) {};

	/****************************************
	Ordinates associated with a Nbr2
	****************************************/

	P2::P2(Nbr2 nbr2) {
		for (auto const &pr: nbr2.in) {
			p[pr.first] = pr.second->get_ordinate();
		};
	};

	/****************************************
	Ordinates associated with a Nbr4
	****************************************/

	P4::P4(Nbr4 nbr4) {
		for (auto const &pr: nbr4.types) {
			if (pr.second == GridPtType::INSIDE) {
				p[pr.first] = nbr4.in[pr.first]->get_ordinate();
			} else {
				p[pr.first] = nbr4.out[pr.first]->get_ordinate();
			};
		};
	};




























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
		std::unordered_map<GridPtKey, GridPt*, hash_gpk> _grid_pts;

		// Outside grid points
		std::unordered_map<GridPtKey, GridPtOut*, hash_gpk> _grid_pts_out;

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

		void _iterate_get_surrounding_2_grid_pts(IdxSet2 &idxs_local, IdxSet &idxs_lower, IdxSet &idxs_upper, Nbr2 &nbr2, int dim) const;

		void _iterate_get_surrounding_4_grid_pts(IdxSet4 &idxs_local, IdxSet &idxs_0, IdxSet &idxs_1, IdxSet &idxs_2, IdxSet &idxs_3, Nbr4 &nbr4, int dim) const;

		/********************
		Interpolate
		********************/

		double _iterate_interpolate(int delta, int d, std::vector<double> &frac_abscissas, IdxSet4 &idxs_j, P4 &p) const;

		/********************
		Deriv wrt p
		********************/

		double _iterate_deriv_pt_value(int delta, int d, std::vector<double> &frac_abscissas, IdxSet4 &idxs_j, P4 &p, IdxSet &idxs_i, const IdxSet4 &idxs_k) const;

		/********************
		Derivative wrt x
		********************/

		double _iterate_deriv_x(int delta, int k, int d, std::vector<double> &frac_abscissas, IdxSet4 &idxs_j, P4 &p) const;

		/********************
		Constructor helpers
		********************/

		void _shared_constructor();
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
		
		const GridPt* get_grid_point(IdxSet idx_set) const;

		const GridPtOut* get_grid_point_outside(IdxSet idx_set) const;

		/*****
		Refs
		*****/

		GridPt& get_grid_point_ref(IdxSet idx_set);

		GridPtOut& get_grid_point_outside_ref(IdxSet idx_set);

		/********************
		Set grid point values
		********************/

		void set_grid_point_ordinate(const GridPt* grid_pt, double val);

		/********************
		Get grid points surrounding a point
		********************/

		Nbr2 get_surrounding_2_grid_pts(std::vector<double> abscissas) const;
		Nbr4 get_surrounding_4_grid_pts(std::vector<double> abscissas) const;

		Nbr2 get_surrounding_2_grid_pts_by_ref(const std::vector<double>& abscissas) const;
		Nbr4 get_surrounding_4_grid_pts_by_ref(const std::vector<double>& abscissas) const;

		/********************
		Get a point by interpolating
		********************/

		double get_val(std::vector<double> abscissas) const;
		double get_val_by_ref(const std::vector<double>& abscissas) const;

		/********************
		Get derivative
		********************/

		double get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet4 idx_set);
		double get_deriv_wrt_pt_value_by_ref(const std::vector<double>& abscissas, const IdxSet4& idxs_k);

		/********************
		Get derivative wrt x
		********************/

		double get_deriv_wrt_x(std::vector<double> abscissas, int k);
		double get_deriv_wrt_x_by_ref(const std::vector<double>& abscissas, int k);

		/********************
		1D funcs
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



















































	/****************************************
	Implementation
	****************************************/

	Grid::Impl::Impl(std::vector<Dimension1D> dims) {
		_dims = dims;
		_shared_constructor();
	};
	void Grid::Impl::_shared_constructor() {
		_dim_grid = _dims.size();
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
	Print
	********************/

	void Grid::Impl::print_grid_pts_inside() const {
		for (auto const &pr: _grid_pts) {
			std::cout << pr.first << " : " << pr.second->get_ordinate() << std::endl;
		};
	};

	/********************
	Helpers for constructors
	********************/

	void Grid::Impl::_clean_up()
	{
		for (auto &pr: _grid_pts) {
			delete pr.second;
			pr.second = nullptr;
		};
		for (auto &pr: _grid_pts_out) {
			delete pr.second;
			pr.second = nullptr;
		};
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
				abscissas.push_back(_dims[dim2].get_pt_by_idx(grid_pt_idxs[dim2]));
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
				// Loop just -1 and _dims[dim].get_no_pts()
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
				abscissas.push_back(_dims[dim2].get_pt_by_idx(grid_pt_idxs[dim2],false));
			};

			// Find the two pts
			IdxSet p1_idxs(_dim_grid), p2_idxs(_dim_grid);
			std::vector<LocInDim> locs;
			for (auto dim2=0; dim2<_dim_grid; dim2++) {
				if (grid_pt_idxs[dim2] == -1) {
					locs.push_back(LocInDim::P0_OUTSIDE);
					p1_idxs[dim2] = grid_pt_idxs[dim2] + 1;
					p2_idxs[dim2] = grid_pt_idxs[dim2] + 2;
				} else if (grid_pt_idxs[dim2] == _dims[dim2].get_no_pts()) {
					locs.push_back(LocInDim::P3_OUTSIDE);
					p1_idxs[dim2] = grid_pt_idxs[dim2] - 1;
					p2_idxs[dim2] = grid_pt_idxs[dim2] - 2;
				} else {
					locs.push_back(LocInDim::INSIDE);
					p1_idxs[dim2] = grid_pt_idxs[dim2];
					p2_idxs[dim2] = grid_pt_idxs[dim2];
				};
			};
			const GridPt *p1 = get_grid_point(p1_idxs);
			const GridPt *p2 = get_grid_point(p2_idxs);

			// Make the outside grid point
			_grid_pts_out[GridPtKey(grid_pt_idxs,_dims)] = new GridPtOut(grid_pt_idxs,abscissas,p1,p2,locs);
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

	/*****
	Const ptrs
	*****/

	const GridPt* Grid::Impl::get_grid_point(IdxSet idx_set) const {
		auto it = _grid_pts.find(GridPtKey(idx_set,_dims));
		if (it == _grid_pts.end()) {
			return nullptr;
		};
		return it->second;
	};

	const GridPtOut* Grid::Impl::get_grid_point_outside(IdxSet idx_set) const {
		auto it = _grid_pts_out.find(GridPtKey(idx_set,_dims));
		if (it == _grid_pts_out.end()) {
			return nullptr;
		};
		return it->second;	
	};

	/*****
	Refs
	*****/

	GridPt& Grid::Impl::get_grid_point_ref(IdxSet idx_set) {
		auto it = _grid_pts.find(GridPtKey(idx_set,_dims));
		if (it == _grid_pts.end()) {
			std::cerr << ">>> Error: Grid::Impl::get_grid_point_ref <<< key not found" << std::endl;
			exit(EXIT_FAILURE);
		};
		return *it->second;
	};

	GridPtOut& Grid::Impl::get_grid_point_outside_ref(IdxSet idx_set) {
		auto it = _grid_pts_out.find(GridPtKey(idx_set,_dims));
		if (it == _grid_pts_out.end()) {
			std::cerr << ">>> Error: Grid::Impl::get_grid_point_outside_ref <<< key not found" << std::endl;
			exit(EXIT_FAILURE);
		};
		return *it->second;
	};

	/********************
	Get grid points surrounding a point
	********************/

	Nbr2 Grid::Impl::get_surrounding_2_grid_pts(std::vector<double> abscissas) const {
		return get_surrounding_2_grid_pts_by_ref(abscissas);
	};
	Nbr4 Grid::Impl::get_surrounding_4_grid_pts(std::vector<double> abscissas) const {
		return get_surrounding_4_grid_pts_by_ref(abscissas);
	};

	Nbr2 Grid::Impl::get_surrounding_2_grid_pts_by_ref(const std::vector<double> &abscissas) const {
		// Check size
		if (abscissas.size() != _dim_grid) {
			std::cerr << ">>> Error:Grid::Impl::get_surrounding_2_grid_pts <<< Abscissa size should equal grid size." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Frac abscissas
		std::vector<double> frac_abscissas;

		// Get bounding idxs
		IdxSet idxs_lower(_dim_grid), idxs_upper(_dim_grid);
		int i;
		for (auto dim=0; dim<_dim_grid; dim++) {
			// Check in dim
			if (!_dims[dim].check_if_pt_is_inside_domain(abscissas[dim])) {
				// Outside grid
				std::cerr << ">>> Error:Grid::Impl::get_surrounding_2_grid_pts <<< Abscissa in dim: " << dim << " value: " << abscissas[dim] << " is outside the grid: " << _dims[dim].get_min_pt() << " to: " << _dims[dim].get_max_pt() << std::endl;
				exit(EXIT_FAILURE);
			};

			i = _dims[dim].get_idxs_surrounding_pt(abscissas[dim]);

			idxs_lower[dim] = i;
			idxs_upper[dim] = i+1;

			// Frac
			frac_abscissas.push_back((abscissas[dim] - _dims[dim].get_pt_by_idx(idxs_lower[dim])) / (_dims[dim].get_pt_by_idx(idxs_upper[dim]) - _dims[dim].get_pt_by_idx(idxs_lower[dim])));
		};

		// Returned
		Nbr2 nbr2(idxs_lower);
		nbr2.frac_abscissas = frac_abscissas;

		// Iterate to fill out the map
		IdxSet2 idxs_local(_dim_grid);
		_iterate_get_surrounding_2_grid_pts(idxs_local,idxs_lower,idxs_upper,nbr2,0);

		return nbr2;
	};

	void Grid::Impl::_iterate_get_surrounding_2_grid_pts(IdxSet2 &idxs_local, IdxSet &idxs_lower, IdxSet &idxs_upper, Nbr2 &nbr2, int dim) const {
		if (dim != _dim_grid) {
			// Deeper!
			// Can be lower (=0) or higher (=+1) in this dim
			idxs_local[dim] = 0;
			_iterate_get_surrounding_2_grid_pts(idxs_local,idxs_lower,idxs_upper,nbr2,dim+1);
			idxs_local[dim] = 1;
			_iterate_get_surrounding_2_grid_pts(idxs_local,idxs_lower,idxs_upper,nbr2,dim+1);

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

			// Add to nbr2
			nbr2.in[idxs_local] = get_grid_point(idxs_grid_pt);
		};
	};

	Nbr4 Grid::Impl::get_surrounding_4_grid_pts_by_ref(const std::vector<double> &abscissas) const {
		// Check size
		if (abscissas.size() != _dim_grid) {
			std::cerr << ">>> Error:Grid::Impl::get_surrounding_4_grid_pts <<< Abscissa size should equal grid size." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Frac abscissas
		std::vector<double> frac_abscissas;

		// Get bounding idxs
		IdxSet idxs_0(_dim_grid), idxs_1(_dim_grid), idxs_2(_dim_grid), idxs_3(_dim_grid);
		int i;
		for (auto dim=0; dim<_dim_grid; dim++) {
			// Check in dim
			if (!_dims[dim].check_if_pt_is_inside_domain(abscissas[dim])) {
				// Outside grid
				std::cerr << ">>> Error:Grid::Impl::get_surrounding_4_grid_pts <<< Abscissa in dim: " << dim << " value: " << abscissas[dim] << " is outside the grid: " << _dims[dim].get_min_pt() << " to: " << _dims[dim].get_max_pt() << std::endl;
				exit(EXIT_FAILURE);
			};

			i = _dims[dim].get_idxs_surrounding_pt(abscissas[dim]);

			idxs_0[dim] = i-1;
			idxs_1[dim] = i;
			idxs_2[dim] = i+1;
			idxs_3[dim] = i+2;

			// Frac
			frac_abscissas.push_back((abscissas[dim] - _dims[dim].get_pt_by_idx(idxs_1[dim])) / (_dims[dim].get_pt_by_idx(idxs_2[dim]) - _dims[dim].get_pt_by_idx(idxs_1[dim])));
		};

		// Returned
		Nbr4 nbr4(idxs_1);
		nbr4.frac_abscissas = frac_abscissas;

		// Iterate to fill out the nbr4
		IdxSet4 idxs_local(_dim_grid);
		_iterate_get_surrounding_4_grid_pts(idxs_local,idxs_0,idxs_1,idxs_2,idxs_3,nbr4,0);

		return nbr4;
	};	

	void Grid::Impl::_iterate_get_surrounding_4_grid_pts(IdxSet4 &idxs_local, IdxSet &idxs_0, IdxSet &idxs_1, IdxSet &idxs_2, IdxSet &idxs_3, Nbr4 &nbr4, int dim) const {
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
			if (inside) {
				nbr4.types[idxs_local] = GridPtType::INSIDE;
				nbr4.in[idxs_local] = get_grid_point(idxs_grid_pt);
			} else {
				nbr4.types[idxs_local] = GridPtType::OUTSIDE;
				nbr4.out[idxs_local] = get_grid_point_outside(idxs_grid_pt);
			};
		};
	};

	/********************
	Get a point by interpolating
	********************/

	double Grid::Impl::_iterate_interpolate(int delta, int d, std::vector<double> &frac_abscissas, IdxSet4 &idxs_j, P4 &p) const {
		double p0,p1,p2,p3;

		if (delta == d) {
			
			// Arrived; return the point itself
			return p.p[idxs_j];

		} else {

			// Deeper
			idxs_j[delta] = 0;
			p0 = _iterate_interpolate(delta+1, d, frac_abscissas, idxs_j, p);

			idxs_j[delta] = 1;
			p1 = _iterate_interpolate(delta+1, d, frac_abscissas, idxs_j, p);

			idxs_j[delta] = 2;
			p2 = _iterate_interpolate(delta+1, d, frac_abscissas, idxs_j, p);

			idxs_j[delta] = 3;
			p3 = _iterate_interpolate(delta+1, d, frac_abscissas, idxs_j, p);

			return f1d_interpolate_by_ref(frac_abscissas[delta],p0,p1,p2,p3);
		};
	};

	double Grid::Impl::get_val(std::vector<double> abscissas) const {
		return get_val_by_ref(abscissas);
	};

	double Grid::Impl::get_val_by_ref(const std::vector<double> &abscissas) const {
		// Stopping
		int d = get_no_dims();

		// Idxs
		IdxSet4 idxs_j(d);

		// nbrs p
		Nbr4 nbr4 = get_surrounding_4_grid_pts_by_ref(abscissas);

		// Frac
		std::vector<double> frac_abscissas = nbr4.frac_abscissas;

		// Convert nbr to p
		P4 p(nbr4);

		// Iterate
		return _iterate_interpolate(0,d,frac_abscissas,idxs_j,p);
	};

	/********************
	Get derivative
	********************/

	double Grid::Impl::_iterate_deriv_pt_value(int delta, int d, std::vector<double> &frac_abscissas, IdxSet4 &idxs_j, P4 &p, IdxSet &idxs_i, const IdxSet4 &idxs_k) const {

		if (delta == d) {
			// Done; evaluate
			// Check main condition
			bool at_least_one_cond_met = false;
			for (auto dim=0; dim<d; dim++) {
				if ((idxs_i[dim] == 0 && idxs_j[dim] == 0) || (idxs_i[dim] == _dims[dim].get_no_pts()-2 && idxs_j[dim] == 3)) {
					at_least_one_cond_met = true;
					break;
				};
			};

			if (at_least_one_cond_met) {

				// Yes, at least one
			

				// Form the new idxs
				IdxSet idxs_j_global = idxs_j - 1 + idxs_i;
				IdxSet idxs_k_global = idxs_k - 1 + idxs_i;

				// Apply M mapping
				IdxSet idxs_m = apply_m_mapping_by_ref(idxs_j_global);

				// Check
				if (idxs_m == idxs_k_global) {
					return 2.0;
				};

				// Apply P mapping
				IdxSet idxs_p = apply_p_mapping_by_ref(idxs_j_global);

				// Check
				if (idxs_p == idxs_k_global) {
					return -1.0;
				};

				// Failed
				// This means that this point = idxs_j is an outside grid point, but does not depend on the point we are differentiating WRT
				return 0.0;

			} else {

				// Not at_least_one_cond_met - evaluate the delta
				
				if (idxs_j == idxs_k) {
					return 1.0;
				} else {
					return 0.0;
				};

			};

		} else {
			// Go deeper

			idxs_j[delta] = 0;
			double dp0dp = _iterate_deriv_pt_value(delta+1,d,frac_abscissas,idxs_j,p,idxs_i,idxs_k);
			double dfdp0 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],0);

			idxs_j[delta] = 1;
			double dp1dp = _iterate_deriv_pt_value(delta+1,d,frac_abscissas,idxs_j,p,idxs_i,idxs_k);
			double dfdp1 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],1);

			idxs_j[delta] = 2;
			double dp2dp = _iterate_deriv_pt_value(delta+1,d,frac_abscissas,idxs_j,p,idxs_i,idxs_k);
			double dfdp2 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],2);

			idxs_j[delta] = 3;
			double dp3dp = _iterate_deriv_pt_value(delta+1,d,frac_abscissas,idxs_j,p,idxs_i,idxs_k);
			double dfdp3 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],3);

			return dfdp0 * dp0dp + dfdp1 * dp1dp + dfdp2 * dp2dp + dfdp3 * dp3dp;
		};
	};

	double Grid::Impl::get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet4 idxs_k) {
		return get_deriv_wrt_pt_value_by_ref(abscissas,idxs_k);
	};

	double Grid::Impl::get_deriv_wrt_pt_value_by_ref(const std::vector<double>& abscissas, const IdxSet4& idxs_k) {
		// d
		int d = get_no_dims();

		// nbrs
		Nbr4 nbr4 = get_surrounding_4_grid_pts_by_ref(abscissas);

		// Frac abscissas
		std::vector<double> frac_abscissas = nbr4.frac_abscissas;

		// Convert to p
		P4 p(nbr4);

		// Check: are there any exterior grid point?
		bool all_inside = true;
		for (auto dim=0; dim<d; dim++) {
			if (nbr4.idxs_i[dim] == 0 || nbr4.idxs_i[dim] == _dims[dim].get_no_pts()-2 ) {
				all_inside = false;
				break;
			};
		};

		if (all_inside) {
			// Case 1: Totally interior point; no dimension near boundary

			double ret=1.0;
			for (auto alpha=0; alpha<d; alpha++) {
				ret *= f1d_deriv_pt_value_by_ref(frac_abscissas[alpha], idxs_k[alpha]);
			};
			return ret;

		} else {
			// Case 2: At least one dimension near boundary

			// Check that we are not taking an illegal derivative
			// i.e. if i=0, j=0 => illegal
			// if i=n-2, j=3 => illegal
			for (auto dim=0; dim<d; dim++) {
				if ((nbr4.idxs_i[dim] == 0 && idxs_k[dim] == 0) || (nbr4.idxs_i[dim] == _dims[dim].get_no_pts()-2 && idxs_k[dim] == 3)) {
					std::cerr << ">>> Error: Grid::Impl::get_deriv_wrt_pt_value <<< The pt specified to take a derivative with respect to is not a real point (an exterior point that is estimated by a linear approx) - this is not allowed!" << std::endl;
					exit(EXIT_FAILURE);
				};
			};

			// Init idx set
			IdxSet4 idxs_j(d);

			// Iterate
			return _iterate_deriv_pt_value(0, d,frac_abscissas, idxs_j, p, nbr4.idxs_i, idxs_k);
		};
	};

	/********************
	Get derivative wrt x
	********************/

	double Grid::Impl::_iterate_deriv_x(int delta, int k, int d, std::vector<double> &frac_abscissas, IdxSet4 &idxs_j, P4 &p) const {
		
		if (delta == k) { // Evaluate the derivative

			idxs_j[k] = 0;
			double p0 = _iterate_interpolate(k+1,d,frac_abscissas,idxs_j,p);

			idxs_j[k] = 1;
			double p1 = _iterate_interpolate(k+1,d,frac_abscissas,idxs_j,p);

			idxs_j[k] = 2;
			double p2 = _iterate_interpolate(k+1,d,frac_abscissas,idxs_j,p);

			idxs_j[k] = 3;
			double p3 = _iterate_interpolate(k+1,d,frac_abscissas,idxs_j,p);

			return f1d_deriv_x_by_ref(frac_abscissas[k],p0,p1,p2,p3);

		} else { // Deeper

			idxs_j[delta] = 0;
			double dp0dxk = _iterate_deriv_x(delta+1,k,d,frac_abscissas,idxs_j,p);
			double dfdp0 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],0);

			idxs_j[delta] = 1;
			double dp1dxk = _iterate_deriv_x(delta+1,k,d,frac_abscissas,idxs_j,p);
			double dfdp1 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],1);

			idxs_j[delta] = 2;
			double dp2dxk = _iterate_deriv_x(delta+1,k,d,frac_abscissas,idxs_j,p);
			double dfdp2 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],2);

			idxs_j[delta] = 3;
			double dp3dxk = _iterate_deriv_x(delta+1,k,d,frac_abscissas,idxs_j,p);
			double dfdp3 = f1d_deriv_pt_value_by_ref(frac_abscissas[delta],3);

			return dfdp0 * dp0dxk + dfdp1 * dp1dxk + dfdp2 * dp2dxk + dfdp3 * dp3dxk;
		};
	};

	double Grid::Impl::get_deriv_wrt_x(std::vector<double> abscissas, int k) {
		return get_deriv_wrt_x_by_ref(abscissas,k);
	};

	double Grid::Impl::get_deriv_wrt_x_by_ref(const std::vector<double> &abscissas, int k) {
		// Stopping
		int d = get_no_dims();

		// Idxs
		IdxSet4 idxs_j(d);

		// nbrs p
		Nbr4 nbr4 = get_surrounding_4_grid_pts_by_ref(abscissas);

		// Frac
		std::vector<double> frac_abscissas = nbr4.frac_abscissas;

		// Convert nbr to p
		P4 p(nbr4);

		// Iterate
		return _iterate_deriv_x(0,k,d,frac_abscissas,idxs_j,p);
	};

	/********************
	1D
	********************/

	double Grid::Impl::f1d_interpolate(double x_frac, double p0, double p1, double p2, double p3) const {
		return f1d_interpolate_by_ref(x_frac,p0,p1,p2,p3);
	};
	double Grid::Impl::f1d_interpolate_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const {
		return (-0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3)*pow(x_frac,3) + (p0 - 2.5*p1 + 2.0*p2 - 0.5*p3)*pow(x_frac,2) + (-0.5*p0 + 0.5*p2)*x_frac + p1;
	};

	// p_idx = 0,1,2,3, depending on loc
	double Grid::Impl::f1d_deriv_pt_value(double x_frac, int p_idx) const {
		return f1d_deriv_pt_value_by_ref(x_frac,p_idx);
	};
	double Grid::Impl::f1d_deriv_pt_value_by_ref(const double &x_frac, int p_idx) const {
		if (p_idx==0) {
			return -0.5*pow(x_frac,3) + pow(x_frac,2) - 0.5*x_frac;
		} else if (p_idx==1) {
			return 1.5*pow(x_frac,3) - 2.5*pow(x_frac,2) + 1.0;
		} else if (p_idx==2) {
			return -1.5*pow(x_frac,3) + 2.0*pow(x_frac,2) + 0.5*x_frac;
		} else if (p_idx==3) {
			return 0.5*pow(x_frac,3) - 0.5*pow(x_frac,2);
		};

		// Never get here
		std::cerr << ">>> Error: Grid::Impl::f1d_deriv_pt_value_by_ref <<< p = 0,1,2,3" << std::endl;
		exit(EXIT_FAILURE);
	};

	double Grid::Impl::f1d_deriv_x(double x_frac, double p0, double p1, double p2, double p3) const {
		return f1d_deriv_x_by_ref(x_frac,p0,p1,p2,p3);
	};
	double Grid::Impl::f1d_deriv_x_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const {
		return 3.0*(-0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3)*pow(x_frac,2) + 2.0*(p0 - 2.5*p1 + 2.0*p2 - 0.5*p3)*x_frac + (-0.5*p0 + 0.5*p2);
	};

    /********************
	Apply M or P mappings
	********************/

    IdxSet Grid::Impl::apply_m_mapping_by_ref(const IdxSet &idxs) const {
		IdxSet idxs_m = idxs;
		for (auto dim=0; dim<get_no_dims(); dim++) {
			if (idxs_m[dim] == -1) {
				idxs_m[dim] += 1;
			} else if (idxs_m[dim] == _dims[dim].get_no_pts()) {
				idxs_m[dim] -= 1;
			};
		};
		return idxs_m;
    };
    IdxSet Grid::Impl::apply_p_mapping_by_ref(const IdxSet &idxs) const {
		IdxSet idxs_p = idxs;
		for (auto dim=0; dim<get_no_dims(); dim++) {
			if (idxs_p[dim] == -1) {
				idxs_p[dim] += 2;
			} else if (idxs_p[dim] == _dims[dim].get_no_pts()) {
				idxs_p[dim] -= 2;
			};
		};
		return idxs_p;
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
	Print
	********************/

	void Grid::print_grid_pts_inside() const {
		_impl->print_grid_pts_inside();
	};

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

	/*****
	Const ptrs
	*****/

	const GridPt* Grid::get_grid_point(std::vector<int> grid_idxs) const {
		return get_grid_point(IdxSet(grid_idxs));
	};
	const GridPt* Grid::get_grid_point(IdxSet idx_set) const {
		return _impl->get_grid_point(idx_set);
	};

	const GridPtOut* Grid::get_grid_point_outside(std::vector<int> grid_idxs) const {
		return get_grid_point_outside(IdxSet(grid_idxs));
	};
	const GridPtOut* Grid::get_grid_point_outside(IdxSet idx_set) const {
		return _impl->get_grid_point_outside(idx_set);
	};

	/*****
	Refs
	*****/

	GridPt& Grid::get_grid_point_ref(std::vector<int> grid_idxs) {
		return get_grid_point_ref(IdxSet(grid_idxs));
	};
	GridPt& Grid::get_grid_point_ref(IdxSet idx_set) {
		return _impl->get_grid_point_ref(idx_set);
	};

	GridPtOut& Grid::get_grid_point_outside_ref(std::vector<int> grid_idxs) {
		return get_grid_point_outside_ref(IdxSet(grid_idxs));
	};
	GridPtOut& Grid::get_grid_point_outside_ref(IdxSet idx_set) {
		return _impl->get_grid_point_outside_ref(idx_set);
	};

	/********************
	Get grid points surrounding a point
	********************/

	Nbr2 Grid::get_surrounding_2_grid_pts(std::vector<double> abscissas) const {
		return _impl->get_surrounding_2_grid_pts(abscissas);
	};
	Nbr4 Grid::get_surrounding_4_grid_pts(std::vector<double> abscissas) const {
		return _impl->get_surrounding_4_grid_pts(abscissas);
	};

	Nbr2 Grid::get_surrounding_2_grid_pts_by_ref(const std::vector<double>& abscissas) const {
		return _impl->get_surrounding_2_grid_pts_by_ref(abscissas);
	};
	Nbr4 Grid::get_surrounding_4_grid_pts_by_ref(const std::vector<double>& abscissas) const {
		return _impl->get_surrounding_4_grid_pts_by_ref(abscissas);
	};

	/********************
	Get a point by interpolating
	********************/

	double Grid::get_val(std::vector<double> abscissas) const {
		return _impl->get_val(abscissas);
	};
	double Grid::get_val_by_ref(const std::vector<double>& abscissas) const {
		return _impl->get_val_by_ref(abscissas);
	};

	/********************
	Get derivative
	********************/

	double Grid::get_deriv_wrt_pt_value(std::vector<double> abscissas, IdxSet4 idx_set) {
		return _impl->get_deriv_wrt_pt_value(abscissas,idx_set);
	};
	double Grid::get_deriv_wrt_pt_value_by_ref(const std::vector<double>& abscissas, const IdxSet4& idxs_k) {
		return _impl->get_deriv_wrt_pt_value_by_ref(abscissas,idxs_k);
	};

	/********************
	Get derivative wrt x
	********************/

	double Grid::get_deriv_wrt_x(std::vector<double> abscissas, int k) {
		return _impl->get_deriv_wrt_x(abscissas,k);
	};
	double Grid::get_deriv_wrt_x_by_ref(const std::vector<double>& abscissas, int k) {
		return _impl->get_deriv_wrt_x_by_ref(abscissas,k);
	};

	/********************
	1D funcs
	********************/

	double Grid::f1d_interpolate(double x_frac, double p0, double p1, double p2, double p3) const {
		return _impl->f1d_interpolate(x_frac,p0,p1,p2,p3);
	};
	double Grid::f1d_interpolate_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const {
		return _impl->f1d_interpolate_by_ref(x_frac,p0,p1,p2,p3);
	};

	double Grid::f1d_deriv_pt_value(double x_frac, int p_idx) const {
		return _impl->f1d_deriv_pt_value(x_frac,p_idx);
	};
	double Grid::f1d_deriv_pt_value_by_ref(const double &x_frac, int p_idx) const {
		return _impl->f1d_deriv_pt_value_by_ref(x_frac,p_idx);
	};
	
	double Grid::f1d_deriv_x(double x_frac, double p0, double p1, double p2, double p3) const {
		return _impl->f1d_deriv_x(x_frac,p0,p1,p2,p3);
	};
	double Grid::f1d_deriv_x_by_ref(const double &x_frac, const double &p0, const double &p1, const double &p2, const double &p3) const {
		return _impl->f1d_deriv_x_by_ref(x_frac,p0,p1,p2,p3);
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

    /********************
	Apply M or P mappings
	********************/

    IdxSet Grid::apply_m_mapping_by_ref(const IdxSet &idxs) const {
    	return apply_m_mapping_by_ref(idxs);
    };
    IdxSet Grid::apply_p_mapping_by_ref(const IdxSet &idxs) const {
    	return apply_p_mapping_by_ref(idxs);
    };


};