#include "../include/dcubic_bits/grid.hpp"

// Other headers
#include "../include/dcubic_bits/dimension_1d.hpp"
#include "../include/dcubic_bits/grid_pt.hpp"
#include "../include/dcubic_bits/neighbors.hpp" // includes idxset

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#define DIAG_PROJ 0

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Gridementation
	****************************************/

	Grid::Grid(int no_dims, Dimension1D** dims) {
		_no_dims = no_dims;
		_dims = new Dimension1D*[_no_dims];
		std::copy(dims,dims+_no_dims,_dims);

		_shared_constructor();
	};
	Grid::Grid(std::vector<Dimension1D*> dims) {
		_no_dims = dims.size();
		_dims = new Dimension1D*[_no_dims];
		for (auto i=0; i<_no_dims; i++) {
			_dims[i] = dims[i];
		};

		_shared_constructor();
	};
	Grid::Grid(const Grid& other) {
		_copy(other);
	};
	Grid::Grid(Grid&& other) {
		_move(other);
	};
    Grid& Grid::operator=(const Grid& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Grid& Grid::operator=(Grid&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	Grid::~Grid()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void Grid::_clean_up()
	{
		if (_dims) {
			for (auto i=0; i<_no_dims; i++) {
				// Grid does not own the Dimensions! Just delete the array
				_dims[i] = nullptr;
			};
			delete[] _dims;
			_dims = nullptr;
		};
		if (_no_pts_in_dim) {
			delete[] _no_pts_in_dim;
			_no_pts_in_dim = nullptr;
		};
		if (_grid_pts_in) {
			for (auto ipt=0; ipt<_no_grid_pts; ipt++) {
				if (_grid_pts_in[ipt]) {
					delete _grid_pts_in[ipt];
					_grid_pts_in[ipt] = nullptr;
				};
			};
			delete[] _grid_pts_in;
			_grid_pts_in = nullptr;
		};
		if (_grid_pts_out) {
			for (auto ipt=0; ipt<_no_grid_pts; ipt++) {
				if (_grid_pts_out[ipt]) {
					delete _grid_pts_out[ipt];
					_grid_pts_out[ipt] = nullptr;
				};
			};
			delete[] _grid_pts_out;
			_grid_pts_out = nullptr;
		};
	};
	void Grid::_copy(const Grid& other)
	{
		_no_dims = other._no_dims;
		_no_grid_pts = other._no_grid_pts;

		_dims = new Dimension1D*[_no_dims];
		std::copy(other._dims,other._dims+_no_dims,_dims);

		_no_pts_in_dim = new int[_no_dims];
		std::copy(other._no_pts_in_dim,other._no_pts_in_dim+_no_dims,_no_pts_in_dim);

		_grid_pts_in = new GridPtIn*[_no_grid_pts];
		std::copy(other._grid_pts_in,other._grid_pts_in+_no_grid_pts,_grid_pts_in);
		_grid_pts_out = new GridPtOut*[_no_grid_pts];
		std::copy(other._grid_pts_out,other._grid_pts_out+_no_grid_pts,_grid_pts_out);
	};
	void Grid::_move(Grid& other)
	{
		_no_dims = other._no_dims;
		_no_grid_pts = other._no_grid_pts;

		_dims = other._dims;

		_no_pts_in_dim = other._no_pts_in_dim;

		_grid_pts_in = other._grid_pts_in;
		_grid_pts_out = other._grid_pts_out;

		// Reset other
		other._no_dims = 0;
		other._no_grid_pts = 0;
		other._dims = nullptr;
		other._no_pts_in_dim = nullptr;
		other._grid_pts_in = nullptr;
		other._grid_pts_out = nullptr;
	};
	void Grid::_shared_constructor() {
		_no_grid_pts = 1;
		_no_pts_in_dim = new int[_no_dims];	
		for (auto dim=0; dim<_no_dims; dim++) {
			_no_pts_in_dim[dim] = _dims[dim]->get_no_pts();
			_no_grid_pts *= _no_pts_in_dim[dim] + 2; // outside left/right
		};

		// Allocate grid pts
		_grid_pts_in = new GridPtIn*[_no_grid_pts];
		std::fill_n(_grid_pts_in,_no_grid_pts,nullptr);
		_grid_pts_out = new GridPtOut*[_no_grid_pts];
		std::fill_n(_grid_pts_out,_no_grid_pts,nullptr);

		// Make grid pts inside
		double* abscissas = new double[_no_dims];
		_iterate_make_grid_pt_inside(0,IdxSet(_no_dims),abscissas);

		// Make grid pts outside
		Loc* locs = new Loc[_no_dims];
		_iterate_make_grid_pt_outside(0,IdxSet(_no_dims),abscissas,locs);

		// Clean up
		delete[] abscissas;
		delete[] locs;
		abscissas = nullptr;
		locs = nullptr;
	};

	/********************
	Make grid pts
	********************/

	void Grid::_iterate_make_grid_pt_inside(int dim, IdxSet idxs, double* abscissas) {
		if (dim != _no_dims) {
			// Deeper!
			for (idxs[dim]=1; idxs[dim]<=_no_pts_in_dim[dim]; idxs[dim]++) { // extra 2 for left/right
				_iterate_make_grid_pt_inside(dim+1, idxs, abscissas);
			};
		} else {
			// Do something

			// Make the abscissa
			for (auto dim2=0; dim2<_no_dims; dim2++) {
				abscissas[dim2] = _dims[dim2]->get_pt_by_idx(idxs[dim2],true);
			};

			// Make the grid pt
			move_grid_point_inside(idxs,new GridPtIn(_no_dims,abscissas));
		};
	};

	void Grid::_iterate_make_grid_pt_outside(int dim, IdxSet idxs, double* abscissas, Loc* locs) {
		if (dim != _no_dims) {
			// Deeper!
			for (idxs[dim]=0; idxs[dim]<=_no_pts_in_dim[dim]+1; idxs[dim]++) { // extra 2 for left/right
				_iterate_make_grid_pt_outside(dim+1, idxs, abscissas, locs);
			};
		} else {
			// Do something

			// Check that at least one dim is outside
			bool inside=true;
			IdxSet idxs_p1(idxs), idxs_p2(idxs);
			for (auto dim2=0; dim2<_no_dims; dim2++) {
				if (idxs[dim2] == 0) {
					inside = false;
					locs[dim2] = Loc::OUTSIDE_LEFT;
					idxs_p1[dim2] += 1;
					idxs_p2[dim2] += 2;
				} else if (idxs[dim2] == _no_pts_in_dim[dim2]+1) {
					inside = false;
					locs[dim2] = Loc::OUTSIDE_RIGHT;
					idxs_p1[dim2] -= 1;
					idxs_p2[dim2] -= 2;
				} else {
					locs[dim2] = Loc::INSIDE;
				};	
			};

			// This pt is inside
			if (inside) {
				return;
			};

			// Make the abscissa
			for (auto dim2=0; dim2<_no_dims; dim2++) {
				abscissas[dim2] = _dims[dim2]->get_pt_by_idx(idxs[dim2],true);
			};

			// Make the grid pt
			move_grid_point_outside(idxs,new GridPtOut(_no_dims,abscissas,get_grid_point_inside(idxs_p1),get_grid_point_inside(idxs_p2),locs));
		};
	};

	/********************
	Print
	********************/

	void Grid::print() const {
		for (auto ipt=0; ipt<_no_grid_pts; ipt++) {
			std::cout << _grid_pts_in[ipt]->print() << std::endl;
		};
	};

	/********************
	Get dims
	********************/

	int Grid::get_no_dims() const {
		return _no_dims;
	};
	Dimension1D* Grid::get_dim(int dim) const {
		return _dims[dim];
	};

	/********************
	Get grid point
	********************/

	GridPt* Grid::get_grid_point(IdxSet idxs) const {
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= _no_pts_in_dim[dim2];
			};

			idx += term;
		};
		if (_grid_pts_in[idx]) {
			return _grid_pts_in[idx];
		} else if (_grid_pts_out[idx]) {
			return _grid_pts_out[idx];
		} else {
			std::cerr << ">>> Error: get_grid_point <<< does not exist at idxs: " << idxs << std::endl;
			exit(EXIT_FAILURE);
		};
	};
	GridPtIn* Grid::get_grid_point_inside(IdxSet idxs) const {
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= _no_pts_in_dim[dim2];
			};

			idx += term;
		};
		if (_grid_pts_in[idx]) {
			return _grid_pts_in[idx];
		} else {
			std::cerr << ">>> Error: get_grid_point_inside <<< does not exist at idxs: " << idxs << std::endl;
			exit(EXIT_FAILURE);
		};
	};
	GridPtOut* Grid::get_grid_point_outside(IdxSet idxs) const {
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= _no_pts_in_dim[dim2];
			};

			idx += term;
		};
		if (_grid_pts_out[idx]) {
			return _grid_pts_out[idx];
		} else {
			std::cerr << ">>> Error: get_grid_point_outside <<< does not exist at idxs: " << idxs << std::endl;
			exit(EXIT_FAILURE);
		};		
	};
	void Grid::move_grid_point_inside(IdxSet idxs, GridPtIn* grid_pt) {
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= _no_pts_in_dim[dim2];
			};

			idx += term;
		};

		// Delete existing grid pt
		if (_grid_pts_in[idx]) {
			delete _grid_pts_in[idx];
		};

		_grid_pts_in[idx] = grid_pt;
	};
	void Grid::move_grid_point_outside(IdxSet idxs, GridPtOut* grid_pt) {
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= _no_pts_in_dim[dim2];
			};

			idx += term;
		};

		// Delete existing grid pt
		if (_grid_pts_out[idx]) {
			delete _grid_pts_out[idx];
		};

		_grid_pts_out[idx] = grid_pt;
	};

	/********************
	Get grid points surrounding a point
	********************/

	Nbr2 Grid::get_surrounding_2_grid_pts(double* abscissas) const {

		// Frac abscissas
		double* frac_abscissas = new double[_no_dims];

		// Get bounding idxs
		IdxSet idxs_lower(_no_dims), idxs_upper(_no_dims);
		for (auto dim=0; dim<_no_dims; dim++) {
			// Check in dim
			if (!_dims[dim]->check_if_pt_is_inside_domain(abscissas[dim])) {
				// Outside grid
				std::cerr << ">>> Error: Grid::get_surrounding_2_grid_pts <<< Abscissa in dim: " << dim << " value: " << abscissas[dim] << " is outside the grid: " << _dims[dim]->get_min_pt() << " to: " << _dims[dim]->get_max_pt() << std::endl;
				exit(EXIT_FAILURE);
			};

			idxs_lower[dim] = _dims[dim]->get_idxs_surrounding_pt(abscissas[dim],true);
			idxs_upper[dim] = idxs_lower[dim]+1;

			// Frac
			frac_abscissas[dim] = _dims[dim]->get_frac_between(abscissas[dim],idxs_lower[dim],true);
		};

		// Returned
		Nbr2 nbr2(idxs_lower, frac_abscissas);

		// Clean up
		delete[] frac_abscissas;
		frac_abscissas = nullptr;

		// Iterate to fill out the map
		_iterate_get_surrounding_2_grid_pts(0,IdxSet(_no_dims),idxs_lower,nbr2);

		return nbr2;
	};

	void Grid::_iterate_get_surrounding_2_grid_pts(int dim, IdxSet idxs_local, const IdxSet &idxs_lower, Nbr2 &nbr2) const {
		if (dim != _no_dims) {
			// Deeper!
			// Can be lower (=0) or higher (=+1) in this dim
			idxs_local[dim] = 0;
			_iterate_get_surrounding_2_grid_pts(dim+1,idxs_local,idxs_lower,nbr2);
			idxs_local[dim] = 1;
			_iterate_get_surrounding_2_grid_pts(dim+1,idxs_local,idxs_lower,nbr2);

		} else {
			// Do something

			// Add to nbr2
			nbr2.set_grid_point(idxs_local, get_grid_point(idxs_lower+idxs_local));
		};
	};

	Nbr4 Grid::get_surrounding_4_grid_pts(double* abscissas) const {

		// Frac abscissas
		double* frac_abscissas = new double[_no_dims];

		// Get bounding idxs
		IdxSet idxs_p0(_no_dims), idxs_p1(_no_dims), idxs_p2(_no_dims), idxs_p3(_no_dims);
		for (auto dim=0; dim<_no_dims; dim++) {
			// Check in dim
			if (!_dims[dim]->check_if_pt_is_inside_domain(abscissas[dim])) {
				// Outside grid
				std::cerr << ">>> Error: Grid::get_surrounding_4_grid_pts <<< Abscissa in dim: " << dim << " value: " << abscissas[dim] << " is outside the grid: " << _dims[dim]->get_min_pt() << " to: " << _dims[dim]->get_max_pt() << std::endl;
				exit(EXIT_FAILURE);
			};

			idxs_p1[dim] = _dims[dim]->get_idxs_surrounding_pt(abscissas[dim],true);
			idxs_p0[dim] = idxs_p1[dim]-1;
			idxs_p2[dim] = idxs_p1[dim]+1;
			idxs_p3[dim] = idxs_p1[dim]+2;

			// Frac
			frac_abscissas[dim] = _dims[dim]->get_frac_between(abscissas[dim],idxs_p1[dim],true);
		};

		// Returned
		Nbr4 nbr4(idxs_p1, frac_abscissas);

		// Clean up
		delete[] frac_abscissas;
		frac_abscissas = nullptr;

		// Iterate to fill out the map
		_iterate_get_surrounding_4_grid_pts(0,IdxSet(_no_dims),idxs_p0,nbr4);

		return nbr4;
	};

	void Grid::_iterate_get_surrounding_4_grid_pts(int dim, IdxSet idxs_local, const IdxSet &idxs_p0, Nbr4 &nbr4) const {
		if (dim != _no_dims) {
			// Deeper!
			// Can be 0,1,2,3
			idxs_local[dim] = 0;
			_iterate_get_surrounding_4_grid_pts(dim+1,idxs_local,idxs_p0,nbr4);
			idxs_local[dim] = 1;
			_iterate_get_surrounding_4_grid_pts(dim+1,idxs_local,idxs_p0,nbr4);
			idxs_local[dim] = 2;
			_iterate_get_surrounding_4_grid_pts(dim+1,idxs_local,idxs_p0,nbr4);
			idxs_local[dim] = 3;
			_iterate_get_surrounding_4_grid_pts(dim+1,idxs_local,idxs_p0,nbr4);

		} else {
			// Do something

			// Add to nbr4
			nbr4.set_grid_point(idxs_local, get_grid_point(idxs_p0+idxs_local));
		};
	};

	/********************
	Get a point by interpolating
	********************/

	double Grid::get_val(double* abscissas) const {
		// nbrs p
		Nbr4 nbr4 = get_surrounding_4_grid_pts(abscissas);

		// Return
		return get_val(nbr4);
	};
	double Grid::get_val(const Nbr4 &nbr4) const {
		// Iterate
		return _iterate_interpolate(0,_no_dims,nbr4,IdxSet(_no_dims));
	};
	double Grid::_iterate_interpolate(int delta, int d, const Nbr4 &nbr4, IdxSet idxs_j) const {
		double p0,p1,p2,p3;

		if (delta == d) {
			
			// Arrived; return the point P itself
			return nbr4.get_grid_point(idxs_j)->get_ordinate();

		} else {

			// Deeper
			idxs_j[delta] = 0;
			p0 = _iterate_interpolate(delta+1, d, nbr4, idxs_j);

			idxs_j[delta] = 1;
			p1 = _iterate_interpolate(delta+1, d, nbr4, idxs_j);

			idxs_j[delta] = 2;
			p2 = _iterate_interpolate(delta+1, d, nbr4, idxs_j);

			idxs_j[delta] = 3;
			p3 = _iterate_interpolate(delta+1, d, nbr4, idxs_j);

			return f1d_interpolate(nbr4.get_frac_abscissa(delta),p0,p1,p2,p3);
		};
	};

	/********************
	Get derivative wrt pt value
	********************/

	double Grid::get_deriv_wrt_pt_value(double* abscissas, IdxSet idxs_k) {

		// nbrs
		Nbr4 nbr4 = get_surrounding_4_grid_pts(abscissas);

		// Get
		return get_deriv_wrt_pt_value(nbr4,idxs_k);
	};
	double Grid::get_deriv_wrt_pt_value(const Nbr4 &nbr4, IdxSet idxs_k) {

		// Check: are there any exterior grid points?
		if (nbr4.check_are_all_pts_inside()) {
			// Case 1: Totally interior point; no dimension near boundary

			double ret=1.0;
			for (auto alpha=0; alpha<_no_dims; alpha++) {
				ret *= f1d_deriv_pt_value(nbr4.get_frac_abscissa(alpha), idxs_k[alpha]);
			};
			return ret;

		} else {
			// Case 2: At least one dimension near boundary

			// Check that we are not taking an illegal derivative with respect to an outside grid pt
			if (nbr4.get_grid_point(idxs_k)) {
				if (nbr4.get_grid_point(idxs_k)->get_type() == GridPtType::OUTSIDE) {
					std::cerr << ">>> Error: Grid::get_deriv_wrt_pt_value <<< The pt specified to take a derivative with respect to is not a real point (an exterior point that is estimated by a linear approx) - this is not allowed!" << std::endl;
					exit(EXIT_FAILURE);
				};
			} else {
				std::cerr << ">>> Error: Grid::get_deriv_wrt_pt_value <<< pt at this idx does not exist" << std::endl;
				exit(EXIT_FAILURE);
			};

			// Iterate
			return _iterate_deriv_pt_value(0, _no_dims, nbr4, IdxSet(_no_dims), idxs_k);
		};
	};

	double Grid::_iterate_deriv_pt_value(int delta, int d, const Nbr4 &nbr4, IdxSet idxs_j, IdxSet idxs_k) const {

		if (delta == d) {
			// Done; evaluate
			// Check main condition
			bool at_least_one_cond_met = false;
			for (auto dim=0; dim<d; dim++) {
				if ((nbr4.get_idx_i(dim) == 1 && idxs_j[dim] == 0) || (nbr4.get_idx_i(dim) == _no_pts_in_dim[dim]-1 && idxs_j[dim] == 3)) {
					at_least_one_cond_met = true;
					break;
				};
			};

			if (at_least_one_cond_met) {

				// Yes, at least one
			

				// Form the new idxs
				IdxSet idxs_j_global = idxs_j - 1 + nbr4.get_idxs_i();
				IdxSet idxs_k_global = idxs_k - 1 + nbr4.get_idxs_i();

				// Apply M mapping
				IdxSet idxs_m = apply_m_mapping(idxs_j_global);

				// Check
				if (idxs_m == idxs_k_global) {
					return 2.0;
				};

				// Apply P mapping
				IdxSet idxs_p = apply_p_mapping(idxs_j_global);

				// Check
				if (idxs_p == idxs_k_global) {
					return -1.0;
				};

				// Failed (but this is ok!)
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
			double dp0dp = _iterate_deriv_pt_value(delta+1,d,nbr4,idxs_j,idxs_k);
			double dfdp0 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),0);

			idxs_j[delta] = 1;
			double dp1dp = _iterate_deriv_pt_value(delta+1,d,nbr4,idxs_j,idxs_k);
			double dfdp1 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),1);

			idxs_j[delta] = 2;
			double dp2dp = _iterate_deriv_pt_value(delta+1,d,nbr4,idxs_j,idxs_k);
			double dfdp2 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),2);

			idxs_j[delta] = 3;
			double dp3dp = _iterate_deriv_pt_value(delta+1,d,nbr4,idxs_j,idxs_k);
			double dfdp3 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),3);

			return dfdp0 * dp0dp + dfdp1 * dp1dp + dfdp2 * dp2dp + dfdp3 * dp3dp;
		};
	};

	/********************
	Get derivative wrt x
	********************/

	double Grid::get_deriv_wrt_x(double* abscissas, int k) {

		// nbrs p
		Nbr4 nbr4 = get_surrounding_4_grid_pts(abscissas);

		// Get
		return get_deriv_wrt_x(nbr4, k);
	};
	double Grid::get_deriv_wrt_x(const Nbr4 &nbr4, int k) {
		// Iterate
		return _iterate_deriv_x(0,k,_no_dims,nbr4,IdxSet(_no_dims));
	};
	double Grid::_iterate_deriv_x(int delta, int k, int d, const Nbr4 &nbr4, IdxSet idxs_j) const {
		
		if (delta == k) { // Evaluate the derivative

			idxs_j[k] = 0;
			double p0 = _iterate_interpolate(k+1,d,nbr4,idxs_j);

			idxs_j[k] = 1;
			double p1 = _iterate_interpolate(k+1,d,nbr4,idxs_j);

			idxs_j[k] = 2;
			double p2 = _iterate_interpolate(k+1,d,nbr4,idxs_j);

			idxs_j[k] = 3;
			double p3 = _iterate_interpolate(k+1,d,nbr4,idxs_j);

			return f1d_deriv_x(nbr4.get_frac_abscissa(k),p0,p1,p2,p3);

		} else { // Deeper

			idxs_j[delta] = 0;
			double dp0dxk = _iterate_deriv_x(delta+1,k,d,nbr4,idxs_j);
			double dfdp0 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),0);

			idxs_j[delta] = 1;
			double dp1dxk = _iterate_deriv_x(delta+1,k,d,nbr4,idxs_j);
			double dfdp1 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),1);

			idxs_j[delta] = 2;
			double dp2dxk = _iterate_deriv_x(delta+1,k,d,nbr4,idxs_j);
			double dfdp2 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),2);

			idxs_j[delta] = 3;
			double dp3dxk = _iterate_deriv_x(delta+1,k,d,nbr4,idxs_j);
			double dfdp3 = f1d_deriv_pt_value(nbr4.get_frac_abscissa(delta),3);

			return dfdp0 * dp0dxk + dfdp1 * dp1dxk + dfdp2 * dp2dxk + dfdp3 * dp3dxk;
		};
	};

	/********************
	1D
	********************/

	double Grid::f1d_interpolate(double x_frac, double p0, double p1, double p2, double p3) const {
		return (-0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3)*pow(x_frac,3) + (p0 - 2.5*p1 + 2.0*p2 - 0.5*p3)*pow(x_frac,2) + (-0.5*p0 + 0.5*p2)*x_frac + p1;
	};

	// p_idx = 0,1,2,3, depending on loc
	double Grid::f1d_deriv_pt_value(double x_frac, int p_idx) const {
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
		std::cerr << ">>> Error: Grid::f1d_deriv_pt_value <<< p =/= 0,1,2,3" << std::endl;
		exit(EXIT_FAILURE);
	};

	double Grid::f1d_deriv_x(double x_frac, double p0, double p1, double p2, double p3) const {
		return 3.0*(-0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3)*pow(x_frac,2) + 2.0*(p0 - 2.5*p1 + 2.0*p2 - 0.5*p3)*x_frac + (-0.5*p0 + 0.5*p2);
	};

    /********************
	Apply M or P mappings
	********************/

    IdxSet Grid::apply_m_mapping(IdxSet idxs) const {
		IdxSet idxs_m = idxs;
		for (auto dim=0; dim<_no_dims; dim++) {
			if (idxs_m[dim] == 0) {
				idxs_m[dim] += 1;
			} else if (idxs_m[dim] == _no_pts_in_dim[dim]+1) {
				idxs_m[dim] -= 1;
			};
		};
		return idxs_m;
    };
    IdxSet Grid::apply_p_mapping(IdxSet idxs) const {
		IdxSet idxs_p = idxs;
		for (auto dim=0; dim<_no_dims; dim++) {
			if (idxs_p[dim] == 0) {
				idxs_p[dim] += 2;
			} else if (idxs_p[dim] == _no_pts_in_dim[dim]+1) {
				idxs_p[dim] -= 2;
			};
		};
		return idxs_p;
    };

	/********************
	Read/write grid
	********************/

	void Grid::read_from_file(std::string fname) {
		std::ifstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Grid::read_from_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Abscissas
		std::string abscissa_str;
		double abscissa;
		IdxSet idx_set(_no_dims);

		// Ordinate
		std::string ordinate_str="";
		double ordinate;

		// Read
		std::string line;
		std::istringstream iss;
		while (getline(f,line)) {
			// Skip empty lines
			if (line == "") { continue; };
			// Line
			iss = std::istringstream(line);
			// Abscissas
			for (auto dim=0; dim<_no_dims; dim++) {
				iss >> abscissa_str;
				abscissa = atof(abscissa_str.c_str());
				idx_set[dim] = _dims[dim]->get_closest_idx(abscissa,true);
			};
			// Ordinate
			iss >> ordinate_str;
			ordinate = atof(ordinate_str.c_str());
			// Set grid pt
			get_grid_point_inside(idx_set)->set_ordinate(ordinate);
			// Reset strs
			abscissa_str = "";
			ordinate_str = "";
		};

		// Close
		f.close();
	};
	void Grid::write_to_file(std::string fname) const {
		std::ofstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			std::cerr << ">>> Error: Grid::write_to_file <<< could not write to file: " << fname << std::endl;
			exit(EXIT_FAILURE);
		};

		// Go through all grid pts
		bool first_line=true;
		for (auto i=0; i<_no_grid_pts; i++) {

			if (!_grid_pts_in[i]) {
				continue;
			};

			if (!first_line) {
				f << "\n";
			} else {
				first_line = false;
			};

			// Write
			for (auto dim=0; dim<_no_dims; dim++) {
				f << _grid_pts_in[i]->get_abscissa(dim) << " ";
			};
			f << _grid_pts_in[i]->get_ordinate();
		};

		// Close
		f.close();
	};
};