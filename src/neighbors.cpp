#include "../include/dcubic_bits/neighbors.hpp"

// Other headers
#include "../include/dcubic_bits/grid_pt.hpp"

#include <cmath>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Neighbor 2
	****************************************/

	/********************
	Constructor
	********************/

	Nbr2::Nbr2(IdxSet idxs_i) : _idxs_i(idxs_i) {
		_no_dims = _idxs_i.size();

		_frac_abscissas = new double[_no_dims];
		std::fill_n(_frac_abscissas,_no_dims,0.0);

		_shared_constructor();
	};
	Nbr2::Nbr2(IdxSet idxs_i, double* frac_abscissas) : _idxs_i(idxs_i) {
		_idxs_i = idxs_i;
		_no_dims = _idxs_i.size();

		_frac_abscissas = new double[_no_dims];
		std::copy(frac_abscissas,frac_abscissas+_no_dims,_frac_abscissas);

		_shared_constructor();
	};
	Nbr2::Nbr2(IdxSet idxs_i, std::vector<double> frac_abscissas) : _idxs_i(idxs_i) {
		_idxs_i = idxs_i;
		_no_dims = _idxs_i.size();

		_frac_abscissas = new double[_no_dims];
		for (auto i=0; i<_no_dims; i++) {
			_frac_abscissas[i] = frac_abscissas[i];
		};

		_shared_constructor();
	};
	Nbr2::Nbr2(const Nbr2& other) : _idxs_i(other._idxs_i) {
		_copy(other);
	};
	Nbr2::Nbr2(Nbr2&& other) : _idxs_i(std::move(other._idxs_i)) {
		_move(other);
	};
    Nbr2& Nbr2::operator=(const Nbr2& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Nbr2& Nbr2::operator=(Nbr2&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	Nbr2::~Nbr2()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void Nbr2::_clean_up()
	{
		if (_grid_pts) {
			for (auto i=0; i<_no_grid_pts; i++) {
				// Nbr4 doesnt own grid pt; just clear
				if (_grid_pts[i]) {
					_grid_pts[i] = nullptr;
				};
			};
			delete[] _grid_pts;
			_grid_pts = nullptr;
		};

		if (_frac_abscissas) {
			delete[] _frac_abscissas;
			_frac_abscissas = nullptr;
		};	
	};
	void Nbr2::_copy(const Nbr2& other)
	{
		_no_dims = other._no_dims;
		_idxs_i = other._idxs_i;
		_no_grid_pts = other._no_grid_pts;

		_grid_pts = new GridPt*[_no_grid_pts];
		std::copy(other._grid_pts,other._grid_pts+_no_grid_pts,_grid_pts);

		_frac_abscissas = new double[_no_dims];
		std::copy(other._frac_abscissas,other._frac_abscissas+_no_dims,_frac_abscissas);
	};
	void Nbr2::_move(Nbr2& other)
	{
		_no_dims = other._no_dims;
		_idxs_i = std::move(other._idxs_i);
		_no_grid_pts = other._no_grid_pts;

		_grid_pts = other._grid_pts;

		_frac_abscissas = other._frac_abscissas;

		// Reset other
		other._no_dims = 0;
		other._no_grid_pts = 0;
		other._grid_pts = nullptr;
		other._frac_abscissas = nullptr;
	};
	void Nbr2::_shared_constructor() {
		_no_grid_pts = pow(2,_no_dims);
		_grid_pts = new GridPt*[_no_grid_pts];
		std::fill_n(_grid_pts,_no_grid_pts,nullptr);
	};

	/********************
	Access
	********************/

	// Grid pts
	void Nbr2::set_grid_point(IdxSet idxs, GridPt *grid_pt) {
		int idx = 0;
		for (auto dim=0; dim<_no_dims; dim++) {
			idx += idxs[dim] * pow(2,dim);
		};
		_grid_pts[idx] = grid_pt;
	};
	GridPt* Nbr2::get_grid_point(IdxSet idxs) const {
		int idx = 0;
		for (auto dim=0; dim<_no_dims; dim++) {
			idx += idxs[dim] * pow(2,dim);
		};
		return _grid_pts[idx];
	};
	GridPt* Nbr2::get_grid_point(int i) const {
		return _grid_pts[i];
	};
	
	// No grid pts
	int Nbr2::get_no_grid_pts() const {
		return _no_grid_pts;
	};

	// Frac abscissas
	double Nbr2::get_frac_abscissa(int dim) const {
		return _frac_abscissas[dim];
	};

	// Center idx
	IdxSet Nbr2::get_idxs_i() const {
		return _idxs_i;
	};
	int Nbr2::get_idx_i(int dim) const {
		return _idxs_i[dim];
	};















































	/****************************************
	Neighbor 4
	****************************************/

	/********************
	Constructor
	********************/

	Nbr4::Nbr4(IdxSet idxs_i) : _idxs_i(idxs_i) {
		_no_dims = _idxs_i.size();

		_frac_abscissas = new double[_no_dims];
		std::fill_n(_frac_abscissas,_no_dims,0.0);

		_shared_constructor();
	};
	Nbr4::Nbr4(IdxSet idxs_i, double* frac_abscissas) : _idxs_i(idxs_i) {
		_idxs_i = idxs_i;
		_no_dims = _idxs_i.size();

		_frac_abscissas = new double[_no_dims];
		std::copy(frac_abscissas,frac_abscissas+_no_dims,_frac_abscissas);

		_shared_constructor();
	};
	Nbr4::Nbr4(IdxSet idxs_i, std::vector<double> frac_abscissas) : _idxs_i(idxs_i) {
		_idxs_i = idxs_i;
		_no_dims = _idxs_i.size();

		_frac_abscissas = new double[_no_dims];
		for (auto i=0; i<_no_dims; i++) {
			_frac_abscissas[i] = frac_abscissas[i];
		};

		_shared_constructor();
	};
	Nbr4::Nbr4(const Nbr4& other) : _idxs_i(other._idxs_i) {
		_copy(other);
	};
	Nbr4::Nbr4(Nbr4&& other) : _idxs_i(std::move(other._idxs_i)) {
		_move(other);
	};
    Nbr4& Nbr4::operator=(const Nbr4& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Nbr4& Nbr4::operator=(Nbr4&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	Nbr4::~Nbr4()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void Nbr4::_clean_up()
	{
		if (_grid_pts) {
			for (auto i=0; i<_no_grid_pts; i++) {
				// Nbr4 doesnt own grid pt; just clear
				if (_grid_pts[i]) {
					_grid_pts[i] = nullptr;
				};
			};
			delete[] _grid_pts;
			_grid_pts = nullptr;
		};

		if (_frac_abscissas) {
			delete[] _frac_abscissas;
			_frac_abscissas = nullptr;
		};
	};
	void Nbr4::_copy(const Nbr4& other)
	{
		_no_dims = other._no_dims;
		_idxs_i = other._idxs_i;
		_no_grid_pts = other._no_grid_pts;

		_grid_pts = new GridPt*[_no_grid_pts];
		std::copy(other._grid_pts,other._grid_pts+_no_grid_pts,_grid_pts);

		_frac_abscissas = new double[_no_dims];
		std::copy(other._frac_abscissas,other._frac_abscissas+_no_dims,_frac_abscissas);
	};
	void Nbr4::_move(Nbr4& other)
	{
		_no_dims = other._no_dims;
		_idxs_i = std::move(other._idxs_i);
		_no_grid_pts = other._no_grid_pts;

		_grid_pts = other._grid_pts;

		_frac_abscissas = other._frac_abscissas;

		// Reset other
		other._no_dims = 0;
		other._no_grid_pts = 0;
		other._grid_pts = nullptr;
		other._frac_abscissas = nullptr;
	};
	void Nbr4::_shared_constructor() {
		_no_grid_pts = pow(4,_no_dims);
		_grid_pts = new GridPt*[_no_grid_pts];
		std::fill_n(_grid_pts,_no_grid_pts,nullptr);
	};

	/********************
	Access
	********************/

	// Grid pts
	void Nbr4::set_grid_point(IdxSet idxs, GridPt *grid_pt) {
		int idx = 0;
		for (auto dim=0; dim<_no_dims; dim++) {
			idx += idxs[dim] * pow(4,dim);
		};
		_grid_pts[idx] = grid_pt;
	};
	GridPt* Nbr4::get_grid_point(IdxSet idxs) const {
		int idx = 0;
		for (auto dim=0; dim<_no_dims; dim++) {
			idx += idxs[dim] * pow(4,dim);
		};
		return _grid_pts[idx];
	};
	GridPt* Nbr4::get_grid_point(int i) const {
		return _grid_pts[i];
	};
	
	// No grid pts
	int Nbr4::get_no_grid_pts() const {
		return _no_grid_pts;
	};

	// Frac abscissas
	double Nbr4::get_frac_abscissa(int dim) const {
		return _frac_abscissas[dim];
	};

	// Center idx
	IdxSet Nbr4::get_idxs_i() const {
		return _idxs_i;
	};
	int Nbr4::get_idx_i(int dim) const {
		return _idxs_i[dim];
	};

	// Check if all pts are interior
	bool Nbr4::check_are_all_pts_inside() const {
		for (auto i=0; i<pow(4,_no_dims); i++) {
			if (_grid_pts[i]) {
				if (_grid_pts[i]->get_type() == GridPtType::OUTSIDE) {
					// At least one is outside
					return false;
				};
			};
		};
		// All inside
		return true;
	};
};