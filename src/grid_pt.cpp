#include "../include/dcubic_bits/grid_pt.hpp"

#include <iostream>
#include <sstream>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Grid pt
	****************************************/

	GridPt::GridPt(int no_dims, double* abscissas) {
		_no_dims = no_dims;
		_abcissas = new double[_no_dims];
		std::fill_n(_abcissas,_no_dims,0.0);
	};
	GridPt::GridPt(std::vector<double> abscissas) {
		_no_dims = abscissas.size();
		_abcissas = new double[_no_dims];
		for (auto i=0; i<_no_dims; i++) {
			_abcissas[i] = abscissas[i];
		};	
	};
	GridPt::GridPt(const GridPt& other) {
		_copy(other);
	};
	GridPt::GridPt(GridPt&& other) {
		_move(other);
	};
    GridPt& GridPt::operator=(const GridPt& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    GridPt& GridPt::operator=(GridPt&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	GridPt::~GridPt()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void GridPt::_clean_up()
	{
		// Nothing...
	};
	void GridPt::_copy(const GridPt& other)
	{
		_no_dims = other._no_dims;
		_abcissas = new double[_no_dims];
		std::copy(other._abcissas,other._abcissas+_no_dims,_abcissas);
	};
	void GridPt::_move(GridPt& other)
	{
		_no_dims = other._no_dims;
		_abcissas = other._abcissas;

		// Reset other
		other._no_dims = 0;
		other._abcissas = nullptr;
	};

	/********************
	Access
	********************/

	// Abscissa
	double GridPt::get_abscissa(int dim) const {
		return _abcissas[dim];
	};

	/********************
	Print
	********************/

	std::string GridPt::print() const {
		std::ostringstream s;
		s << "(";
		for (auto dim=0; dim<_no_dims; dim++) {
			s << _abcissas[dim];
			if (dim != _no_dims-1) {
				s << " ";
			};
		};
		s << "): " << get_ordinate();
		return s.str();
	};




































	/****************************************
	Grid pt interior
	****************************************/

	GridPtIn::GridPtIn(int no_dims, double* abscissas) : GridPt(no_dims,abscissas) {
		_ordinate = 0.0;
	};
	GridPtIn::GridPtIn(std::vector<double> abscissas) : GridPt(abscissas) {
		_ordinate = 0.0;
	};
	GridPtIn::GridPtIn(const GridPtIn& other) : GridPt(other) {
		_copy(other);
	};
	GridPtIn::GridPtIn(GridPtIn&& other) : GridPt(std::move(other)) {
		_move(other);
	};
    GridPtIn& GridPtIn::operator=(const GridPtIn& other) {
		if (this != &other) {
	        GridPt::operator=(other);
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    GridPtIn& GridPtIn::operator=(GridPtIn&& other) {
		if (this != &other) {
	        GridPt::operator=(std::move(other));
			_clean_up();
			_move(other);
		};
		return *this;
    };
	GridPtIn::~GridPtIn()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void GridPtIn::_clean_up()
	{
		// Nothing...
	};
	void GridPtIn::_copy(const GridPtIn& other)
	{
		_ordinate = other._ordinate;
	};
	void GridPtIn::_move(GridPtIn& other)
	{
		_ordinate = other._ordinate;

		// Reset other
		other._ordinate = 0.0;
	};

	/********************
	Access
	********************/

	// Abscissa
	double GridPtIn::get_ordinate() const {
		return _ordinate;
	};
	const double& GridPtIn::get_ordinate_const_ref() const {
		return _ordinate;
	};
	void GridPtIn::set_ordinate(double val) {
		_ordinate = val;
	};
	void GridPtIn::increment_ordinate(double val) {
		_ordinate += val;
	};




































	/****************************************
	Grid pt exterior
	****************************************/

	GridPtOut::GridPtOut(int no_dims, double* abscissas, const GridPtIn* p1, const GridPtIn* p2, Loc* locs) : GridPt(no_dims,abscissas) {
		_p1 = p1;
		_p2 = p2;
		_locs = new Loc[_no_dims];
		std::copy(locs,locs+_no_dims,_locs);
	};
	GridPtOut::GridPtOut(std::vector<double> abscissas, const GridPtIn* p1, const GridPtIn* p2, std::vector<Loc> locs) : GridPt(abscissas) {
		_p1 = p1;
		_p2 = p2;
		_locs = new Loc[_no_dims];
		for (int i=0; i<_no_dims; i++) {
			_locs[i] = locs[i];
		};
	};
	GridPtOut::GridPtOut(const GridPtOut& other) : GridPt(other) {
		_copy(other);
	};
	GridPtOut::GridPtOut(GridPtOut&& other) : GridPt(std::move(other)) {
		_move(other);
	};
    GridPtOut& GridPtOut::operator=(const GridPtOut& other) {
		if (this != &other) {
	        GridPt::operator=(other);
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    GridPtOut& GridPtOut::operator=(GridPtOut&& other) {
		if (this != &other) {
	        GridPt::operator=(std::move(other));
			_clean_up();
			_move(other);
		};
		return *this;
    };
	GridPtOut::~GridPtOut()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void GridPtOut::_clean_up()
	{
		if (_locs) {
			delete _locs;
			_locs = nullptr;
		};
	};
	void GridPtOut::_copy(const GridPtOut& other)
	{
		_p1 = other._p1;
		_p2 = other._p2;

		_locs = new Loc[_no_dims];
		std::copy(other._locs,other._locs+_no_dims,_locs);
	};
	void GridPtOut::_move(GridPtOut& other)
	{
		_p1 = other._p1;
		_p2 = other._p2;
		_locs = other._locs;

		// Reset other
		other._p1 = nullptr;
		other._p2 = nullptr;
		other._locs = nullptr;
	};

	/********************
	Access
	********************/

	// Abscissa
	double GridPtOut::get_ordinate() const {
		return 2.0 * _p1->get_ordinate() - _p2->get_ordinate(); 
	};

	// Get dependent points
	const GridPtIn* GridPtOut::get_dep_p1() const {
		return _p1;
	};
	const GridPtIn* GridPtOut::get_dep_p2() const {
		return _p2;
	};

	// Which dims are outside
	Loc GridPtOut::get_loc_in_dim(int dim) const {
		if (dim >= _no_dims || dim < 0) {
			std::cout << ">>> Error: GridPtOut::get_loc_in_dim <<< requested dim: " << dim << " but there are " << _no_dims << " dims." << std::endl;
			exit(EXIT_FAILURE);
		};
		return _locs[dim];
	};
};