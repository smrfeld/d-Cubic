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

	GridPt::GridPt(std::vector<double> abscissas) {
		_abcissas = abscissas;
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
		_abcissas = other._abcissas;
	};
	void GridPt::_move(GridPt& other)
	{
		_abcissas = other._abcissas;

		// Reset other
		other._abcissas.clear();
	};

	/********************
	Access
	********************/

	// Abscissa
	double GridPt::get_abscissa(int dim) const {
		return _abcissas[dim];
	};
	std::vector<double> GridPt::get_abscissas() const {
		return _abcissas;
	};

	/********************
	Print
	********************/

	std::string GridPt::print() const {
		std::ostringstream s;
		s << "(";
		for (auto dim=0; dim<_abcissas.size(); dim++) {
			s << _abcissas[dim];
			if (dim != _abcissas.size()-1) {
				s << " ";
			};
		};
		s << "): " << get_ordinate();
		return s.str();
	};




































	/****************************************
	Grid pt interior
	****************************************/

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
};