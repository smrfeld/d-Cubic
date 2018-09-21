#include "../include/dcubic_bits/grid_pt.hpp"

// Other headers
#include "../include/dcubic_bits/idx_set.hpp"
#include "../include/dcubic_bits/grid_pt_out.hpp"

#include <iostream>
#include <sstream>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Implementation header
	****************************************/

	class GridPt::Impl {

	private:

		// Abscissa values
		std::vector<double> _abcissas;

		// Ordinate
		double _ordinate;

		// Indexes
		IdxSet _idxs;

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _move(Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(IdxSet idxs, std::vector<double> abscissas);
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl &other);
		Impl& operator=(Impl &&other);
		~Impl();

		/********************
		Access
		********************/

		// Abscissa
		double get_abscissa(int dim) const;
		std::vector<double> get_abscissas() const;

		// Ordinate
		double get_ordinate() const;
		void set_ordinate(double val);

		// Idxs
		int get_idx(int dim) const;
		IdxSet get_idxs() const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;

	};






















	/****************************************
	Implementation
	****************************************/

	GridPt::Impl::Impl(IdxSet idxs, std::vector<double> abscissas) : _idxs(idxs) {
		// Check lengths
		if (idxs.size() != abscissas.size()) {
			std::cout << ">>> Error: GridPt::Impl::Impl <<< Sizes must match." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Store
		_abcissas = abscissas;

		// Init
		_ordinate = 0.0;
	};
	GridPt::Impl::Impl(const Impl& other) : _idxs(other._idxs) {
		_copy(other);
	};
	GridPt::Impl::Impl(Impl&& other) : _idxs(std::move(other._idxs)) {
		_move(other);
	};
    GridPt::Impl& GridPt::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    GridPt::Impl& GridPt::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	GridPt::Impl::~Impl()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void GridPt::Impl::_clean_up()
	{
		// Nothing...
	};
	void GridPt::Impl::_copy(const Impl& other)
	{
		_abcissas = other._abcissas;
		_ordinate = other._ordinate;
		_idxs = other._idxs;
	};
	void GridPt::Impl::_move(Impl& other)
	{
		_abcissas = other._abcissas;
		_ordinate = other._ordinate;
		_idxs = other._idxs;

		// Reset other
		other._abcissas.clear();
		other._ordinate = 0.0;
		// other._idxs = IdxSet();
	};

	/********************
	Access
	********************/

	// Abscissa
	double GridPt::Impl::get_abscissa(int dim) const {
		return _abcissas[dim];
	};
	std::vector<double> GridPt::Impl::get_abscissas() const {
		return _abcissas;
	};

	// Ordinate
	double GridPt::Impl::get_ordinate() const {
		return _ordinate;
	};
	void GridPt::Impl::set_ordinate(double val) {
		_ordinate = val;
	};

	// Idxs
	int GridPt::Impl::get_idx(int dim) const {
		return _idxs[dim];
	};
	IdxSet GridPt::Impl::get_idxs() const {
		return _idxs;
	};

	/********************
	Print
	********************/

	std::string GridPt::Impl::print_abscissa() const {
		std::ostringstream s;
		s << "(";
		for (auto dim=0; dim<_abcissas.size(); dim++) {
			s << _abcissas[dim];
			if (dim != _abcissas.size()-1) {
				s << " ";
			};
		};
		s << ")";
		return s.str();
	};





















	/****************************************
	Impl forwards
	****************************************/

	/********************
	Constructor
	********************/

	GridPt::GridPt(IdxSet idxs, std::vector<double> abscissas) : _impl(new Impl(idxs, abscissas)) {};
	GridPt::GridPt(const GridPt& other) : _impl(new Impl(*other._impl)) {};
	GridPt::GridPt(GridPt&& other) : _impl(std::move(other._impl)) {};
	GridPt& GridPt::operator=(const GridPt &other) {
        _impl.reset(new Impl(*other._impl));
        return *this; 
	};
	GridPt& GridPt::operator=(GridPt &&other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	GridPt::~GridPt() = default;

	/********************
	Access
	********************/

	// Abscissa
	double GridPt::get_abscissa(int dim) const {
		return _impl->get_abscissa(dim);
	};
	std::vector<double> GridPt::get_abscissas() const {
		return _impl->get_abscissas();
	};

	// Ordinate
	double GridPt::get_ordinate() const {
		return _impl->get_ordinate();
	};
	void GridPt::set_ordinate(double val) {
		_impl->set_ordinate(val);
	};

	// Idxs
	int GridPt::get_idx(int dim) const {
		return _impl->get_idx(dim);
	};
	IdxSet GridPt::get_idxs() const {
		return _impl->get_idxs();
	};

	/********************
	Print
	********************/

	std::string GridPt::print_abscissa() const {
		return _impl->print_abscissa();
	};

};