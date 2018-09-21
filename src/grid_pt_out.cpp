#include "../include/dCubic_bits/grid_pt_out.hpp"

// Other headers
#include "../include/dCubic_bits/grid_pt.hpp"
#include "../include/dCubic_bits/idx_set.hpp"

#include <iostream>
#include <sstream>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Implementation header
	****************************************/

	class GridPtOut::Impl {

	private:

		// Idxs
		IdxSet _idxs;

		// Abscissa values
		std::vector<double> _abcissas;

		// Dependencies
		const GridPt* _p1, *_p2;

		// Which dims are outside
		std::vector<LocInDim> _locs;
		
		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _move(Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(IdxSet idxs, std::vector<double> abscissas, const GridPt* p1, const GridPt* p2, std::vector<LocInDim> loc);
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

		// Idxs
		int get_idx(int dim) const;
		IdxSet get_idxs() const;

		// Get dependent points
		const GridPt* get_dep_p1() const;
		const GridPt* get_dep_p2() const;

		// Which dims are outside
		const std::vector<LocInDim>& get_locs() const;
		LocInDim get_loc_in_dim(int dim) const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;
	};






















	/****************************************
	Implementation
	****************************************/

	GridPtOut::Impl::Impl(IdxSet idxs, std::vector<double> abscissas, const GridPt* p1, const GridPt* p2, std::vector<LocInDim> loc) : _p1(p1), _p2(p2), _idxs(idxs) {
		// Check lengths
		if (idxs.size() != abscissas.size()) {
			std::cout << ">>> Error: GridPt::Impl::Impl <<< Sizes must match." << std::endl;
			exit(EXIT_FAILURE);
		};

		// Store
		_abcissas = abscissas;
		_locs = loc;
	};
	GridPtOut::Impl::Impl(const Impl& other) : _idxs(other._idxs) {
		_copy(other);
	};
	GridPtOut::Impl::Impl(Impl&& other) : _idxs(std::move(other._idxs)) {
		_move(other);
	};
    GridPtOut::Impl& GridPtOut::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    GridPtOut::Impl& GridPtOut::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	GridPtOut::Impl::~Impl()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void GridPtOut::Impl::_clean_up()
	{
		// Nothing...
	};
	void GridPtOut::Impl::_copy(const Impl& other)
	{
		_idxs = other._idxs;
		_abcissas = other._abcissas;
		_p1 = other._p1;
		_p2 = other._p2;
		_locs = other._locs;
	};
	void GridPtOut::Impl::_move(Impl& other)
	{
		_copy(other);

		// Reset other
		// other._idxs = IdxSet();
		other._abcissas.clear();
		other._p1 = nullptr;
		other._p2 = nullptr;
		other._locs.clear();
	};

	/********************
	Access
	********************/

	// Abscissa
	double GridPtOut::Impl::get_abscissa(int dim) const {
		return _abcissas[dim];
	};
	std::vector<double> GridPtOut::Impl::get_abscissas() const {
		return _abcissas;
	};

	// Ordinate
	double GridPtOut::Impl::get_ordinate() const {
		return 2.0 * _p1->get_ordinate() - _p2->get_ordinate(); 
	};

	// Idxs
	int GridPtOut::Impl::get_idx(int dim) const {
		return _idxs[dim];
	};
	IdxSet GridPtOut::Impl::get_idxs() const {
		return _idxs;
	};

	// Get dependent points
	const GridPt* GridPtOut::Impl::get_dep_p1() const {
		return _p1;
	};
	const GridPt* GridPtOut::Impl::get_dep_p2() const {
		return _p2;
	};

	// Which dims are outside
	const std::vector<LocInDim>& GridPtOut::Impl::get_locs() const {
		return _locs;
	};
	LocInDim GridPtOut::Impl::get_loc_in_dim(int dim) const {
		if (dim >= _locs.size() || dim < 0) {
			std::cout << ">>> Error: GridPtOut::Impl::get_loc_in_dim <<< dim: " << dim << " is outside 0-" << _locs.size() << std::endl;
			exit(EXIT_FAILURE);
		};
		return _locs[dim];
	};

	/********************
	Print
	********************/

	std::string GridPtOut::Impl::print_abscissa() const {
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

	GridPtOut::GridPtOut(IdxSet idxs, std::vector<double> abscissas, const GridPt* p1, const GridPt* p2,  std::vector<LocInDim> loc) : _impl(new Impl(idxs,abscissas, p1, p2, loc)) {};
	GridPtOut::GridPtOut(const GridPtOut& other) : _impl(new Impl(*other._impl)) {};
	GridPtOut::GridPtOut(GridPtOut&& other) : _impl(std::move(other._impl)) {};
	GridPtOut& GridPtOut::operator=(const GridPtOut &other) {
        _impl.reset(new Impl(*other._impl));
        return *this; 
	};
	GridPtOut& GridPtOut::operator=(GridPtOut &&other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	GridPtOut::~GridPtOut() = default;

	/********************
	Access
	********************/

	// Abscissa
	double GridPtOut::get_abscissa(int dim) const {
		return _impl->get_abscissa(dim);
	};
	std::vector<double> GridPtOut::get_abscissas() const {
		return _impl->get_abscissas();
	};

	// Ordinate
	double GridPtOut::get_ordinate() const {
		return _impl->get_ordinate();
	};

	// Idxs
	int GridPtOut::get_idx(int dim) const {
		return _impl->get_idx(dim);
	};
	IdxSet GridPtOut::get_idxs() const {
		return _impl->get_idxs();
	};

	// Get dependent points
	const GridPt* GridPtOut::get_dep_p1() const {
		return _impl->get_dep_p1();
	};
	const GridPt* GridPtOut::get_dep_p2() const {
		return _impl->get_dep_p2();
	};

	// Which dims are outside
	const std::vector<LocInDim>& GridPtOut::get_locs() const {
		return _impl->get_locs();
	};
	LocInDim GridPtOut::get_loc_in_dim(int dim) const {
		return _impl->get_loc_in_dim(dim);
	};

	/********************
	Print
	********************/

	std::string GridPtOut::print_abscissa() const {
		return _impl->print_abscissa();
	};

};