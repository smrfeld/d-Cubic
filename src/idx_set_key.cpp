#include "idx_set_key.hpp"
#include "../include/dCubic_bits/dimension_1d.hpp"
#include "../include/dCubic_bits/idx_set.hpp"

#include <iostream>
#include <sstream>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Declaration
	****************************************/

	class IdxSetKey::Impl {

	private:

		// Inside or outside
		GridPtType _type;

		// Idxs
		IdxSet _idxs;
		std::vector<int> _no_idxs_possible;

		// Constructor helpers
		void _clean_up();
		void _copy(const Impl& other);
		void _move(Impl& other);

	public:

		/********************
		Constructor
		********************/

		Impl(IdxSet idxs, int no_idxs_possible);
		Impl(IdxSet idxs, std::vector<std::shared_ptr<Dimension1D>> dims);
		Impl(const Impl& other);
		Impl(Impl&& other);
		Impl& operator=(const Impl &other);
		Impl& operator=(Impl &&other);
		~Impl();

		/********************
		Accessors
		********************/

		// Type
		GridPtType get_type() const;

		// Linear idx
		int get_linear() const;

		// Set from linear
		void set_from_linear(int idx_linear);

		std::string print() const;
	};






























	/****************************************
	IdxSetKey
	****************************************/

	/********************
	Constructor
	********************/

	IdxSetKey::Impl::Impl(IdxSet idxs, int no_idxs_possible) {
		for (auto i=0; i<idxs.size(); i++) {
			_no_idxs_possible.push_back(no_idxs_possible);
		};

		// Check if outside
		_type = GridPtType::INSIDE;
		for (auto dim=0; dim<idxs.size(); dim++) {
			if (idxs[dim] < 0 || idxs[dim] > _no_idxs_possible[dim]-1) {
				_type = GridPtType::OUTSIDE;
				break;
			};
		};

		// Set vals
		_idxs = idxs;
	};
	IdxSetKey::Impl::Impl(IdxSet idxs, std::vector<std::shared_ptr<Dimension1D>> dims) {
		for (auto i=0; i<dims.size(); i++) {
			// Add 2 to take into account outside dims
			_no_idxs_possible.push_back(dims[i]->get_no_pts()+2);
		};

		// Check if grid pt is outside
		_type = GridPtType::INSIDE;
		for (auto dim=0; dim<idxs.size(); dim++) {
			if (idxs[dim] < 0 || idxs[dim] > _no_idxs_possible[dim]-1) {
				_type = GridPtType::OUTSIDE;
				break;
			};
		};

		// Set vals
		_idxs = idxs;
	};
	IdxSetKey::Impl::Impl(const Impl& other) {
		_copy(other);
	};
	IdxSetKey::Impl::Impl(Impl&& other) {
		_move(other);
	};
    IdxSetKey::Impl& IdxSetKey::Impl::operator=(const Impl& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    IdxSetKey::Impl& IdxSetKey::Impl::operator=(Impl&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	IdxSetKey::Impl::~Impl()
	{
		_clean_up();
	};

	/********************
	Helpers for constructors
	********************/

	void IdxSetKey::Impl::_clean_up()
	{
		// Nothing...
	};
	void IdxSetKey::Impl::_copy(const Impl& other)
	{
		_type = other._type;
		_idxs = other._idxs;
		_no_idxs_possible = other._no_idxs_possible;
	};
	void IdxSetKey::Impl::_move(Impl& other)
	{
		_type = other._type;
		_idxs = other._idxs;
		_no_idxs_possible = other._no_idxs_possible;

		// Reset other
		other._idxs = IdxSet();
		other._no_idxs_possible.clear();
	};

	/********************
	Accessors
	********************/

	// Type
	GridPtType IdxSetKey::Impl::get_type() const {
		return _type;
	};

	// Linear idx
	int IdxSetKey::Impl::get_linear() const {
		int ret=0,add=0;
		for (auto i=0; i<_idxs.size(); i++) {
			add = _idxs[i];

			// Add 1 for outside because starts at -1 instead of 0
			add += 1;

			for (auto j=i+1; j<_idxs.size(); j++) {
				add *= _no_idxs_possible[j];
			};
			ret += add;
		};

		return ret;
	};

	// Set from linear
	void IdxSetKey::Impl::set_from_linear(int idx_linear) {
		int idx_working = idx_linear;

		// Determine the idxs
		int pwr;
		for (auto dim=0; dim<_no_idxs_possible.size(); dim++) {
			pwr = 1;
			for (auto dim2=dim+1; dim2<_no_idxs_possible.size(); dim2++) {
				pwr *= _no_idxs_possible[dim2];
			};

			// Add
			_idxs[dim] = int(idx_working/pwr);

			// Remove
			idx_working -= _idxs[dim]*pwr;

			// Remove 1 for outside because starts at -1 instead of 0
			_idxs[dim] -= 1;
		};
	};
 
 	// Print
 	std::string IdxSetKey::Impl::print() const {
 		return _idxs.print();
 	};





































	/****************************************
	Forwards
	****************************************/

	/********************
	Constructor
	********************/

	IdxSetKey::IdxSetKey(IdxSet idxs, int no_idxs_possible) : _impl(new Impl(idxs,no_idxs_possible)) {};
	IdxSetKey::IdxSetKey(IdxSet idxs, std::vector<std::shared_ptr<Dimension1D>> dims) : _impl(new Impl(idxs,dims)) {};
	IdxSetKey::IdxSetKey(const IdxSetKey& other) : _impl(new Impl(*other._impl)) {};
	IdxSetKey::IdxSetKey(IdxSetKey&& other) : _impl(std::move(other._impl)) {};
	IdxSetKey& IdxSetKey::operator=(const IdxSetKey &other)  {
        _impl.reset(new Impl(*other._impl));
        return *this; 
	};	
	IdxSetKey& IdxSetKey::operator=(IdxSetKey &&other) {
        _impl = std::move(other._impl);
        return *this; 
	};
	IdxSetKey::~IdxSetKey() = default;

	/********************
	Accessors
	********************/

	// Type
	GridPtType IdxSetKey::get_type() const {
		return _impl->get_type();
	};

	// Linear idx
	int IdxSetKey::get_linear() const {
		return _impl->get_linear();
	};

	// Set from linear
	void IdxSetKey::set_from_linear(int idx_linear) {
		_impl->set_from_linear(idx_linear);
	};

	std::string IdxSetKey::print() const {
		return _impl->print();
	};

	/********************
	Comparators
	********************/

	// Printing
	std::ostream& operator<<(std::ostream& stream, const IdxSetKey& idxs) {
		stream << idxs.print();
		return stream;
	 };

	// Comparator
	bool operator <(const IdxSetKey& x, const IdxSetKey& y) {
		return x.get_linear() < y.get_linear();
	};
	bool operator ==(const IdxSetKey& x, const IdxSetKey& y) {
		return x.get_linear() == y.get_linear();
	};
};