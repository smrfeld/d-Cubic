#include "../include/dCubic_bits/idx_set.hpp"

#include <iostream>
#include <sstream>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Index set
	****************************************/

	// Constructors
	IdxSet::IdxSet() {};
	IdxSet::IdxSet(int no_idxs) {
		for (auto i=0; i<no_idxs; i++) {
			_idxs.push_back(0);
		};
	};
	IdxSet::IdxSet(std::vector<int> idxs) {
		_idxs = idxs;
	};
	IdxSet::IdxSet(const IdxSet& other) {
		_copy(other);
	};
	IdxSet::IdxSet(IdxSet&& other) {
		_move(other);
	};
    IdxSet& IdxSet::operator=(const IdxSet& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    IdxSet& IdxSet::operator=(IdxSet&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	IdxSet::~IdxSet()
	{
		_clean_up();
	};

	// Helpers
	void IdxSet::_clean_up()
	{
		// Nothing...
	};
	void IdxSet::_copy(const IdxSet& other)
	{
		_idxs = other._idxs;
	};
	void IdxSet::_move(IdxSet& other)
	{
		_idxs = std::move(other._idxs);
	};

	// Accessors
	int IdxSet::operator [](int idx) const {
		return _idxs[idx];
	};
	int & IdxSet::operator [](int idx) {
		return _idxs[idx];
	};

	// Size
	int IdxSet::size() const {
		return _idxs.size();
	};

	// Find
	bool IdxSet::find(int val) {
		auto it = std::find(_idxs.begin(), _idxs.end(), val);
		if (it == _idxs.end()) {
			return false;
		} else {
			return true;
		};
	};

 	// Print
 	std::string IdxSet::print() const {
 		std::ostringstream stream;
		stream << "(";
		for (auto i=0; i<_idxs.size(); i++) {
			stream << _idxs[i];
			if (i != _idxs.size()-1) {
				stream << " ";
			};
		};
		stream << ")";
	    return stream.str();
 	};

	// Printing
	std::ostream& operator<<(std::ostream& stream, const IdxSet& idxs) {
		stream << idxs.print();
		return stream;
	 };
};