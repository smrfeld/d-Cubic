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
	bool IdxSet::find(int val) const {
		auto it = std::find(_idxs.begin(), _idxs.end(), val);
		if (it == _idxs.end()) {
			return false;
		} else {
			return true;
		};
	};

	// Get vector
	const std::vector<int>& IdxSet::get_vector_idxs() const {
		return _idxs;
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

    // Comparator
    bool operator==(const IdxSet &lhs, const IdxSet &rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		if (lhs[i] != rhs[i]) {
    			return false;
    		};
    	};
    	return true;
    };

	// Math
    IdxSet operator+(IdxSet lhs, const IdxSet& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs[i];
    	};
    	return lhs;
    };
    IdxSet operator-(IdxSet lhs, const IdxSet& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs[i];
    	};
    	return lhs;
    };
    IdxSet operator+(IdxSet lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs;
    	};
    	return lhs;
    };
    IdxSet operator-(IdxSet lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs;
    	};
    	return lhs;
    };

































	/****************************************
	IdxSet2
	****************************************/

    int IdxSet2::get_linear() const {
		int ret=0,add=0;
		for (auto i=0; i<_idxs.size(); i++) {
			add = _idxs[i];

			for (auto j=i+1; j<_idxs.size(); j++) {
				add *= 2;
			};
			ret += add;
		};

		return ret;
    };

	void IdxSet2::set_from_linear(int idx_linear) {
		int idx_working = idx_linear;

		// Determine the idxs
		int pwr;
		for (auto dim=0; dim<_idxs.size(); dim++) {
			pwr = 1;
			for (auto dim2=dim+1; dim2<_idxs.size(); dim2++) {
				pwr *= 2;
			};

			// Add
			_idxs[dim] = int(idx_working/pwr);

			// Remove
			idx_working -= _idxs[dim]*pwr;
		};
	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const IdxSet2& idxs) {
		stream << idxs.print();
		return stream;
    };

    // Comparator
    bool operator==(const IdxSet2 &lhs, const IdxSet2 &rhs) {
    	return lhs.get_linear() == rhs.get_linear();
    };
    bool operator<(const IdxSet2 &lhs, const IdxSet2 &rhs) {
    	return lhs.get_linear() < rhs.get_linear();
    };

    // Math
    IdxSet2 operator+(IdxSet2 lhs, const IdxSet2& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs[i];
    	};
    	return lhs;
    };
    IdxSet2 operator-(IdxSet2 lhs, const IdxSet2& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs[i];
    	};
    	return lhs;
    };
    IdxSet2 operator+(IdxSet2 lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs;
    	};
    	return lhs;
    };
    IdxSet2 operator-(IdxSet2 lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs;
    	};
    	return lhs;
    };


































	/****************************************
	IdxSet4
	****************************************/

    int IdxSet4::get_linear() const {
		int ret=0,add=0;
		for (auto i=0; i<_idxs.size(); i++) {
			add = _idxs[i];

			for (auto j=i+1; j<_idxs.size(); j++) {
				add *= 4;
			};
			ret += add;
		};

		return ret;
    };

	void IdxSet4::set_from_linear(int idx_linear) {
		int idx_working = idx_linear;

		// Determine the idxs
		int pwr;
		for (auto dim=0; dim<_idxs.size(); dim++) {
			pwr = 1;
			for (auto dim2=dim+1; dim2<_idxs.size(); dim2++) {
				pwr *= 4;
			};

			// Add
			_idxs[dim] = int(idx_working/pwr);

			// Remove
			idx_working -= _idxs[dim]*pwr;
		};
	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const IdxSet4& idxs) {
		stream << idxs.print();
		return stream;
    };

    // Comparator
    bool operator==(const IdxSet4 &lhs, const IdxSet4 &rhs) {
    	return lhs.get_linear() == rhs.get_linear();
    };
    bool operator<(const IdxSet4 &lhs, const IdxSet4 &rhs) {
    	return lhs.get_linear() < rhs.get_linear();
    };
    
    // Math
    IdxSet4 operator+(IdxSet4 lhs, const IdxSet4& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs[i];
    	};
    	return lhs;
    };
    IdxSet4 operator-(IdxSet4 lhs, const IdxSet4& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs[i];
    	};
    	return lhs;
    };
    IdxSet4 operator+(IdxSet4 lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs;
    	};
    	return lhs;
    };
    IdxSet4 operator-(IdxSet4 lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs;
    	};
    	return lhs;
    };


};