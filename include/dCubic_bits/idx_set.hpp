#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Index set
	****************************************/

	class IdxSet {

	private:
		
		// Idxs
		std::vector<int> _idxs;

	public:

		/********************
		Constructor
		********************/

		IdxSet();
		IdxSet(int no_idxs);
		IdxSet(std::vector<int> idxs);

		/********************
		Accessors
		********************/

		int operator [](int idx) const;
		int & operator [](int idx);

		int size() const;

		bool find(int val);

		std::string print() const;
	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const IdxSet& idxs);
};