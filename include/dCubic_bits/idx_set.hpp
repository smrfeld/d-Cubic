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

		// Helpers
		void _clean_up();
		void _copy(const IdxSet& other);
		void _move(IdxSet& other);

	public:

		/********************
		Constructor
		********************/

		IdxSet();
		IdxSet(int no_idxs);
		IdxSet(std::vector<int> idxs);
		IdxSet(const IdxSet& other);
		IdxSet(IdxSet&& other);
		IdxSet& operator=(const IdxSet &other);
		IdxSet& operator=(IdxSet &&other);
		~IdxSet();

		/********************
		Accessors
		********************/

		int operator [](int idx) const;
		int & operator [](int idx);

		int size() const;

		bool find(int val) const;

		std::vector<int> get_vector_idxs() const;

		std::string print() const;
	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const IdxSet& idxs);

    // Comparator
    bool operator==(const IdxSet &lhs, const IdxSet &rhs);

    // Math
    IdxSet operator+(IdxSet lhs, const IdxSet& rhs);
    IdxSet operator-(IdxSet lhs, const IdxSet& rhs);
    IdxSet operator+(IdxSet lhs, int rhs);
    IdxSet operator-(IdxSet lhs, int rhs);

};

