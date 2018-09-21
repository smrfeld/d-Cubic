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

		// Helpers
		void _clean_up();
		void _copy(const IdxSet& other);
		void _move(IdxSet& other);

	protected:
		
		// Idxs
		std::vector<int> _idxs;

	public:

		/********************
		Constructor
		********************/

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

		const std::vector<int>& get_vector_idxs() const;

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




























	/****************************************
	Index set 2
	****************************************/

	class IdxSet2 : public IdxSet {

	public:

		/********************
		Constructor
		********************/

		IdxSet2(int no_idxs);
		IdxSet2(std::vector<int> idxs);

		/********************
		Accessors
		********************/

		int get_linear() const;
		void set_from_linear(int idx_linear);

	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const IdxSet2& idxs);

    // Comparator
    bool operator==(const IdxSet2 &lhs, const IdxSet2 &rhs);
    bool operator<(const IdxSet2 &lhs, const IdxSet2 &rhs);

    // Math
    IdxSet2 operator+(IdxSet2 lhs, const IdxSet2& rhs);
    IdxSet2 operator-(IdxSet2 lhs, const IdxSet2& rhs);
    IdxSet2 operator+(IdxSet2 lhs, int rhs);
    IdxSet2 operator-(IdxSet2 lhs, int rhs);























	/****************************************
	Index set 4
	****************************************/

	class IdxSet4 : public IdxSet {

	public:

		/********************
		Constructor
		********************/

		IdxSet4(int no_idxs);
		IdxSet4(std::vector<int> idxs);

		/********************
		Accessors
		********************/

		int get_linear() const;
		void set_from_linear(int idx_linear);

	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const IdxSet4& idxs);

    // Comparator
    bool operator==(const IdxSet4 &lhs, const IdxSet4 &rhs);
    bool operator<(const IdxSet4 &lhs, const IdxSet4 &rhs);

    // Math
    IdxSet4 operator+(IdxSet4 lhs, const IdxSet4& rhs);
    IdxSet4 operator-(IdxSet4 lhs, const IdxSet4& rhs);
    IdxSet4 operator+(IdxSet4 lhs, int rhs);
    IdxSet4 operator-(IdxSet4 lhs, int rhs);





};

