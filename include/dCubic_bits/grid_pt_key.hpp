#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forward
	class IdxSet;
	class Dimension1D;
	enum class GridPtType: unsigned int;

	/****************************************
	Grid pt key
	****************************************/

	class GridPtKey {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		GridPtKey(IdxSet idxs, int no_idxs_possible_in_each_dim);
		GridPtKey(IdxSet idxs, const std::vector<Dimension1D*> dims);
		GridPtKey(IdxSet idxs, const std::vector<Dimension1D>& dims);
		GridPtKey(const GridPtKey& other);
		GridPtKey(GridPtKey&& other);
		GridPtKey& operator=(const GridPtKey &other);
		GridPtKey& operator=(GridPtKey &&other);
		~GridPtKey();

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

	/********************
	Comparators
	********************/

	// Printing
    std::ostream& operator<< (std::ostream& stream, const GridPtKey& idxs);

	// Comparator
	bool operator <(const GridPtKey& x, const GridPtKey& y);
	bool operator ==(const GridPtKey& x, const GridPtKey& y);
};