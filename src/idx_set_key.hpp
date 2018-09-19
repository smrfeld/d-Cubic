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
	enum class GridPtType: unsigned int { INSIDE, OUTSIDE };

	/****************************************
	Index set key
	****************************************/

	class IdxSetKey {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		IdxSetKey(IdxSet idxs, int no_idxs_possible);
		IdxSetKey(IdxSet idxs, std::vector<std::shared_ptr<Dimension1D>> dims);
		IdxSetKey(const IdxSetKey& other);
		IdxSetKey(IdxSetKey&& other);
		IdxSetKey& operator=(const IdxSetKey &other);
		IdxSetKey& operator=(IdxSetKey &&other);
		~IdxSetKey();

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
    std::ostream& operator<< (std::ostream& stream, const IdxSetKey& idxs);

	// Comparator
	bool operator <(const IdxSetKey& x, const IdxSetKey& y);
	bool operator ==(const IdxSetKey& x, const IdxSetKey& y);
};