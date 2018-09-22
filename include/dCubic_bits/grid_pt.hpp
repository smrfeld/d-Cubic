#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef UNORDERED_MAP_H
#define UNORDERED_MAP_H
#include <unordered_map>
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forwards
	class IdxSet;
	class GridPt;
	class GridPtOut;

	/****************************************
	Grid pt type
	****************************************/

	enum class GridPtType: unsigned int { INSIDE, OUTSIDE };

	/****************************************
	Interior grid pt
	****************************************/

	class GridPt {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		GridPt(IdxSet idxs, std::vector<double> abscissas);
		GridPt(const GridPt& other);
		GridPt(GridPt&& other);
		GridPt& operator=(const GridPt &other);
		GridPt& operator=(GridPt &&other);
		~GridPt();

		/********************
		Access
		********************/

		// Abscissa
		double get_abscissa(int dim) const;
		std::vector<double> get_abscissas() const;

		// Ordinate
		double get_ordinate() const;
		void set_ordinate(double val);
		void increment_ordinate(double val);

		// Idxs
		int get_idx(int dim) const;
		IdxSet get_idxs() const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;
	};
};