#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forward
	class GridPt;
	class IdxSet;

	/****************************************
	Exterior grid pt
	****************************************/

	class GridPtOut {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		// Value = 2 p1 - p2
		GridPtOut(IdxSet idxs, std::vector<double> abscissas, std::shared_ptr<GridPt> p1, std::shared_ptr<GridPt> p2);
		GridPtOut(const GridPtOut& other);
		GridPtOut(GridPtOut&& other);
		GridPtOut& operator=(const GridPtOut &other);
		GridPtOut& operator=(GridPtOut &&other);
		~GridPtOut();

		/********************
		Accessors
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
		std::shared_ptr<GridPt> get_dep_p1() const;
		std::shared_ptr<GridPt> get_dep_p2() const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;
		
	};

};