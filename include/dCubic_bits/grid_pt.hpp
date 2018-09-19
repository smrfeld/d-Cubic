#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

#ifndef MAP_H
#define MAP_H
#include <map>
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forwards
	class IdxSet;
	class GridPtKey;
	class GridPt;
	class GridPtOut;

	/****************************************
	Grid pt type
	****************************************/

	enum class GridPtType: unsigned int { INSIDE, OUTSIDE };

	/****************************************
	Neighborhood of points surrounding a point, 2 in each dim
	****************************************/

	struct Nbr2 {
		std::map<GridPtKey, const GridPt*> in;
	};

	/****************************************
	Neighborhood of points surrounding a point, 4 in each dim
	****************************************/

	struct Nbr4 {
		std::map<GridPtKey, GridPtType> types;
		std::map<GridPtKey, const GridPt*> in;
		std::map<GridPtKey, const GridPtOut*> out;
	};

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

		// Idxs
		int get_idx(int dim) const;
		IdxSet get_idxs() const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;
	};
};