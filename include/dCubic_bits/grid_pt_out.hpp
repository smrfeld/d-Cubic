#include <vector>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	// Forward
	class GridPt;
	class IdxSet;

	/****************************************
	Location of pt in each dim
	****************************************/

	enum class LocInDim: unsigned int {P0_OUTSIDE, INSIDE, P3_OUTSIDE};

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
		GridPtOut(IdxSet idxs, std::vector<double> abscissas, const GridPt* p1, const GridPt* p2, std::vector<LocInDim> loc);
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
		const GridPt* get_dep_p1() const;
		const GridPt* get_dep_p2() const;

		// Which dims are outside
		const std::vector<LocInDim>& get_locs() const;
		LocInDim get_loc_in_dim(int dim) const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;
		
	};

};