#include <vector>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Type
	****************************************/

	enum class GridPtType: unsigned int {INSIDE, OUTSIDE};

	/****************************************
	Grid pt
	--- CAUTION: abstract base!!! ---
	****************************************/

	class GridPt {

	protected:

		// No dims
		int _no_dims;

	private:

		// Abscissa values
		double* _abcissas;

		// Constructor helpers
		void _clean_up();
		void _copy(const GridPt& other);
		void _move(GridPt& other);

	public:

		/********************
		Constructor
		********************/

		GridPt(int no_dims, double* abscissas);
		GridPt(std::vector<double> abscissas);
		GridPt(const GridPt& other);
		GridPt(GridPt&& other);
		GridPt& operator=(const GridPt &other);
		GridPt& operator=(GridPt &&other);
		~GridPt();

		/********************
		Access
		********************/

		// Type
		virtual GridPtType get_type() const = 0;

		// Abscissa
		double get_abscissa(int dim) const;

		// Ordinate
		// --- PURE ---
		virtual double get_ordinate() const = 0;

		/********************
		Print
		********************/

		std::string print() const;
	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const GridPt& pt);


































	/****************************************
	Interior grid pt
	****************************************/

	class GridPtIn : public GridPt {

	private:

		// Ordinate
		double _ordinate;

		// Constructor helpers
		void _clean_up();
		void _copy(const GridPtIn& other);
		void _move(GridPtIn& other);

	public:

		/********************
		Constructor
		********************/

		GridPtIn(int no_dims, double* abscissas);
		GridPtIn(std::vector<double> abscissas);
		GridPtIn(const GridPtIn& other);
		GridPtIn(GridPtIn&& other);
		GridPtIn& operator=(const GridPtIn &other);
		GridPtIn& operator=(GridPtIn &&other);
		~GridPtIn();

		/********************
		Access
		********************/

		// Type
		GridPtType get_type() const;

		// Ordinate
		double get_ordinate() const;
		const double& get_ordinate_const_ref() const;
		void set_ordinate(double val);
		void increment_ordinate(double val);

	};





































	/****************************************
	Exterior grid pt
	****************************************/

	// Location of pt in each dim
	enum class Loc: unsigned int {OUTSIDE_LEFT, INSIDE, OUTSIDE_RIGHT};

	class GridPtOut : public GridPt {

	private:

		// Dependent pts
		const GridPtIn* _p1;
		const GridPtIn* _p2;

		// Locs
		Loc* _locs;

		// Constructor helpers
		void _clean_up();
		void _copy(const GridPtOut& other);
		void _move(GridPtOut& other);

	public:

		/********************
		Constructor
		********************/

		GridPtOut(int no_dims, double* abscissas, const GridPtIn* p1, const GridPtIn* p2, Loc* locs);
		GridPtOut(std::vector<double> abscissas, const GridPtIn* p1, const GridPtIn* p2, std::vector<Loc> locs);
		GridPtOut(const GridPtOut& other);
		GridPtOut(GridPtOut&& other);
		GridPtOut& operator=(const GridPtOut &other);
		GridPtOut& operator=(GridPtOut &&other);
		~GridPtOut();

		/********************
		Access
		********************/

		// Type
		GridPtType get_type() const;

		// Ordinate
		double get_ordinate() const;

		// Get dependent points
		const GridPtIn* get_dep_p1() const;
		const GridPtIn* get_dep_p2() const;

		// Which dims are outside
		Loc get_loc_in_dim(int dim) const;
	};
};