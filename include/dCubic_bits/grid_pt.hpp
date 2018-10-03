#include <vector>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Grid pt
	--- CAUTION: abstract base!!! ---
	****************************************/

	class GridPt {

	private:

		// Abscissa values
		std::vector<double> _abcissas;

		// Constructor helpers
		void _clean_up();
		void _copy(const GridPt& other);
		void _move(GridPt& other);

	public:

		/********************
		Constructor
		********************/

		GridPt(std::vector<double> abscissas);
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
		// --- PURE ---
		virtual double get_ordinate() const = 0;

		/********************
		Print
		********************/

		std::string print() const;
	};



































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

		GridPtIn(std::vector<double> abscissas);
		GridPtIn(const GridPtIn& other);
		GridPtIn(GridPtIn&& other);
		GridPtIn& operator=(const GridPtIn &other);
		GridPtIn& operator=(GridPtIn &&other);
		~GridPtIn();

		/********************
		Access
		********************/

		// Ordinate
		double get_ordinate() const;
		const double& get_ordinate_const_ref() const;
		void set_ordinate(double val);
		void increment_ordinate(double val);

	};
};