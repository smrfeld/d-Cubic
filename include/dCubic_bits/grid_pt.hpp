#include <vector>

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
	Interior grid pt
	****************************************/

	class GridPt {

	private:

		// Abscissa values
		std::vector<double> _abcissas;

		// Ordinate
		double _ordinate;
		
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
		const double& get_ordinate_const_ref() const;
		void set_ordinate(double val);
		void increment_ordinate(double val);

		// Update
		void set_update(double val);
		void increment_update(double val);
		void multiply_update(double val);
		void reset_update();
		double get_update() const;
		void committ_update(); // automatically resets

		// Idxs
		int get_idx(int dim) const;
		IdxSet get_idxs() const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;
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
		const double& get_ordinate_const_ref() const;
		void set_ordinate(double val);
		void increment_ordinate(double val);

		// Update
		void set_update(double val);
		void increment_update(double val);
		void multiply_update(double val);
		void reset_update();
		double get_update() const;
		void committ_update(); // automatically resets

		// Idxs
		int get_idx(int dim) const;
		IdxSet get_idxs() const;

		/********************
		Print
		********************/

		std::string print_abscissa() const;
	};
};