#ifndef STRING_H
#define STRING_H
#include <string>
#endif

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
	class Dimension1D;
	class GridPt;
	class GridPtOut;
	class GridPtKey;
	class IdxSet;
	struct Nbr4;
	struct Nbr2;

	/****************************************
	Grid
	****************************************/

	class Grid {

	private:

		class Impl;
		std::unique_ptr<Impl> _impl;

	public:

		/********************
		Constructor
		********************/

		Grid(std::vector<std::shared_ptr<Dimension1D>> dims);
		Grid(const Grid& other);
		Grid(Grid&& other);
		Grid& operator=(const Grid &other);
		Grid& operator=(Grid &&other);
		~Grid();

		/********************
		Get dims
		********************/

		int get_no_dims() const;
		std::vector<std::shared_ptr<Dimension1D>> get_dims() const;

		/********************
		Get grid points
		********************/

		std::map<GridPtKey, std::shared_ptr<GridPt>> get_grid_points() const;
		std::shared_ptr<GridPt> get_grid_point(std::vector<int> grid_idxs) const;
		std::shared_ptr<GridPt> get_grid_point(IdxSet grid_idxs) const;
		std::shared_ptr<GridPt> get_grid_point(GridPtKey key) const;

		std::map<GridPtKey, std::shared_ptr<GridPtOut>> get_grid_points_outside() const;
		std::shared_ptr<GridPtOut> get_grid_point_outside(std::vector<int> grid_idxs) const;
		std::shared_ptr<GridPtOut> get_grid_point_outside(IdxSet grid_idxs) const;
		std::shared_ptr<GridPtOut> get_grid_point_outside(GridPtKey key) const;

		/********************
		Get grid points surrounding a point
		********************/

		Nbr2 get_surrounding_2(std::vector<double> abscissas) const;
		Nbr4 get_surrounding_4(std::vector<double> abscissas) const;

		/********************
		Project
		********************/

		void project();

		/********************
		Write
		********************/

		void write_solution(std::string fname) const;

	};

};