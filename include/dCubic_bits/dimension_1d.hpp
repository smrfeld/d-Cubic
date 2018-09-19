#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Dimension1D
	****************************************/

	class Dimension1D {

	private:

		double _start_pt;
		double _end_pt;
		double _spacing;
		int _no_pts;
		std::vector<double> _pts;

	public:

		Dimension1D(double start_pt, double end_pt, int no_pts);

		// Accessors
		int get_no_pts() const;
		double get_start_pt() const;
		double get_end_pt() const;
		double get_spacing() const;

		double get_pt_at_idx(int idx) const;
		std::vector<double> get_pts() const;

		std::pair<bool,std::pair<int,int>> get_surrounding_idxs(double pt) const;
	};
};