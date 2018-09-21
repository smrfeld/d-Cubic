#include "../include/dcubic_bits/dimension_1d.hpp"

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Dimension1D
	****************************************/

	Dimension1D::Dimension1D(double start_pt, double end_pt, int no_pts) {
		_start_pt = start_pt;
		_end_pt = end_pt;
		_no_pts = no_pts;
		_spacing = (_end_pt - _start_pt) / (_no_pts - 1.0);
		for (int i=0; i<_no_pts; i++) {
			_pts.push_back(_start_pt + i*_spacing);
		};
	};

	// Accessors
	int Dimension1D::get_no_pts() const {
		return _no_pts;
	};
	double Dimension1D::get_start_pt() const {
		return _start_pt;
	};
	double Dimension1D::get_end_pt() const {
		return _end_pt;
	};
	double Dimension1D::get_spacing() const {
		return _spacing;
	};

	double Dimension1D::get_pt_at_idx(int idx) const {
		if (idx >= 0 && idx <= _no_pts-1) {
			return _pts[idx];
		} else {
			return _start_pt + idx*_spacing;
		};
	};
	std::vector<double> Dimension1D::get_pts() const {
		return _pts;
	};

	std::pair<bool,std::pair<int,int>> Dimension1D::get_surrounding_idxs(double pt) const {
		if (pt < _start_pt || pt > _end_pt) {
			// Outside
			return std::make_pair(false,std::make_pair(0,0));
		};

		int idx;
		if (pt == _end_pt) {
			// At end
			idx = _no_pts-2;
		} else {
			idx = int((pt - _start_pt) / _spacing);
		};

		return std::make_pair(true,std::make_pair(idx,idx+1));
	};
};