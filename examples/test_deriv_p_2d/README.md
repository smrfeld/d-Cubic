# Take derivative with respect to p in 2D

See `test_deriv_p_2d.cpp` and for validation the Mathematica file `test_deriv_p_2d.nb`.

Make two dimensions:
```
Dimension1D dim_1(0.0,1.0,30);
Dimension1D dim_2(0.0,1.0,23);
```

Make the grid:
```
Grid grid({dim_1,dim_2});
```

Fill randomly:
```
IdxSet v(2);
for (v[0]=1; v[0]<=30; v[0]++) {
	for (v[1]=1; v[1]<=17; v[1]++) {
		grid.get_grid_point_inside(v)->set_ordinate(fRand(0.0,5.0));
	};
};
```

Write to file:
```
grid.write_to_file("test_deriv_p_2d.txt");
```

Point to evaluate at
```
double* x = new double[2];
x[0] = 0.71;
x[1] = 0.33;
```

Derivs
```
IdxSet local_idxs(2);
double x_deriv;
for (local_idxs[0]=0; local_idxs[0]<=3; local_idxs[0]++) {
	for (local_idxs[1]=0; local_idxs[1]<=3; local_idxs[1]++) {
		x_deriv = grid.get_deriv_wrt_pt_value(x,local_idxs);
		std::cout << "deriv @ " << x[0] << "," << x[1] << " wrt p" << local_idxs[0] << local_idxs[1] << " = " << x_deriv << std::endl;
	};
};
```