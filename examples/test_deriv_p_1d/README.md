# Take derivative with respect to p in 1D

See `test_deriv_p_1d.cpp` and for validation the Mathematica file `test_deriv_p_1d.nb`.

Make two dimensions:
```
Dimension1D dim(0.0,1.0,30);
```

Make the grid:
```
Grid grid({&dim});
```

Fill randomly:
```
IdxSet v(1);
for (v[0]=1; v[0]<=30; v[0]++) {
	grid.get_grid_point_inside(v)->set_ordinate(fRand(0.0,5.0));
};
```

Write to file:
```
grid.write_to_file("test_deriv_p_1d.txt");
```

Point to evaluate at
```
double* x = new double[1];
x[0] = 0.71;
```

Derivs
```
IdxSet local_idxs(1);
double x_deriv;
for (local_idxs[0]=0; local_idxs[0]<=3; local_idxs[0]++) {
	x_deriv = grid.get_deriv_wrt_pt_value(x,local_idxs);
	std::cout << "deriv @ " << x[0] << " wrt p" << local_idxs[0] << " = " << x_deriv << std::endl;
};
```