# Interpolate in 1D

See `test_interp_1d.cpp` and for validation the Mathematica file `test_interp_1d.nb`.

Make a dimension:
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
grid.write_to_file("test_interp_1d.txt");
```

Get val:
```
double* x = new double[1];
x[0] = 0.71;
double y = grid.get_val(x);
std::cout << "val @ 0.71 = " << y << std::endl;
```