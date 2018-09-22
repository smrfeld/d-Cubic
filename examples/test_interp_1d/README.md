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
std::vector<int> v({0});
for (int i=0; i<30; i++) {
	v[0] = i;
	grid.get_grid_point_ref(v).set_ordinate(fRand(0.0,5.0));
};
```

Write to file:
```
grid.write_to_file("test_interp_1d.txt");
```

Get val:
```
double x = grid.get_val({0.71});
std::cout << "val @ 0.71 = " << x << std::endl;
```