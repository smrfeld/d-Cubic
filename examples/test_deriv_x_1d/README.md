# Take derivative with respect to x in 1D

See `test_deriv_x_1d.cpp` and for validation the Mathematica file `test_deriv_x_1d.nb`.

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
std::vector<int> v({0});
for (int i=0; i<30; i++) {
	v[0] = i;
	grid.get_grid_point_ref(v).set_ordinate(fRand(-4.0,-2.0));
};
```

Write to file:
```
grid.write_to_file("test_deriv_x_1d.txt");
```

Point to evaluate at
```
std::vector<double> abcissa({0.46});
```

Derivs
```
double x_deriv = grid.get_deriv_wrt_x(abcissa,0);
std::cout << "deriv @ " << abcissa[0] << " wrt x = " << x_deriv << std::endl;
```
