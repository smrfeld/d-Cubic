# Take derivative with respect to x in 1D

See `test_deriv_x_1d.cpp` and for validation the Mathematica file `test_deriv_x_1d.nb`.

Make two dimensions:
```
Dimension1D dim(0.0,1.0,30);
```

Make the grid:
```
Grid grid({dim});
```

Fill randomly:
```
IdxSet v(1);
for (v[0]=1; v[0]<=30; v[0]++) {
	grid.get_grid_point_inside(v)->set_ordinate(fRand(-4.0,-2.0));
};
```

Write to file:
```
grid.write_to_file("test_deriv_x_1d.txt");
```

Point to evaluate at
```
double* x = new double[1];
x[0] = 0.46;
```

Derivs
```
double x_deriv = grid.get_deriv_wrt_x(x,0);
std::cout << "deriv @ " << x[0] << " wrt x = " << x_deriv << std::endl;
```
