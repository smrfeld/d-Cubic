# Take derivative with respect to p in 2D at the boundary

See `test_deriv_p_boundary_2d.cpp` and for validation the Mathematica file `test_deriv_p_boundary_2d.nb`.

Make two dimensions:
```
Dimension1D dim_1(0.0,1.0,11);
Dimension1D dim_2(0.0,1.0,21);
```

Make the grid:
```
Grid grid({&dim_1,&dim_2});
```

Fill randomly
```
std::vector<int> v({0,0});
for (int i=0; i<11; i++) {
	for (int j=0; j<21; j++) {
		v[0] = i;
		v[1] = j;
		grid.get_grid_point_ref(v).set_ordinate(fRand(0.0,5.0));
	};
};
```

Write to file:
```
grid.write_to_file("test_deriv_p_boundary_2d.txt");
```

Point to evaluate at
```
std::vector<double> abcissa({0.023,0.993});
```

NOTE: derivative wrt p0x or px3 will throw an error, since it is out of the grid
```
IdxSet4 local_idxs({0,0});
double x_deriv;
for (int i=1; i<=3; i++) {
	for (int j=0; j<=2; j++) {
		local_idxs[0] = i;
		local_idxs[1] = j;
		x_deriv = grid.get_deriv_wrt_pt_value(abcissa,local_idxs);
		std::cout << "deriv @ " << abcissa[0] << "," << abcissa[1] << " wrt p" << local_idxs[0] << local_idxs[1] << " = " << x_deriv << std::endl;
	};
};
```