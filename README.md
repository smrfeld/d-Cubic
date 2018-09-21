# d-Cubic
Cubic Interpolation in d-Dimensions and Derivatives

## Documentation

See [documentation](doc/cubic_interp_and_derivatives.pdf).

## Installation

Use the makefile: `make` and `make install` to build the shared library. The default location is `/usr/local/lib` for the library and `/usr/local/include` for the headers. 

There is a convenient header `/include/dcubic` such that we can just use `#include <dcubic>` after installation.

## Linking

It is a shared library - use `g++ -std=c++14 -ldcubic my_program.cpp -o my_program.o`.

## Namespace

The namespace is: `dcu`.

## Examples

For code, see the examples folder.

* [Interpolation in 1D](examples/test_interp_1d/)
* [Interpolation in 2D](examples/test_interp_2d/)
* [Derivative with respect to p in 1D](examples/test_deriv_p_1d/)
* [Derivative with respect to p in 2D](examples/test_deriv_p_2d/)
* [Derivative with respect to p in 1D, near the boundary](examples/test_deriv_p_boundary_1d/)
* [Derivative with respect to p in 2D, near the boundary](examples/test_deriv_p_boundary_2d/)
* [Derivative with respect to x in 1D](examples/test_deriv_x_1d/)
* [Derivative with respect to x in 2D](examples/test_deriv_x_2d/)

For some comparable analytic formula, see this [Mathematica file](doc/formula_1d_2d.nb).