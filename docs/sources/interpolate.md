# Interpolate {#docinterpolate}

This submodule has several related features:

- A module @ref polynomial with basic features to work with polynomials
- A module @ref fitpack to performe interpolation/fitting using an expansion in B-splines
- A module @ref csplines to work with cubic splines, represented as piecewise cubic polynomials


## CubicSplines

 `CubicSplines` gives a frame to interpolate data using Cubic Splines. It uses an internal representation as a list of piece-wise cubic polynomials. It has a functional interface and an Object oriented interface. Both modes of use are equivalent and they rely on the same implementation.
 
### `CubicSpline` Object and methods

The use of the object oriented CubicSplines is very convenient, as it is illustrated in the following example

@include{lineno} ex_interp_csplines_oo.f90

 This outputs:
```{.shell-session}
Value of the function at 1.5     = 0.993493
Value of first derivative at 1.5 = 0.066617
Value of integral from 0 to 1.5  = 0.927102
# There are 7 roots. In units of π:
-0.0000000025
 0.9999594164
 1.9999193799
 2.9998803876
 3.9998428449
 4.9998041860
 5.9995842765
```

and 

```{.shell-session}
$> cat data/ex_interp_cspline1.dat
 # x                sin(x)             evaluate
  0.00000100000000  0.00000100000000  0.00000100000000
  0.25303896202532  0.25034730035496  0.24863008190056
  0.50607692405063  0.48474965603025  0.48261414919792
  0.75911488607595  0.68827961288887  0.68730825131711
  1.01215284810127  0.84797489664940  0.84806743768314
 ...
 18.97784815189873  0.12794059510705  0.12747484922469
 19.23088611392405  0.37215544351995  0.36755693669675
 19.48392407594936  0.59266871498231  0.58580060371012
 19.73696203797468  0.77543651368883  0.77015926631456
 19.99000000000000  0.90881885124069  0.90858634055980

$> cat data/ex_interp_cspline2.dat
 # x                cos(x)             derivative
  0.00000100000000  0.99999999999950  0.99222241991253
  0.25303896202532  0.96815609754057  0.96328391891808
  0.50607692405063  0.87465294316006  0.87646892201067
  0.75911488607595  0.72544550069702  0.73177742919030
  1.01215284810127  0.53003638993227  0.52920944045695
 ...
 18.97784815189873  0.99178183292681  0.97608188013021
 19.23088611392405  0.92817041854310  0.91358086898952
 19.48392407594936  0.80544633234078  0.80347204495850
 19.73696203797468  0.63142554053358  0.64575540803714
 19.99000000000000  0.41719095823083  0.44043095822545
```

Note that in order to get the interpolated values only three lines are relevant:
 - line 3:  where the object is declared.
 - line 20: where the interpolation to the data is determined.
 - line 22: where the interpolated values are obtained.

The evaluation of the data in the new points results in:

@image html ex_interp_cspline1.png

while the derivative gives

@image html ex_interp_cspline2.png

@note that both interpolations are particularly good at both ends because we gave good values for the second derivative.

### `CubicSpline` Functions and subroutines

Interpolation using cubic splines may be accomplished in a very similar manner using the functional interface to `CubicSplines`. Its use is very similar to the object oriented aproach, as shown in the translation of the above example to this interface

@include{lineno} ex_interp_csplines_fp.f90


## Fitpack

This module allows to either fit or interpolate data using an expansion in B-splinebasis of order \f$k \le 5\f$. The implementation is based in the FITPACK set of routines by P. Diercxx.


### Use of `splrep()` with `splev()`

The first example illustrates how to use B-Splines to interpolate a function `y=f(x)` from a few data points:

@include ex_interp_splrep.f90

Prints:

```{.shell-session}
$> cat data/ex_interp_splrep1.dat
 # x                   y
 0.00000000000000  0.00000000000000
 0.11219973762821  0.11401345780438
 0.22439947525641  0.22522911715938
 0.33659921288462  0.33269883705930
 0.44879895051283  0.43547447649845
  ...
 2.69279370307697  0.43547447649845
 2.80499344070517  0.33269883705930
 2.91719317833338  0.22522911715938
 3.02939291596159  0.11401345780438
 3.14159265358979  0.00000000000000
```

Notice that we force interpolation by using `s=0`,

@image html ex_interp_splrep1.png


### Use of `splprep()` with `splev()`
This second example illustrates how to use B-Splines for a parametric curve in the plane.

@include ex_interp_splprep.f90


 Prints:
 
```{.shell-session}
$> cat data/ex_interp_splprep1.dat
 # u                    x                   y
 0.00000000000000  1.50000000000000  0.00000000000000
 0.03616230734577  1.46779335229253  0.23853963732962
 0.07211603380768  1.37398960267532  0.45870512901973
 0.10765440708122  1.22676038619218  0.64385351896871
 0.14257428865429  1.03883011365998  0.78063018793468
 ....
 0.85742571134571  1.03883011365998 -0.78063018793468
 0.89234559291878  1.22676038619218 -0.64385351896871
 0.92788396619232  1.37398960267532 -0.45870512901973
 0.96383769265423  1.46779335229253 -0.23853963732962
 1.00000000000000  1.50000000000000 -0.00000000000000 
```

Notice that:
  1. We force interpolation by using `s=0`,
  2. The parameterization array `u` is generated automatically.

@image html ex_interp_splprep1.png
