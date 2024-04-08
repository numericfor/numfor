# Submodule random {#submodule-random}

This module contains routines for

  * Random [number generation](@ref randomnumber)
  * Random [distributions](@ref randist)
  
## @anchor randomnumber Random Number generation

There are implemented random number generators using the 64-bit version of MT19937
originally coded by Takuji Nishimura and Makoto Matsumoto (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html) and the Fortran translation by Rémi Piatek (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/fortran.html)

The random number generator may be initialized with `seed()` and then a random number is generated:

@include ex_random.f90

## @anchor randist Random Distributions

Implementations of random distributions. Currently are coded:

- @ref uniform Random distributions
- @ref gauss distributions and 
- @ref exponential distributions 
