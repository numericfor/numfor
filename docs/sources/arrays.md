# Module Arrays {#docarrays}

This module contains routines to generate and operate on arrays

  * Grid creation: There are a few functions to quickly create arrays. They are described in @subpage modgrids.
  * Utilities for array manipulation and information, described in @subpage modarrut
  * Search and sort: Functions for search an element in a sorted array, and sort arrays.
  * Generation of histograms. Do we want histogram2d?

## Grids

This submodule provides convenience routines to create commonly occurring grids, somewhat mimicking those appearing in `Numpy`:

  * `linspace()` returns evenly (linearly) spaced numbers over a specified interval.
  * `logspace()` returns logarithmically evenly spaced numbers over a specified interval.
  * `arange()` returns an array of **integer** numbers from a given interval.

### Details ###

The signature of `linspace` is:

```
linspace(start, end, num[, endpoint=.True.] [, retstep] )
```

This means (compare to documentation of `numpy.linspace`):
  * The name of the function is `linspace`
  * `start`: The first value of the sequence desired
  * `end`: The end value of the sequence, unless endpoint is set to False. In that case, the sequence consists of all but the last of num + 1 evenly spaced samples, so that stop is excluded. Note that the step size changes when endpoint is False.
  * `num`: Number of samples to generate. **Note**: it is required **not optional**.
  * `endpoint`: (boolean, optional). If .True., `end` is the last sample. Otherwise, it is not included. Default is .True.
  * `retstep`: (optional). If present will have the value of `step` on output.

@include ex_linspace.f90

Prints

```{.shell-session}
2.00 2.25 2.50 2.75 3.00
2.00 2.20 2.40 2.60 2.80
2.00 2.25 2.50 2.75 3.00
--------------------
 0.00000  0.20000  0.40000  0.70000  1.10000
 0.05000  0.25000  0.45000  0.80000  1.20000
 0.10000  0.30000  0.50000  0.90000  1.30000
 0.15000  0.35000  0.60000  1.00000  1.40000
``` 

@image html ex_linspace.png

The signature of `logspace` is:

```
logspace(start, end, num = 50[, endpoint=.True.] [, base=10._dp] )
```

where:
  * The interval spans between `base**start` and `base**end`
  * `endpoint` indicates if the final point `end` will be included, like in `linspace`.
  * `base` of the log space. By default is 10.
  
@include ex_logspace.f90

Prints
```{.shell-session}
  100.00   215.44   464.16  1000.00
    4.00     5.04     6.35     8.00
 1000.00   100.00    10.00     1.00
--------------------
    1.00    10.00   100.00  1000.00
-1000.00  -100.00   -10.00    -1.00
--------------------
 0.000010  0.000183  0.003360  0.061585  1.128838
 0.000021  0.000379  0.006952  0.127427  2.335721
 0.000043  0.000785  0.014384  0.263665  4.832930
 0.000089  0.001624  0.029764  0.545559 10.000000
``` 

@image html ex_logspace.png


