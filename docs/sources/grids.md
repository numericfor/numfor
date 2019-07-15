# Submodule grids {#modgrids}

This module provides convenience routines to create commonly occurring grids, somewhat mimicking those appearing in `Numpy`:

  * `linspace()` returns evenly (linearly) spaced numbers over a specified interval.
  * `logspace()` returns logarithmically evenly spaced numbers over a specified interval.
  * `arange()` returns an array of **integer** numbers from a given interval.

## Details

The signature of `linspace` is:

```
linspace(start, end, num = 50[, endpoint=.True.] [, retstep] )
```

This means (compare to documentation of `numpy.linspace`):
  * The name of the function is `linspace`
  * `start`: The first value of the sequence desired
  * `end`: The end value of the sequence, unless endpoint is set to False. In that case, the sequence consists of all but the last of num + 1 evenly spaced samples, so that stop is excluded. Note that the step size changes when endpoint is False.
  * `num`: Number of samples to generate. **Note**: it is required **not optional**.
  * `endpoint`: (boolean, optional). If .True., `end` is the last sample. Otherwise, it is not included. Default is .True.
  * `retstep`: (optional). If present will have the value of `step` on output.


The signature of `logspace` is:

```
logspace(start, end, num = 50[, endpoint=.True.] [, base=10._dp] )
```

where:
  * The interval spans between `base**start` and `base**end`
  * `endpoint` indicates if the final point `end` will be included, like in `linspace`.
  * `base` of the log space. By default is 10.
  
# Examples of use

## Examples of linspace


```{.f90}
  use basic, only: dp, Zero
  use grids, only: linspace
  implicit none

  integer, parameter :: N = 5
  real(dp), dimension(N) :: a
  real(dp), dimension(:), allocatable :: b
  real(dp) :: step
  real(dp) :: x1, x2

  x1 = 1._dp
  x2 = 3._dp

  a = linspace(Zero, x1, N, endpoint=.False., retstep=step)
  b = linspace(x1, x2, nint((x2 - x1) / step) + 1)
  print *, 'size(a)', size(a)
  print '(5(f5.2,1x))', a
  print '(A)', repeat('-', 65)
  print *, 'size(b)', size(b)
  print '(11(f5.2,1x))', b

```
which gives
```{.sh}
 size(a)           5
 0.00  0.20  0.40  0.60  0.80
-----------------------------------------------------------------
 size(b)          11
 1.00  1.20  1.40  1.60  1.80  2.00  2.20  2.40  2.60  2.80  3.00
 ```

Note here that the first call to `linspace()` gives `N` points between 0 and 1, not including the last number. Also gives in the variable `step` the value of the stepsize used. The second call uses a number of points, such that the step is preserved.

## Examples of logspace


```{.f90}
  use basic, only: dp, Zero
  use grids
  implicit none

  integer, parameter :: N = 10
  real(dp), dimension(2 * N) :: a
  real(dp) :: x1, x2

  x1 = -3._dp
  x2 = 3._dp

  a(:N) = logspace(x1, Zero, N, endpoint=.False.)
  a(N:) = linspace(1._dp, x2, N)
  print *, 'size(a)', size(a)
  print '(5(g8.2,1x))', a
```
