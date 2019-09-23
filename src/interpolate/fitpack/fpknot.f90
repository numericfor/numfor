subroutine fpknot(x, m, t, n, fpint, nrdata, nrint, nest, istart)
  implicit none
!  subroutine fpknot locates an additional knot for a spline of degree
!  k and adjusts the corresponding parameters,i.e.
!    t     : the position of the knots.
!    n     : the number of knots.
!    nrint : the number of knotintervals.
!    fpint : the sum of squares of residual right hand sides
!            for each knot interval.
!    nrdata: the number of data points inside each knot interval.
!  istart indicates that the smallest data point at which the new knot
!  may be added is x(istart+1)
!  ..
!  ..scalar arguments..
  integer :: m, n, nrint, nest, istart
!  ..array arguments..
  real(8) :: x(m), t(nest), fpint(nest)
  integer :: nrdata(nest)
!  ..local scalars..
  real(8) :: an, am, fpmax
  integer :: ihalf, j, jbegin, jj, jk, jpoint, k, maxbeg, maxpt, next, nrx, number
!  ..
  k = (n - nrint - 1) / 2
!  search for knot interval t(number+k) <= x <= t(number+k+1) where
!  fpint(number) is maximal on the condition that nrdata(number)
!  not equals zero.
  fpmax = 0.
  jbegin = istart
  do j = 1, nrint
    jpoint = nrdata(j)
    if (fpmax < fpint(j) .and. jpoint /= 0) then
      fpmax = fpint(j)
      number = j
      maxpt = jpoint
      maxbeg = jbegin
    end if
    jbegin = jbegin + jpoint + 1
  end do

!  let coincide the new knot t(number+k+1) with a data point x(nrx)
!  inside the old knot interval t(number+k) <= x <= t(number+k+1).
  ihalf = maxpt / 2 + 1
  nrx = maxbeg + ihalf
  next = number + 1
  if (next <= nrint) then
!  adjust the different parameters.
    do j = next, nrint
      jj = next + nrint - j
      fpint(jj + 1) = fpint(jj)
      nrdata(jj + 1) = nrdata(jj)
      jk = jj + k
      t(jk + 1) = t(jk)
    end do
  end if

  nrdata(number) = ihalf - 1
  nrdata(next) = maxpt - ihalf
  am = maxpt
  an = nrdata(number)
  fpint(number) = fpmax * an / am
  an = nrdata(next)
  fpint(next) = fpmax * an / am
  jk = next + k
  t(jk) = x(nrx)
  n = n + 1
  nrint = nrint + 1
  return
end
