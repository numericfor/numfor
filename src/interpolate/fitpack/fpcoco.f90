subroutine fpcoco(iopt, m, x, y, w, v, s, nest, maxtr, maxbin, n, t, c, sq, sx, bind, e, &
        &wrk, lwrk, iwrk, kwrk, ier)
!  ..scalar arguments..
  real(8) :: s, sq
  integer :: iopt, m, nest, maxtr, maxbin, n, lwrk, kwrk, ier
!  ..array arguments..
  integer :: iwrk(kwrk)
  real(8) :: x(m), y(m), w(m), v(m), t(nest), c(nest), sx(m), e(nest), wrk(lwrk)
  logical :: bind(nest)
!  ..local scalars..
  integer :: i, ia, ib, ic, iq, iu, iz, izz, i1, j, k, l, l1, m1, nmax, nr, n4, n6, n8, &
    &ji, jib, jjb, jl, jr, ju, mb, nm
  real(8) :: sql, sqmax, term, tj, xi, half
!  ..subroutine references..
!    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno
!  ..
!  set constant
  half = 0.5e0
!  determine the maximal admissible number of knots.
  nmax = m + 4
!  the initial choice of knots depends on the value of iopt.
!    if iopt=0 the program starts with the minimal number of knots
!    so that can be guarantied that the concavity/convexity constraints
!    will be satisfied.
!    if iopt = 1 the program will continue from the point on where she
!    left at the foregoing call.
  if (iopt > 0) go to 80
!  find the minimal number of knots.
!  a knot is located at the data point x(i), i=2,3,...m-1 if
!    1) v(i) ^= 0    and
!    2) v(i)*v(i-1) <= 0  or  v(i)*v(i+1) <= 0.
  m1 = m - 1
  n = 4
  do i = 2, m1
    if (v(i) == 0._8 .or. (v(i) * v(i - 1) > 0._8 .and. (v(i) * v(i + 1) > 0._8))) cycle
    n = n + 1
!  test whether the required storage space exceeds the available one.
    if (n + 4 > nest) then
      ier = 4
      return
    end if

    t(n) = x(i)
  end do

!  find the position of the knots t(1),...t(4) and t(n-3),...t(n) which
!  are needed for the b-spline representation of s(x).
  do i = 1, 4
    t(i) = x(1)
    n = n + 1
    t(n) = x(m)
  end do

!  test whether the minimum number of knots exceeds the maximum number.
  if (n > nmax) then
    ier = 5
    return
  end if

!  main loop for the different sets of knots.
!  find corresponding values e(j) to the knots t(j+3),j=1,2,...n-6
!    e(j) will take the value -1,1, or 0 according to the requirement
!    that s(x) must be locally convex or concave at t(j+3) or that the
!    sign of s''(x) is unrestricted at that point.
40 i = 1
  xi = x(1)
  j = 4
  tj = t(4)
  n6 = n - 6
  do l = 1, n6
    do while (xi /= tj)
      i = i + 1
      xi = x(i)
    end do

    e(l) = v(i)
    j = j + 1
    tj = t(j)
  end do

!  we partition the working space
  nm = n + maxbin
  mb = maxbin + 1
  ia = 1
  ib = ia + 4 * n
  ic = ib + nm * maxbin
  iz = ic + n
  izz = iz + n
  iu = izz + n
  iq = iu + maxbin
  ji = 1
  ju = ji + maxtr
  jl = ju + maxtr
  jr = jl + maxtr
  jjb = jr + maxtr
  jib = jjb + mb
!  given the set of knots t(j),j=1,2,...n, find the least-squares cubic
!  spline which satisfies the imposed concavity/convexity constraints.
  call fpcosp(m, x, y, w, n, t, e, maxtr, maxbin, c, sq, sx, bind, nm, mb, wrk(ia),&
    &wrk(ib), wrk(ic), wrk(iz), wrk(izz), wrk(iu), wrk(iq), iwrk(ji), iwrk(ju), &
    &iwrk(jl), iwrk(jr), iwrk(jjb), iwrk(jib), ier)
!  if sq <= s or in case of abnormal exit from fpcosp, control is
!  repassed to the driver program.
  if (sq <= s .or. ier > 0) return
!  calculate for each knot interval t(l-1) <= xi <= t(l) the
!  sum((wi*(yi-s(xi)))**2).
!  find the interval t(k-1) <= x <= t(k) for which this sum is maximal
!  on the condition that this interval contains at least one interior
!  data point x(nr) and that s(x) is not given there by a straight line.
80 sqmax = 0.
  sql = 0.
  l = 5
  nr = 0
  i1 = 1
  n4 = n - 4
  do i = 1, m
    term = (w(i) * (sx(i) - y(i)))**2
    if (x(i) < t(l) .or. l > n4) go to 100
    term = term * half
    sql = sql + term
    if (i - i1 <= 1 .or. (bind(l - 4) .and. bind(l - 3))) go to 90
    if (sql <= sqmax) go to 90
    k = l
    sqmax = sql
    nr = i1 + (i - i1) / 2
90  l = l + 1
    i1 = i
    sql = 0.
100 sql = sql + term
  end do

  if (m - i1 <= 1 .or. (bind(l - 4) .and. bind(l - 3))) go to 120
  if (sql <= sqmax) go to 120
  k = l
  nr = i1 + (m - i1) / 2
!  if no such interval is found, control is repassed to the driver
!  program (ier = -1).
120 if (nr == 0) then
    ier = -1
    return
  end if

!  if s(x) is given by the same straight line in two succeeding knot
!  intervals t(l-1) <= x <= t(l) and t(l) <= x <= t(l+1),delete t(l)
  n8 = n - 8
  l1 = 0
  if (n8 <= 0) go to 150
  do i = 1, n8
    if (.not. (bind(i) .and. bind(i + 1) .and. bind(i + 2))) cycle
    l = i + 4 - l1
    if (k > l) k = k - 1
    n = n - 1
    l1 = l1 + 1
    do j = l, n
      t(j) = t(j + 1)
    end do
  end do

!  test whether we cannot further increase the number of knots.
150 if (n == nmax) then
    ier = -2
    return
  end if

  if (n == nest) then
    ier = -3
    return
  end if

!  locate an additional knot at the point x(nr).
  j = n
  do i = k, n
    t(j + 1) = t(j)
    j = j - 1
  end do

  t(k) = x(nr)
  n = n + 1
!  restart the computations with the new set of knots.
  go to 40
!  error codes and messages.
end subroutine fpcoco

