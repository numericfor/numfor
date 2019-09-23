subroutine sproot(t, n, c, zero, mest, m, ier)
  !  subroutine sproot finds the zeros of a cubic spline s(x),which is
  !  given in its normalized b-spline representation.
  !
  !  calling sequence:
  !     call sproot(t,n,c,zero,mest,m,ier)
  !
  !  input parameters:
  !    t    : real array,length n, containing the knots of s(x).
  !    n    : integer, containing the number of knots.  n>=8
  !    c    : real array,length n, containing the b-spline coefficients.
  !    mest : integer, specifying the dimension of array zero.
  !
  !  output parameters:
  !    zero : real array,length mest, containing the zeros of s(x).
  !    m    : integer,giving the number of zeros.
  !    ier  : error flag:
  !      ier = 0: normal return.
  !      ier = 1: the number of zeros exceeds mest.
  !      ier =10: invalid input data (see restrictions).
  !
  !  other subroutines required: fpcuro
  !
  !  restrictions:
  !    1) n>= 8.
  !    2) t(4) < t(5) < ... < t(n-4) < t(n-3).
  !       t(1) <= t(2) <= t(3) <= t(4)
  !       t(n-3) <= t(n-2) <= t(n-1) <= t(n)
  !
  !  author :
  !    p.dierckx
  !    dept. computer science, k.u.leuven
  !    celestijnenlaan 200a, b-3001 heverlee, belgium.
  !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
  !
  !  latest update : march 1987
  !
  ! ..
  ! ..scalar arguments..
  integer :: n, mest, m, ier
  !  ..array arguments..
  real(8) :: t(n), c(n), zero(mest)
  !  ..local scalars..
  integer :: i, j, j1, l, n4
  real(8) :: ah, a0, a1, a2, a3, bh, b0, b1, c1, c2, c3, c4, c5, d4, d5, h1, h2,&
    &three, two, t1, t2, t3, t4, t5, zz
  logical z0, z1, z2, z3, z4, nz0, nz1, nz2, nz3, nz4
  !  ..local array..
  real(8) :: y(3)
  !  ..
  !  set some constants
  two = 0.2d+01
  three = 0.3d+01
  !  before starting computations a data check is made. if the input data
  !  are invalid, control is immediately repassed to the calling program.
  n4 = n - 4
  ier = 10
  if (n < 8) go to 800
  j = n
  do i = 1, 3
    if (t(i) > t(i + 1)) go to 800
    if (t(j) < t(j - 1)) go to 800
    j = j - 1
  end do

  do i = 4, n4
    if (t(i) >= t(i + 1)) go to 800
  end do

  !  the problem considered reduces to finding the zeros of the cubic
  !  polynomials pl(x) which define the cubic spline in each knot
  !  interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
  !  the condition that it belongs to the knot interval.
  !  the cubic polynomial pl(x) is determined by computing s(t(l)),
  !  s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
  !  s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
  !  splines and their derivatives, the value of s(t(l)) and s'(t(l))
  !  is already known from the foregoing knot interval.
  ier = 0
  !  evaluate some constants for the first knot interval
  h1 = t(4) - t(3)
  h2 = t(5) - t(4)
  t1 = t(4) - t(2)
  t2 = t(5) - t(3)
  t3 = t(6) - t(4)
  t4 = t(5) - t(2)
  t5 = t(6) - t(3)
  !  calculate a0 = s(t(4)) and ah = s'(t(4)).
  c1 = c(1)
  c2 = c(2)
  c3 = c(3)
  c4 = (c2 - c1) / t4
  c5 = (c3 - c2) / t5
  d4 = (h2 * c1 + t1 * c2) / t4
  d5 = (t3 * c2 + h1 * c3) / t5
  a0 = (h2 * d4 + h1 * d5) / t2
  ah = three * (h2 * c4 + h1 * c5) / t2
  z1 = .true.
  if (ah < 0.0d0) z1 = .false.
  nz1 = .not. z1
  m = 0
  !  main loop for the different knot intervals.
  do l = 4, n4
    !  evaluate some constants for the knot interval t(l) <= x <= t(l+1).
    h1 = h2
    h2 = t(l + 2) - t(l + 1)
    t1 = t2
    t2 = t3
    t3 = t(l + 3) - t(l + 1)
    t4 = t5
    t5 = t(l + 3) - t(l)
    !  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
    c1 = c2
    c2 = c3
    c3 = c(l)
    c4 = c5
    c5 = (c3 - c2) / t5
    d4 = (h2 * c1 + t1 * c2) / t4
    d5 = (h1 * c3 + t3 * c2) / t5
    b0 = (h2 * d4 + h1 * d5) / t2
    bh = three * (h2 * c4 + h1 * c5) / t2
    !  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
    !  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
    a1 = ah * h1
    b1 = bh * h1
    a2 = three * (b0 - a0) - b1 - two * a1
    a3 = two * (a0 - b0) + b1 + a1
    !  test whether or not pl(x) could have a zero in the range
    !  t(l) <= x <= t(l+1).
    z3 = .true.
    if (b1 < 0.0d0) z3 = .false.
    nz3 = .not. z3
    if (a0 * b0 <= 0.0d0) go to 100
    z0 = .true.
    if (a0 < 0.0d0) z0 = .false.
    nz0 = .not. z0
    z2 = .true.
    if (a2 < 0.) z2 = .false.
    nz2 = .not. z2
    z4 = .true.
    if (3.0d0 * a3 + a2 < 0.0d0) z4 = .false.
    nz4 = .not. z4
    if (.not. ((z0 .and. (nz1 .and. (z3 .or. z2 .and. nz4) .or. nz2 .and.&
      &z3 .and. z4) .or. nz0 .and. (z1 .and. (nz3 .or. nz2 .and. z4) .or. z2 .and.&
      &nz3 .and. nz4)))) go to 200
    !  find the zeros of ql(y).
100 call fpcuro(a3, a2, a1, a0, y, j)
    if (j == 0) go to 200
    !  find which zeros of pl(x) are zeros of s(x).
    do i = 1, j
      if (y(i) < 0.0d0 .or. y(i) > 1.0d0) go to 150
      !  test whether the number of zeros of s(x) exceeds mest.
      if (m >= mest) go to 700
      m = m + 1
      zero(m) = t(l) + h1 * y(i)
150   continue
    end do

200 a0 = b0
    ah = bh
    z1 = z3
    nz1 = nz3
  end do

  !  the zeros of s(x) are arranged in increasing order.
  if (m < 2) go to 800
  do i = 2, m
    j = i
350 j1 = j - 1
    if (j1 == 0) go to 400
    if (zero(j) >= zero(j1)) go to 400
    zz = zero(j)
    zero(j) = zero(j1)
    zero(j1) = zz
    j = j1
    go to 350
400 continue
  end do

  j = m
  m = 1
  do i = 2, j
    if (zero(i) == zero(m)) go to 500
    m = m + 1
    zero(m) = zero(i)
500 continue
  end do

  go to 800
700 ier = 1
800 return
end subroutine sproot
