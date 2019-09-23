subroutine fpsurf(iopt, m, x, y, z, w, xb, xe, yb, ye, kxx, kyy, s, nxest,&
  &nyest, eta, tol, maxit, nmax, km1, km2, ib1, ib3, nc, intest, nrest,&
  &nx0, tx, ny0, ty, c, fp, fp0, fpint, coord, f, ff, a, q, bx, by, spx, spy, h,&
  &index, nummer, wrk, lwrk, ier)
!  ..
!  ..scalar arguments..
  real(8) :: xb, xe, yb, ye, s, eta, tol, fp, fp0
  integer :: iopt, m, kxx, kyy, nxest, nyest, maxit, nmax, km1, km2, ib1, ib3,&
    &nc, intest, nrest, nx0, ny0, lwrk, ier
!  ..array arguments..
  real(8) :: x(m), y(m), z(m), w(m), tx(nmax), ty(nmax), c(nc), fpint(intest),&
    &coord(intest), f(nc), ff(nc), a(nc, ib1), q(nc, ib3), bx(nmax, km2),&
    &by(nmax, km2), spx(m, km1), spy(m, km1), h(ib3), wrk(lwrk)
  integer :: index(nrest), nummer(m)
!  ..local scalars..
  real(8) :: acc, arg, cos, dmax, fac1, fac2, fpmax, fpms, f1, f2, f3, hxi, p, pinv,&
    &piv, p1, p2, p3, sigma, sin, sq, store, wi, x0, x1, y0, y1, zi, eps,&
    &rn, one, con1, con9, con4, half, ten
  integer :: i, iband, iband1, iband3, iband4, ibb, ichang, ich1, ich3, ii,&
    &in, irot, iter, i1, i2, i3, j, jrot, jxy, j1, kx, kx1, kx2, ky, ky1, ky2, l,&
    &la, lf, lh, lwest, lx, ly, l1, l2, n, ncof, nk1x, nk1y, nminx, nminy, nreg,&
    &nrint, num, num1, nx, nxe, nxx, ny, nye, nyy, n1, rank
!  ..local arrays..
  real(8) :: hx(6), hy(6)
!  ..function references..
  real(8) :: abs, fprati, sqrt
  integer :: min0
!  ..subroutine references..
!    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota
!  ..
!  set constants
  one = 0.1e+01
  con1 = 0.1e0
  con9 = 0.9e0
  con4 = 0.4e-01
  half = 0.5e0
  ten = 0.1e+02
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 1: determination of the number of knots and their position.     c
! ****************************************************************     c
! given a set of knots we compute the least-squares spline sinf(x,y),  c
! and the corresponding weighted sum of squared residuals fp=f(p=inf). c
! if iopt=-1  sinf(x,y) is the requested approximation.                c
! if iopt=0 or iopt=1 we check whether we can accept the knots:        c
!   if fp <=s we will continue with the current set of knots.          c
!   if fp > s we will increase the number of knots and compute the     c
!      corresponding least-squares spline until finally  fp<=s.        c
! the initial choice of knots depends on the value of s and iopt.      c
!   if iopt=0 we first compute the least-squares polynomial of degree  c
!     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c
!     fp0=f(0) denotes the corresponding weighted sum of squared       c
!     residuals                                                        c
!   if iopt=1 we start with the knots found at the last call of the    c
!     routine, except for the case that s>=fp0; then we can compute    c
!     the least-squares polynomial directly.                           c
! eventually the independent variables x and y (and the corresponding  c
! parameters) will be switched if this can reduce the bandwidth of the c
! system to be solved.                                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ichang denotes whether(1) or not(-1) the directions have been inter-
!  changed.
  ichang = -1
  x0 = xb
  x1 = xe
  y0 = yb
  y1 = ye
  kx = kxx
  ky = kyy
  kx1 = kx + 1
  ky1 = ky + 1
  nxe = nxest
  nye = nyest
  eps = sqrt(eta)
  if (iopt < 0) go to 20
!  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  if (iopt == 0) go to 10
  if (fp0 > s) go to 20
!  initialization for the least-squares polynomial.
10 nminx = 2 * kx1
  nminy = 2 * ky1
  nx = nminx
  ny = nminy
  ier = -2
  go to 30
20 nx = nx0
  ny = ny0
!  main loop for the different sets of knots. m is a save upper bound
!  for the number of trials.
30 do iter = 1, m
!  find the position of the additional knots which are needed for the
!  b-spline representation of s(x,y).
    l = nx
    do i = 1, kx1
      tx(i) = x0
      tx(l) = x1
      l = l - 1
    end do

    l = ny
    do i = 1, ky1
      ty(i) = y0
      ty(l) = y1
      l = l - 1
    end do

!  find nrint, the total number of knot intervals and nreg, the number
!  of panels in which the approximation domain is subdivided by the
!  intersection of knots.
    nxx = nx - 2 * kx1 + 1
    nyy = ny - 2 * ky1 + 1
    nrint = nxx + nyy
    nreg = nxx * nyy
!  find the bandwidth of the observation matrix a.
!  if necessary, interchange the variables x and y, in order to obtain
!  a minimal bandwidth.
    iband1 = kx * (ny - ky1) + ky
    l = ky * (nx - kx1) + kx
    if (iband1 <= l) go to 130
    iband1 = l
    ichang = -ichang
    do i = 1, m
      store = x(i)
      x(i) = y(i)
      y(i) = store
    end do

    store = x0
    x0 = y0
    y0 = store
    store = x1
    x1 = y1
    y1 = store
    n = min0(nx, ny)
    do i = 1, n
      store = tx(i)
      tx(i) = ty(i)
      ty(i) = store
    end do

    n1 = n + 1
    if (nx < ny) go to 80
    if (nx == ny) go to 120
    go to 100
80  tx(n1:ny) = ty(n1:ny)

    go to 120
100 ty(n1:nx) = tx(n1:nx)

120 l = nx
    nx = ny
    ny = l
    l = nxe
    nxe = nye
    nye = l
    l = nxx
    nxx = nyy
    nyy = l
    l = kx
    kx = ky
    ky = l
    kx1 = kx + 1
    ky1 = ky + 1
130 iband = iband1 + 1
!  arrange the data points according to the panel they belong to.
    call fporde(x, y, m, kx, ky, tx, nx, ty, ny, nummer, index, nreg)
!  find ncof, the number of b-spline coefficients.
    nk1x = nx - kx1
    nk1y = ny - ky1
    ncof = nk1x * nk1y
!  initialize the observation matrix a.

    f(:ncof) = 0.

    a(:ncof, :iband) = 0.

!  initialize the sum of squared residuals.
    fp = 0.
!  fetch the data points in the new order. main loop for the
!  different panels.
    do num = 1, nreg
!  fix certain constants for the current panel; jrot records the column
!  number of the first non-zero element in a row of the observation
!  matrix according to a data point of the panel.
      num1 = num - 1
      lx = num1 / nyy
      l1 = lx + kx1
      ly = num1 - lx * nyy
      l2 = ly + ky1
      jrot = lx * nk1y + ly
!  test whether there are still data points in the panel.
      in = index(num)
150   if (in == 0) go to 250
!  fetch a new data point.
      wi = w(in)
      zi = z(in) * wi
!  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
      call fpbspl(tx, nx, kx, x(in), l1, hx)
!  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
      call fpbspl(ty, ny, ky, y(in), l2, hy)
!  store the value of these b-splines in spx and spy respectively.
      spx(in, :kx1) = hx(:kx1)
      spy(in, :ky1) = hy(:ky1)

!  initialize the new row of observation matrix.
      h(:iband) = 0._8

!  calculate the non-zero elements of the new row by making the cross
!  products of the non-zero b-splines in x- and y-direction.
      i1 = 0
      do i = 1, kx1
        hxi = hx(i)
        j1 = i1
        do j = 1, ky1
          j1 = j1 + 1
          h(j1) = hxi * hy(j) * wi
        end do

        i1 = i1 + nk1y
      end do

!  rotate the row into triangle by givens transformations .
      irot = jrot
      do i = 1, iband
        irot = irot + 1
        piv = h(i)
        if (piv == 0.) cycle
!  calculate the parameters of the givens transformation.
        call fpgivs(piv, a(irot, 1), cos, sin)
!  apply that transformation to the right hand side.
        call fprota(cos, sin, zi, f(irot))
        if (i == iband) go to 230
!  apply that transformation to the left hand side.
        i2 = 1
        i3 = i + 1
        do j = i3, iband
          i2 = i2 + 1
          call fprota(cos, sin, h(j), a(irot, i2))
        end do

      end do

!  add the contribution of the row to the sum of squares of residual
!  right hand sides.
230   fp = fp + zi**2
!  find the number of the next data point in the panel.
      in = nummer(in)
      go to 150
250   continue
    end do

!  find dmax, the maximum value for the diagonal elements in the reduced
!  triangle.
    dmax = 0.
    do i = 1, ncof
      if (a(i, 1) <= dmax) cycle
      dmax = a(i, 1)
    end do

!  check whether the observation matrix is rank deficient.
    sigma = eps * dmax
    do i = 1, ncof
      if (a(i, 1) <= sigma) go to 280
    end do

!  backward substitution in case of full rank.
    call fpback(a, f, ncof, iband, c, nc)
    rank = ncof

    q(:ncof, 1) = a(:ncof, 1) / dmax

    go to 300
!  in case of rank deficiency, find the minimum norm solution.
!  check whether there is sufficient working space
280 lwest = ncof * iband + ncof + iband
    if (lwrk < lwest) go to 780

    ff(:ncof) = f(:ncof)
    q(:ncof, :iband) = a(:ncof, :iband)

    lf = 1
    lh = lf + ncof
    la = lh + iband
    call fprank(q, ff, ncof, iband, nc, sigma, c, sq, rank, wrk(la),&
      &wrk(lf), wrk(lh))

    q(:ncof, 1) = q(:ncof, 1) / dmax

!  add to the sum of squared residuals, the contribution of reducing
!  the rank.
    fp = fp + sq
300 if (ier == (-2)) fp0 = fp
!  test whether the least-squares spline is an acceptable solution.
    if (iopt < 0) go to 820
    fpms = fp - s
    if (abs(fpms) <= acc) then
      if (fp <= 0) go to 815
      go to 820
    endif
!  test whether we can accept the choice of knots.
    if (fpms < 0.) go to 430
!  test whether we cannot further increase the number of knots.
    if (ncof > m) go to 790
    ier = 0
!  search where to add a new knot.
!  find for each interval the sum of squared residuals fpint for the
!  data points having the coordinate belonging to that knot interval.
!  calculate also coord which is the same sum, weighted by the position
!  of the data points considered.

    fpint(:nrint) = 0.
    coord(:nrint) = 0.

    do num = 1, nreg
      num1 = num - 1
      lx = num1 / nyy
      l1 = lx + 1
      ly = num1 - lx * nyy
      l2 = ly + 1 + nxx
      jrot = lx * nk1y + ly
      in = index(num)
330   if (in == 0) cycle
      store = 0.
      i1 = jrot
      do i = 1, kx1
        hxi = spx(in, i)
        j1 = i1
        do j = 1, ky1
          j1 = j1 + 1
          store = store + hxi * spy(in, j) * c(j1)
        end do

        i1 = i1 + nk1y
      end do

      store = (w(in) * (z(in) - store))**2
      fpint(l1) = fpint(l1) + store
      coord(l1) = coord(l1) + store * x(in)
      fpint(l2) = fpint(l2) + store
      coord(l2) = coord(l2) + store * y(in)
      in = nummer(in)
      go to 330
    end do

!  find the interval for which fpint is maximal on the condition that
!  there still can be added a knot.
370 l = 0
    fpmax = 0.
    l1 = 1
    l2 = nrint
    if (nx == nxe) l1 = nxx + 1
    if (ny == nye) l2 = nxx
    if (l1 > l2) go to 810
    do i = l1, l2
      if (fpmax >= fpint(i)) cycle
      l = i
      fpmax = fpint(i)
    end do

!  test whether we cannot further increase the number of knots.
    if (l == 0) go to 785
!  calculate the position of the new knot.
    arg = coord(l) / fpint(l)
!  test in what direction the new knot is going to be added.
    if (l > nxx) go to 400
!  addition in the x-direction.
    jxy = l + kx1
    fpint(l) = 0.
    fac1 = tx(jxy) - arg
    fac2 = arg - tx(jxy - 1)
    if (fac1 > (ten * fac2) .or. fac2 > (ten * fac1)) go to 370
    j = nx
    do i = jxy, nx
      tx(j + 1) = tx(j)
      j = j - 1
    end do

    tx(jxy) = arg
    nx = nx + 1
    cycle
!  addition in the y-direction.
400 jxy = l + ky1 - nxx
    fpint(l) = 0.
    fac1 = ty(jxy) - arg
    fac2 = arg - ty(jxy - 1)
    if (fac1 > (ten * fac2) .or. fac2 > (ten * fac1)) go to 370
    j = ny
    do i = jxy, ny
      ty(j + 1) = ty(j)
      j = j - 1
    end do

    ty(jxy) = arg
    ny = ny + 1
!  restart the computations with the new set of knots.

  end do

!  test whether the least-squares polynomial is a solution of our
!  approximation problem.
430 if (ier == (-2)) go to 830
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 2: determination of the smoothing spline sp(x,y)                c
! *****************************************************                c
! we have determined the number of knots and their position. we now    c
! compute the b-spline coefficients of the smoothing spline sp(x,y).   c
! the observation matrix a is extended by the rows of a matrix,        c
! expressing that sp(x,y) must be a polynomial of degree kx in x and   c
! ky in y. the corresponding weights of these additional rows are set  c
! to 1./p.  iteratively we than have to determine the value of p       c
! such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c
! we already know that the least-squares polynomial corresponds to     c
! p=0  and that the least-squares spline corresponds to p=infinity.    c
! the iteration process which is proposed here makes use of rational   c
! interpolation. since f(p) is a convex and strictly decreasing        c
! function of p, it can be approximated by a rational function r(p)=   c
! (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
! of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
! new value of p such that r(p)=s. convergence is guaranteed by taking c
! f1 > 0 and f3 < 0.                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  kx2 = kx1 + 1
!  test whether there are interior knots in the x-direction.
  if (nk1x == kx1) go to 440
!  evaluate the discotinuity jumps of the kx-th order derivative of
!  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
  call fpdisc(tx, nx, kx2, bx, nmax)
440 ky2 = ky1 + 1
!  test whether there are interior knots in the y-direction.
  if (nk1y == ky1) go to 450
!  evaluate the discontinuity jumps of the ky-th order derivative of
!  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
  call fpdisc(ty, ny, ky2, by, nmax)
!  initial value for p.
450 p1 = 0.
  f1 = fp0 - s
  p3 = -one
  f3 = fpms
  p = sum(a(:ncof, 1))
  rn = ncof
  p = rn / p
!  find the bandwidth of the extended observation matrix.
  iband3 = kx1 * nk1y
  iband4 = iband3 + 1
  ich1 = 0
  ich3 = 0
!  iteration process to find the root of f(p)=s.
  do iter = 1, maxit
    pinv = one / p
!  store the triangularized observation matrix into q.

    ff(:ncof) = f(:ncof)
    q(:ncof, :iband) = a(:ncof, :iband)
    ibb = iband + 1
    q(:ncof, ibb:iband4) = 0.

    if (nk1y == ky1) go to 560
!  extend the observation matrix with the rows of a matrix, expressing
!  that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
    do i = ky2, nk1y
      ii = i - ky1
      do j = 1, nk1x
!  initialize the new row.
        h(:iband) = 0.

!  fill in the non-zero elements of the row. jrot records the column
!  number of the first non-zero element in the row.

        h(:ky2) = by(ii, :ky2) * pinv
        zi = 0.
        jrot = (j - 1) * nk1y + ii
!  rotate the new row into triangle by givens transformations without
!  square roots.
        do irot = jrot, ncof
          piv = h(1)
          i2 = min0(iband1, ncof - irot)
          if (piv == 0.) then
            if (i2 <= 0) go to 550
            go to 520
          endif
!  calculate the parameters of the givens transformation.
          call fpgivs(piv, q(irot, 1), cos, sin)
!  apply that givens transformation to the right hand side.
          call fprota(cos, sin, zi, ff(irot))
          if (i2 == 0) go to 550
!  apply that givens transformation to the left hand side.
          do l = 1, i2
            l1 = l + 1
            call fprota(cos, sin, h(l1), q(irot, l1))
          end do

520       do l = 1, i2
            h(l) = h(l + 1)
          end do

          h(i2 + 1) = 0.
        end do

550     continue
      end do
    end do

560 if (nk1x == kx1) go to 640
!  extend the observation matrix with the rows of a matrix expressing
!  that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
    do i = kx2, nk1x
      ii = i - kx1
      do j = 1, nk1y
!  initialize the new row
        h(:iband4) = 0.

!  fill in the non-zero elements of the row. jrot records the column
!  number of the first non-zero element in the row.
        j1 = 1
        do l = 1, kx2
          h(j1) = bx(ii, l) * pinv
          j1 = j1 + nk1y
        end do

        zi = 0.
        jrot = (i - kx2) * nk1y + j
!  rotate the new row into triangle by givens transformations .
        do irot = jrot, ncof
          piv = h(1)
          i2 = min0(iband3, ncof - irot)
          if (piv == 0.) then
            if (i2 <= 0) go to 630
            go to 600
          endif
!  calculate the parameters of the givens transformation.
          call fpgivs(piv, q(irot, 1), cos, sin)
!  apply that givens transformation to the right hand side.
          call fprota(cos, sin, zi, ff(irot))
          if (i2 == 0) go to 630
!  apply that givens transformation to the left hand side.
          do l = 1, i2
            l1 = l + 1
            call fprota(cos, sin, h(l1), q(irot, l1))
          end do

600       do l = 1, i2
            h(l) = h(l + 1)
          end do

          h(i2 + 1) = 0.
        end do

630     continue
      end do
    end do

!  find dmax, the maximum value for the diagonal elements in the
!  reduced triangle.
640 dmax = 0.
    do i = 1, ncof
      if (q(i, 1) <= dmax) cycle
      dmax = q(i, 1)
    end do

!  check whether the matrix is rank deficient.
    sigma = eps * dmax
    do i = 1, ncof
      if (q(i, 1) <= sigma) go to 670
    end do

!  backward substitution in case of full rank.
    call fpback(q, ff, ncof, iband4, c, nc)
    rank = ncof
    go to 675
!  in case of rank deficiency, find the minimum norm solution.
670 lwest = ncof * iband4 + ncof + iband4
    if (lwrk < lwest) go to 780
    lf = 1
    lh = lf + ncof
    la = lh + iband4
    call fprank(q, ff, ncof, iband4, nc, sigma, c, sq, rank, wrk(la), wrk(lf), wrk(lh))
675 q(:ncof, 1) = q(:ncof, 1) / dmax

!  compute f(p).
    fp = 0.
    do num = 1, nreg
      num1 = num - 1
      lx = num1 / nyy
      ly = num1 - lx * nyy
      jrot = lx * nk1y + ly
      in = index(num)
690   if (in == 0) cycle
      store = 0.
      i1 = jrot
      do i = 1, kx1
        hxi = spx(in, i)
        j1 = i1
        do j = 1, ky1
          j1 = j1 + 1
          store = store + hxi * spy(in, j) * c(j1)
        end do

        i1 = i1 + nk1y
      end do

      fp = fp + (w(in) * (z(in) - store))**2
      in = nummer(in)
      go to 690
    end do

!  test whether the approximation sp(x,y) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) <= acc) go to 820
!  test whether the maximum allowable number of iterations has been
!  reached.
    if (iter == maxit) go to 795
!  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3 /= 0) go to 740
    if ((f2 - f3) > acc) go to 730
!  our initial choice of p is too large.
    p3 = p2
    f3 = f2
    p = p * con4
    if (p <= p1) p = p1 * con9 + p2 * con1
    go to 770
730 if (f2 < 0.) ich3 = 1
740 if (ich1 /= 0) go to 760
    if ((f1 - f2) > acc) go to 750
!  our initial choice of p is too small
    p1 = p2
    f1 = f2
    p = p / con4
    if (p3 < 0.) go to 770
    if (p >= p3) p = p2 * con1 + p3 * con9
    go to 770
750 if (f2 > 0.) ich1 = 1
!  test whether the iteration process proceeds as theoretically
!  expected.
760 if (f2 >= f1 .or. f2 <= f3) go to 800
!  find the new value of p.
    p = fprati(p1, f1, p2, f2, p3, f3)
770 continue
  end do

!  error codes and messages.
780 ier = lwest
  go to 830
785 ier = 5
  go to 830
790 ier = 4
  go to 830
795 ier = 3
  go to 830
800 ier = 2
  go to 830
810 ier = 1
  go to 830
815 ier = -1
  fp = 0.
820 if (ncof /= rank) ier = -rank
!  test whether x and y are in the original order.
830 if (ichang < 0) go to 930
!  if not, interchange x and y once more.
  l1 = 1
  do i = 1, nk1x
    l2 = i
    do j = 1, nk1y
      f(l2) = c(l1)
      l1 = l1 + 1
      l2 = l2 + nk1x
    end do
  end do

  c(:ncof) = f(:ncof)

  do i = 1, m
    store = x(i)
    x(i) = y(i)
    y(i) = store
  end do

  n = min0(nx, ny)
  do i = 1, n
    store = tx(i)
    tx(i) = ty(i)
    ty(i) = store
  end do

  n1 = n + 1
  if (nx < ny) go to 880
  if (nx == ny) go to 920
  go to 900
880 do 890 i = n1, ny
    tx(i) = ty(i)
890 continue
    go to 920
900 do 910 i = n1, nx
      ty(i) = tx(i)
910   continue
920   l = nx
      nx = ny
      ny = l
930   if (iopt < 0) go to 940
      nx0 = nx
      ny0 = ny
940   return
    end

