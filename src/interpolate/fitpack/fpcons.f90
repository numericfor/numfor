subroutine fpcons(iopt, idim, m, u, mx, x, w, ib, ie, k, s, nest, tol, maxit,&
        &k1, k2, n, t, nc, c, fp, fpint, z, a, b, g, q, nrdata, ier)
!  ..
!  ..scalar arguments..
  real(8) :: s, tol, fp
  integer :: iopt, idim, m, mx, ib, ie, k, nest, maxit, k1, k2, n, nc, ier
!  ..array arguments..
  real(8) :: u(m), x(mx), w(m), t(nest), c(nc), fpint(nest),&
    &z(nc), a(nest, k1), b(nest, k2), g(nest, k2), q(m, k1)
  integer :: nrdata(nest)
!  ..local scalars..
  real(8) :: acc, con1, con4, con9, cos, fac, fpart, fpms, fpold, fp0, f1, f2, f3,&
    &half, one, p, pinv, piv, p1, p2, p3, rn, sin, store, term, ui, wi
  integer :: i, ich1, ich3, it, iter, i1, i2, i3, j, jb, je, jj, j1, j2, j3, kbe,&
    &l, li, lj, l0, mb, me, mm, new, nk1, nmax, nmin, nn, nplus, npl1, nrint, n8
!  ..local arrays..
  real(8) :: h(7), xi(10)
!  ..function references
  real(8) :: abs, fprati
  integer :: max0, min0
!  ..subroutine references..
!    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
!  ..
!  set constants
  one = 0.1e+01_8
  con1 = 0.1e0_8
  con9 = 0.9e0_8
  con4 = 0.4e-01_8
  half = 0.5e0_8
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 1: determination of the number of knots and their position     c
!  **************************************************************      c
!  given a set of knots we compute the least-squares curve sinf(u),    c
!  and the corresponding sum of squared residuals fp=f(p=inf).         c
!  if iopt=-1 sinf(u) is the requested curve.                          c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!    if fp <=s we will continue with the current set of knots.         c
!    if fp > s we will increase the number of knots and compute the    c
!       corresponding least-squares curve until finally fp<=s.         c
!    the initial choice of knots depends on the value of s and iopt.   c
!    if s=0 we have spline interpolation; in that case the number of   c
!    knots equals nmax = m+k+1-max(0,ib-1)-max(0,ie-1)                 c
!    if s > 0 and                                                      c
!      iopt=0 we first compute the least-squares polynomial curve of   c
!      degree k; n = nmin = 2*k+2                                      c
!      iopt=1 we start with the set of knots found at the last         c
!      call of the routine, except for the case that s > fp0; then     c
!      we compute directly the polynomial curve of degree k.           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  determine nmin, the number of knots for polynomial approximation.
  nmin = 2 * k1
!  find which data points are to be considered.
  mb = 2
  jb = ib
  if (ib > 0) go to 10
  mb = 1
  jb = 1
10 me = m - 1
  je = ie
  if (ie > 0) go to 20
  me = m
  je = 1
20 if (iopt < 0) go to 60
!  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
!  determine nmax, the number of knots for spline interpolation.
  kbe = k1 - jb - je
  mmin = kbe + 2
  mm = m - mmin
  nmax = nmin + mm
  if (s > 0.) go to 40
!  if s=0, s(u) is an interpolating curve.
!  test whether the required storage space exceeds the available one.
  n = nmax
  if (nmax > nest) go to 420
!  find the position of the interior knots in case of interpolation.
  if (mm == 0) go to 60
25 i = k2
  j = 3 - jb + k / 2
  do l = 1, mm
    t(i) = u(j)
    i = i + 1
    j = j + 1
  end do

  go to 60
!  if s>0 our initial choice of knots depends on the value of iopt.
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial curve which is a spline curve without interior knots.
!  if iopt=1 and fp0>s we start computing the least squares spline curve
!  according to the set of knots found at the last call of the routine.
40 if (iopt == 0) go to 50
  if (n == nmin) go to 50
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 60
50 n = nmin
  fpold = 0.
  nplus = 0
  nrdata(1) = m - 2
!  main loop for the different sets of knots. m is a save upper bound
!  for the number of trials.
60 do iter = 1, m
    if (n == nmin) ier = -2
!  find nrint, tne number of knot intervals.
    nrint = n - nmin + 1
!  find the position of the additional knots which are needed for
!  the b-spline representation of s(u).
    nk1 = n - k1
    i = n
    do j = 1, k1
      t(j) = u(1)
      t(i) = u(m)
      i = i - 1
    end do

!  compute the b-spline coefficients of the least-squares spline curve
!  sinf(u). the observation matrix a is built up row by row and
!  reduced to upper triangular form by givens transformations.
!  at the same time fp=f(p=inf) is computed.
    fp = 0.
!  nn denotes the dimension of the splines
    nn = nk1 - ib - ie
!  initialize the b-spline coefficients and the observation matrix a.

    z(:nc) = 0._8
    c(:nc) = 0._8
    if (me < mb) go to 134
    if (nn /= 0) a(:nn, ::k1) = 0._8
    l = k1
    jj = (mb - 1) * idim
    do it = mb, me
!  fetch the current data point u(it),x(it).
      ui = u(it)
      wi = w(it)
      do j = 1, idim
        jj = jj + 1
        xi(j) = x(jj) * wi
      end do

!  search for knot interval t(l) <= ui < t(l+1).
86    if (ui < t(l + 1) .or. l == nk1) go to 90
      l = l + 1
      go to 86
!  evaluate the (k+1) non-zero b-splines at ui and store them in q.
90    call fpbspl(t, n, k, ui, l, h)
      do i = 1, k1
        q(it, i) = h(i)
        h(i) = h(i) * wi
      end do

!  take into account that certain b-spline coefficients must be zero.
      lj = k1
      j = nk1 - l - ie
      if (j >= 0) go to 94
      lj = lj + j
94    li = 1
      j = l - k1 - ib
      if (j >= 0) go to 96
      li = li - j
      j = 0
96    if (li > lj) exit
!  rotate the new row of the observation matrix into triangle.
      do i = li, lj
        j = j + 1
        piv = h(i)
        if (piv == 0.) cycle
!  calculate the parameters of the givens transformation.
        call fpgivs(piv, a(j, 1), cos, sin)
!  transformations to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

        if (i == lj) exit
        i2 = 1
        i3 = i + 1
        do i1 = i3, lj
          i2 = i2 + 1
!  transformations to left hand side.
          call fprota(cos, sin, h(i1), a(j, i2))
        end do

      end do

!  add contribution of this row to the sum of squares of residual
!  right hand sides.
      do j2 = 1, idim
        fp = fp + xi(j2)**2
      end do

    end do

    if (ier == (-2)) fp0 = fp
    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
!  backward substitution to obtain the b-spline coefficients.
    if (nn == 0) go to 134
    j1 = 1
    do j2 = 1, idim
      j3 = j1 + ib
      call fpback(a, z(j1), nn, k1, c(j3), nest)
      j1 = j1 + n
    end do

!  test whether the approximation sinf(u) is an acceptable solution.
134 if (iopt < 0) go to 440
    fpms = fp - s
    if (abs(fpms) < acc) go to 440
!  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0.) go to 250
!  if n = nmax, sinf(u) is an interpolating spline curve.
    if (n == nmax) go to 430
!  increase the number of knots.
!  if n=nest we cannot increase the number of knots because of
!  the storage capacity limitation.
    if (n == nest) go to 420
!  determine the number of knots nplus we are going to add.
    if (ier == 0) go to 140
    nplus = 1
    ier = 0
    go to 150
140 npl1 = nplus * 2
    rn = nplus
    if (fpold - fp > acc) npl1 = int(rn * fpms / (fpold - fp))
    nplus = min0(nplus * 2, max0(npl1, nplus / 2, 1))
150 fpold = fp
!  compute the sum of squared residuals for each knot interval
!  t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0.
    i = 1
    l = k2
    new = 0
    jj = (mb - 1) * idim
    do it = mb, me
      if (u(it) < t(l) .or. l > nk1) go to 160
      new = 1
      l = l + 1
160   term = 0.
      l0 = l - k2
      do j2 = 1, idim
        fac = 0.
        j1 = l0
        do j = 1, k1
          j1 = j1 + 1
          fac = fac + c(j1) * q(it, j)
        end do

        jj = jj + 1
        term = term + (w(it) * (fac - x(jj)))**2
        l0 = l0 + n
      end do

      fpart = fpart + term
      if (new == 0) cycle
      store = term * half
      fpint(i) = fpart - store
      i = i + 1
      fpart = store
      new = 0
    end do

    fpint(nrint) = fpart
    do l = 1, nplus
!  add a new knot.
      call fpknot(u, m, t, n, fpint, nrdata, nrint, nest, 1)
!  if n=nmax we locate the knots as for interpolation
      if (n == nmax) go to 25
!  test whether we cannot further increase the number of knots.
      if (n == nest) exit
    end do

!  restart the computations with the new set of knots.
  end do

!  test whether the least-squares kth degree polynomial curve is a
!  solution of our approximation problem.
250 if (ier == (-2)) go to 440
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 2: determination of the smoothing spline curve sp(u).          c
!  **********************************************************          c
!  we have determined the number of knots and their position.          c
!  we now compute the b-spline coefficients of the smoothing curve     c
!  sp(u). the observation matrix a is extended by the rows of matrix   c
!  b expressing that the kth derivative discontinuities of sp(u) at    c
!  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
!  ponding weights of these additional rows are set to 1/p.            c
!  iteratively we then have to determine the value of p such that f(p),c
!  the sum of squared residuals be = s. we already know that the least c
!  squares kth degree polynomial curve corresponds to p=0, and that    c
!  the least-squares spline curve corresponds to p=infinity. the       c
!  iteration process which is proposed here, makes use of rational     c
!  interpolation. since f(p) is a convex and strictly decreasing       c
!  function of p, it can be approximated by a rational function        c
!  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
!  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
!  to calculate the new value of p such that r(p)=s. convergence is    c
!  guaranteed by taking f1>0 and f3<0.                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  evaluate the discontinuity jump of the kth derivative of the
!  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
  call fpdisc(t, n, k2, b, nest)
!  initial value for p.
  p1 = 0.
  f1 = fp0 - s
  p3 = -one
  f3 = fpms

  p = sum(a(:nn, 1))
  rn = nn
  p = rn / p
  ich1 = 0
  ich3 = 0
  n8 = n - nmin
!  iteration process to find the root of f(p) = s.
  do iter = 1, maxit
!  the rows of matrix b with weight 1/p are rotated into the
!  triangularised observation matrix a which is stored in g.
    pinv = one / p

    c(:nc) = z(:nc)
    g(:nn, k2) = 0._8
    g(:nn, :k1) = a(:nn, :k1)
    do it = 1, n8
!  the row of matrix b is rotated into triangle by given transformation

      h(:k2) = b(it, :k2) * pinv
      xi(:idim) = 0._8
!  take into account that certain b-spline coefficients must be zero.
      if (it > ib) go to 274
      j1 = ib - it + 2
      j2 = 1
      do i = j1, k2
        h(j2) = h(i)
        j2 = j2 + 1
      end do

      h(j2:k2) = 0._8
274   jj = max0(1, it - ib)
      do j = jj, nn
        piv = h(1)
!  calculate the parameters of the givens transformation.
        call fpgivs(piv, g(j, 1), cos, sin)
!  transformations to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), c(j1))
          j1 = j1 + n
        end do

        if (j == nn) exit
        i2 = min0(nn - j, k1)
        do i = 1, i2
!  transformations to left hand side.
          i1 = i + 1
          call fprota(cos, sin, h(i1), g(j, i1))
          h(i) = h(i1)
        end do

        h(i2 + 1) = 0.
      end do

    end do

!  backward substitution to obtain the b-spline coefficients.
    j1 = 1
    do j2 = 1, idim
      j3 = j1 + ib
      call fpback(g, c(j1), nn, k2, c(j3), nest)
      if (ib /= 0) c(j1:j3) = 0._8
      j1 = j1 + n
    end do

!  computation of f(p).
    fp = 0.
    l = k2
    jj = (mb - 1) * idim
    do it = mb, me
      ! if (u(it) < t(l) .or. l > nk1) go to 310
!                   l = l + 1
! 310               l0 = l - k2
      if (u(it) >= t(l) .and. l <= nk1) l = l + 1
      l0 = l - k2
      term = 0.
      do j2 = 1, idim
        fac = sum(c(l0 + 1:l0 + 1 + k1) * q(it, :k1))

        ! fac = 0.
        ! j1 = l0
        ! do j = 1, k1
        !   j1 = j1 + 1
        !   fac = fac + c(j1) * q(it, j)
        ! end do

        jj = jj + 1
        term = term + (fac - x(jj))**2
        l0 = l0 + n
      end do

      fp = fp + term * w(it)**2
    end do

!  test whether the approximation sp(u) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) go to 440
!  test whether the maximal number of iterations is reached.
    if (iter == maxit) go to 400
!  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3 /= 0) go to 340
    if ((f2 - f3) > acc) go to 335
!  our initial choice of p is too large.
    p3 = p2
    f3 = f2
    p = p * con4
    if (p <= p1) p = p1 * con9 + p2 * con1
    cycle
335 if (f2 < 0.) ich3 = 1
340 if (ich1 /= 0) go to 350
    if ((f1 - f2) > acc) go to 345
!  our initial choice of p is too small
    p1 = p2
    f1 = f2
    p = p / con4
    if (p3 < 0.) cycle
    if (p >= p3) p = p2 * con1 + p3 * con9
    cycle
345 if (f2 > 0.) ich1 = 1
!  test whether the iteration process proceeds as theoretically
!  expected.
350 if (f2 >= f1 .or. f2 <= f3) go to 410
!  find the new value for p.
    p = fprati(p1, f1, p2, f2, p3, f3)
  end do

!  error codes and messages.
400 ier = 3
  go to 440
410 ier = 2
  go to 440
420 ier = 1
  go to 440
430 ier = -1
440 return
end
