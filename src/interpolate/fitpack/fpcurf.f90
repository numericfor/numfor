subroutine fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2,&
  &n, t, c, fp, fpint, z, a, b, g, q, nrdata, ier)
  !  ..
  !  ..scalar arguments..
  real(8) :: xb, xe, s, tol, fp
  integer :: iopt, m, k, nest, maxit, k1, k2, n, ier
  !  ..array arguments..
  real(8) :: x(m), y(m), w(m), t(nest), c(nest), fpint(nest),&
    &z(nest), a(nest, k1), b(nest, k2), g(nest, k2), q(m, k1)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: acc, con1, con4, con9, cos, half, fpart, fpms, fpold, fp0, f1, f2, f3
  real(8) :: one, p, pinv, piv, p1, p2, p3, rn, sin, store, term, wi, xi, yi
  integer :: i, ich1, ich3, it, iter, i1, i2, i3, j, k3, l, l0
  integer :: mk1, new, nk1, nmax, nmin, nplus, npl1, nrint, n8
  !  ..local arrays..
  real(8) :: h(7)
  !  ..function references
  real(8) :: abs, fprati
  ! integer :: max0, min0
  !  ..subroutine references..
  !    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
  !  ..
  !  set constants
  one = 1._8
  con1 = 0.1_8
  con9 = 0.9_8
  con4 = 0.04_8
  half = 0.5_8
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 1: determination of the number of knots and their position     c
  !  **************************************************************      c
  !  given a set of knots we compute the least-squares spline sinf(x),   c
  !  and the corresponding sum of squared residuals fp=f(p=inf).         c
  !  if iopt=-1 sinf(x) is the requested approximation.                  c
  !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
  !    if fp <=s we will continue with the current set of knots.         c
  !    if fp > s we will increase the number of knots and compute the    c
  !       corresponding least-squares spline until finally fp<=s.        c
  !    the initial choice of knots depends on the value of s and iopt.   c
  !    if s=0 we have spline interpolation; in that case the number of   c
  !    knots equals nmax = m+k+1.                                        c
  !    if s > 0 and                                                      c
  !      iopt=0 we first compute the least-squares polynomial of         c
  !      degree k; n = nmin = 2*k+2                                      c
  !      iopt=1 we start with the set of knots found at the last         c
  !      call of the routine, except for the case that s > fp0; then     c
  !      we compute directly the least-squares polynomial of degree k.   c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  determine nmin, the number of knots for polynomial approximation.
  nmin = 2 * k1
  if (iopt < 0) go to 60
  !  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  !  determine nmax, the number of knots for spline interpolation.
  nmax = m + k1
  if (s > 0._8) go to 45
  !  if s=0, s(x) is an interpolating spline.
  !  test whether the required storage space exceeds the available one.
  n = nmax
  if (nmax > nest) go to 420
  !  find the position of the interior knots in case of interpolation.
10 mk1 = m - k1
  if (mk1 == 0) go to 60
  k3 = k / 2
  i = k2
  j = k3 + 2
  if (k3 * 2 == k) go to 30
  do l = 1, mk1
    t(i) = x(j)
    i = i + 1
    j = j + 1
  end do

  go to 60
30 do l = 1, mk1
    t(i) = (x(j) + x(j - 1)) * half
    i = i + 1
    j = j + 1
  end do
  go to 60
  !  if s>0 our initial choice of knots depends on the value of iopt.
  !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
  !  polynomial of degree k which is a spline without interior knots.
  !  if iopt=1 and fp0>s we start computing the least squares spline
  !  according to the set of knots found at the last call of the routine.
45 if (iopt == 0) go to 50
  if (n == nmin) go to 50
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 60
50 n = nmin
  fpold = 0._8
  nplus = 0
  nrdata(1) = m - 2
  !  main loop for the different sets of knots. m is a save upper bound
  !  for the number of trials.
60 do iter = 1, m
    if (n == nmin) ier = -2
    !  find nrint, tne number of knot intervals.
    nrint = n - nmin + 1
    !  find the position of the additional knots which are needed for
    !  the b-spline representation of s(x).
    nk1 = n - k1
    i = n
    do j = 1, k1
      t(j) = xb
      t(i) = xe
      i = i - 1
    end do

    !  compute the b-spline coefficients of the least-squares spline
    !  sinf(x). the observation matrix a is built up row by row and
    !  reduced to upper triangular form by givens transformations.
    !  at the same time fp=f(p=inf) is computed.
    fp = 0._8
    !  initialize the observation matrix a.

    z(:nk1) = 0._8
    a(:nk1, :k1) = 0._8
    l = k1
    do it = 1, m
      !  fetch the current data point x(it),y(it).
      xi = x(it)
      wi = w(it)
      yi = y(it) * wi
      !  search for knot interval t(l) <= xi < t(l+1).
      do while (t(l + 1) <= xi .and. l /= nk1)
        l = l + 1
      end do
      !  evaluate the (k+1) non-zero b-splines at xi and store them in q.
      call fpbspl(t, n, k, xi, l, h)
      q(it, :k1) = h(:k1)
      h(:k1) = h(:k1) * wi
      !  rotate the new row of the observation matrix into triangle.
      j = l - k1
      do i = 1, k1
        j = j + 1
        piv = h(i)
        if (piv == 0._8) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a(j, 1), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, z(j))
        if (i == k1) exit
        i2 = 1
        i3 = i + 1
        do i1 = i3, k1
          i2 = i2 + 1
          !  transformations to left hand side.
          call fprota(cos, sin, h(i1), a(j, i2))
        end do
      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
      fp = fp + yi * yi
    end do

    if (ier == (-2)) fp0 = fp
    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
    !  backward substitution to obtain the b-spline coefficients.
    call fpback(a, z, nk1, k1, c, nest)
    !  test whether the approximation sinf(x) is an acceptable solution.
    if (iopt < 0) go to 440
    fpms = fp - s
    if (abs(fpms) < acc) go to 440
    !  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0._8) go to 250
    !  if n = nmax, sinf(x) is an interpolating spline.
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
    nplus = min(nplus * 2, max(npl1, nplus / 2, 1))
150 fpold = fp
    !  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
    !  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0._8
    i = 1
    l = k2
    new = 0
    do it = 1, m
      if (x(it) >= t(l) .and. l <= nk1) then
        new = 1
        l = l + 1
      end if

      term = 0._8
      l0 = l - k2
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      term = (w(it) * (term - y(it)))**2
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
      call fpknot(x, m, t, n, fpint, nrdata, nrint, nest, 1)
      !  if n=nmax we locate the knots as for interpolation.
      if (n == nmax) go to 10
      !  test whether we cannot further increase the number of knots.
      if (n == nest) exit
    end do

    !  restart the computations with the new set of knots.
    print "(A)", repeat('-', 60)
    print "(5(f0.3,1x))", t
    print "(A)", repeat('-', 60)
  end do

  !  test whether the least-squares kth degree polynomial is a solution
  !  of our approximation problem.
250 if (ier == (-2)) go to 440
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 2: determination of the smoothing spline sp(x).                c
  !  ***************************************************                 c
  !  we have determined the number of knots and their position.          c
  !  we now compute the b-spline coefficients of the smoothing spline    c
  !  sp(x). the observation matrix a is extended by the rows of matrix   c
  !  b expressing that the kth derivative discontinuities of sp(x) at    c
  !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
  !  ponding weights of these additional rows are set to 1/p.            c
  !  iteratively we then have to determine the value of p such that      c
  !  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
  !  the least-squares kth degree polynomial corresponds to p=0, and     c
  !  that the least-squares spline corresponds to p=infinity. the        c
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
  p1 = 0._8
  f1 = fp0 - s
  p3 = -one
  f3 = fpms
  p = sum(a(:nk1, 1))

  rn = nk1
  p = rn / p
  ich1 = 0
  ich3 = 0
  n8 = n - nmin
  !  iteration process to find the root of f(p) = s.
  do iter = 1, maxit
    !  the rows of matrix b with weight 1/p are rotated into the
    !  triangularised observation matrix a which is stored in g.
    pinv = one / p
    c(:nk1) = z(:nk1)
    g(:nk1, k2) = 0._8
    g(:nk1, :k1) = a(:nk1, :k1)
    do it = 1, n8
      !  the row of matrix b is rotated into triangle by givens transformation
      h(:k2) = b(it, :k2) * pinv
      yi = 0._8
      do j = it, nk1
        piv = h(1)
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, g(j, 1), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, c(j))
        if (j == nk1) exit
        i2 = k1
        if (j > n8) i2 = nk1 - j
        do i = 1, i2
          !  transformations to left hand side.
          i1 = i + 1
          call fprota(cos, sin, h(i1), g(j, i1))
          h(i) = h(i1)
        end do

        h(i2 + 1) = 0._8
      end do
    end do

    !  backward substitution to obtain the b-spline coefficients.
    call fpback(g, c, nk1, k2, c, nest)
    !  computation of f(p).
    fp = 0._8
    l = k2
    do it = 1, m
      if (x(it) >= t(l) .and. l <= nk1) l = l + 1
      l0 = l - k2
      term = 0._8
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      fp = fp + (w(it) * (term - y(it)))**2
    end do

    !  test whether the approximation sp(x) is an acceptable solution.
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
335 if (f2 < 0._8) ich3 = 1
340 if (ich1 /= 0) go to 350
    if ((f1 - f2) > acc) go to 345
    !  our initial choice of p is too small
    p1 = p2
    f1 = f2
    p = p / con4
    if (p3 < 0.) cycle
    if (p >= p3) p = p2 * con1 + p3 * con9
    cycle
345 if (f2 > 0._8) ich1 = 1
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
end subroutine fpcurf

