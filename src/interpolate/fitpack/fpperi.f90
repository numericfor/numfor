subroutine fpperi(iopt, x, y, w, m, k, s, nest, tol, maxit, k1, k2, n, t, c,&
  &fp, fpint, z, a1, a2, b, g1, g2, q, nrdata, ier)
  !  ..
  !  ..scalar arguments..
  real(8) :: s, tol, fp
  integer :: iopt, m, k, nest, maxit, k1, k2, n, ier
  !  ..array arguments..
  real(8) :: x(m), y(m), w(m), t(nest), c(nest), fpint(nest), z(nest),&
    &a1(nest, k1), a2(nest, k), b(nest, k2), g1(nest, k2), g2(nest, k1),&
    &q(m, k1)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: acc, cos, c1, d1, fpart, fpms, fpold, fp0, f1, f2, f3, p, per, pinv, piv,&
    &p1, p2, p3, sin, store, term, wi, xi, yi, rn, one, con1, con4, con9, half
  integer :: i, ich1, ich3, ij, ik, it, iter, i1, i2, i3, j, jk, jper, j1, j2, kk,&
    &kk1, k3, l, l0, l1, l5, mm, m1, new, nk1, nk2, nmax, nmin, nplus, npl1,&
    &nrint, n10, n11, n7, n8
  !  ..local arrays..
  real(8) :: h(6), h1(7), h2(6)
  !  ..function references..
  ! real(8) :: abs, fprati
  integer :: max0, min0
  !  ..subroutine references..
  !    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
  !  ..
  !  set constants
  one = 0.1e+01
  con1 = 0.1e0
  con9 = 0.9e0
  con4 = 0.4e-01
  half = 0.5e0
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 1: determination of the number of knots and their position     c
  !  **************************************************************      c
  !  given a set of knots we compute the least-squares periodic spline   c
  !  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c
  !  the initial choice of knots depends on the value of s and iopt.     c
  !    if s=0 we have spline interpolation; in that case the number of   c
  !    knots equals nmax = m+2*k.                                        c
  !    if s > 0 and                                                      c
  !      iopt=0 we first compute the least-squares polynomial of         c
  !      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c
  !      find that s(x) is a constant function.                          c
  !      iopt=1 we start with the set of knots found at the last         c
  !      call of the routine, except for the case that s > fp0; then     c
  !      we compute directly the least-squares periodic polynomial.      c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  m1 = m - 1
  kk = k
  kk1 = k1
  k3 = 3 * k + 1
  nmin = 2 * k1
  !  determine the length of the period of s(x).
  per = x(m) - x(1)
  if (iopt < 0) go to 50
  !  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  !  determine nmax, the number of knots for periodic spline interpolation
  nmax = m + 2 * k
  if (s > 0. .or. nmax == nmin) go to 30
  !  if s=0, s(x) is an interpolating spline.
  n = nmax
  !  test whether the required storage space exceeds the available one.
  if (n > nest) go to 620
  !  find the position of the interior knots in case of interpolation.
5 if ((k / 2) * 2 == k) go to 20
  do i = 2, m1
    j = i + k
    t(j) = x(i)
  end do

  if (s > 0.) go to 50
  kk = k - 1
  kk1 = k
  if (kk > 0) go to 50
  t(1) = t(m) - per
  t(2) = x(1)
  t(m + 1) = x(m)
  t(m + 2) = t(3) + per
  c(:m1) = y(:m1)
  c(m) = c(1)
  fp = 0.
  fpint(n) = fp0
  fpint(n - 1) = 0.
  nrdata(n) = 0
  go to 630
20 do i = 2, m1
    j = i + k
    t(j) = (x(i) + x(i - 1)) * half
  end do

  go to 50
  !  if s > 0 our initial choice depends on the value of iopt.
  !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
  !  periodic polynomial. (i.e. a constant function).
  !  if iopt=1 and fp0>s we start computing the least-squares periodic
  !  spline according the set of knots found at the last call of the
  !  routine.
30 if (iopt == 0) go to 35
  if (n == nmin) go to 35
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 50
  !  the case that s(x) is a constant function is treated separetely.
  !  find the least-squares constant c1 and compute fp0 at the same time.
35 fp0 = 0.
  d1 = 0.
  c1 = 0.
  do it = 1, m1
    wi = w(it)
    yi = y(it) * wi
    call fpgivs(wi, d1, cos, sin)
    call fprota(cos, sin, yi, c1)
    fp0 = fp0 + yi**2
  end do

  c1 = c1 / d1
  !  test whether that constant function is a solution of our problem.
  fpms = fp0 - s
  if (fpms < acc .or. nmax == nmin) go to 640
  fpold = fp0
  !  test whether the required storage space exceeds the available one.
  if (nmin >= nest) go to 620
  !  start computing the least-squares periodic spline with one
  !  interior knot.
  nplus = 1
  n = nmin + 1
  mm = (m + 1) / 2
  t(k2) = x(mm)
  nrdata(1) = mm - 2
  nrdata(2) = m1 - mm
  !  main loop for the different sets of knots. m is a save upper
  !  bound for the number of trials.
50 do iter = 1, m
    !  find nrint, the number of knot intervals.
    nrint = n - nmin + 1
    !  find the position of the additional knots which are needed for
    !  the b-spline representation of s(x). if we take
    !      t(k+1) = x(1), t(n-k) = x(m)
    !      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
    !      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
    !  then s(x) is a periodic spline with period per if the b-spline
    !  coefficients satisfy the following conditions
    !      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
    t(k1) = x(1)
    nk1 = n - k1
    nk2 = nk1 + 1
    t(nk2) = x(m)
    do j = 1, k
      i1 = nk2 + j
      i2 = nk2 - j
      j1 = k1 + j
      j2 = k1 - j
      t(i1) = t(j1) + per
      t(j2) = t(i2) - per
    end do

    !  compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
    !  periodic spline sinf(x). the observation matrix a is built up row
    !  by row while taking into account condition (**) and is reduced to
    !  triangular form by givens transformations .
    !  at the same time fp=f(p=inf) is computed.
    !  the n7 x n7 triangularised upper matrix a has the form
    !            ! a1 '    !
    !        a = !    ' a2 !
    !            ! 0  '    !
    !  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
    !  matrix of bandwidth k+1 ( n10 = n7-k).
    !  initialization.
    z(:nk1) = 0._8
    a1(:nk1, :kk1) = 0._8
    n7 = nk1 - k
    n10 = n7 - kk
    jper = 0
    fp = 0.
    l = k1
    do it = 1, m1
      !  fetch the current data point x(it),y(it)
      xi = x(it)
      wi = w(it)
      yi = y(it) * wi
      !  search for knot interval t(l) <= xi < t(l+1).
      do while (xi >= t(l + 1))
        l = l + 1
      end do

      !  evaluate the (k+1) non-zero b-splines at xi and store them in q.
      call fpbspl(t, n, k, xi, l, h)
      do i = 1, k1
        q(it, i) = h(i)
        h(i) = h(i) * wi
      end do

      l5 = l - k1
      !  test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
      if (l5 < n10) go to 285
      if (jper /= 0) go to 160
      !  initialize the matrix a2.

      a2(:n7, :kk) = 0._8

      jk = n10 + 1
      do i = 1, kk
        ik = jk
        do j = 1, kk1
          if (ik <= 0) exit
          a2(ik, i) = a1(ik, j)
          ik = ik - 1
        end do

        jk = jk + 1
      end do

      jper = 1
      !  if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
      !  we take account of condition (**) for setting up the new row
      !  of the observation matrix a. this row is stored in the arrays h1
      !  (the part with respect to a1) and h2 (the part with
      !  respect to a2).
160   h1(:kk1) = 0._8
      h2(:kk) = 0._8

      j = l5 - n10
      do i = 1, kk1
        j = j + 1
        l0 = j
180     l1 = l0 - kk
        if (l1 <= 0) go to 200
        if (l1 <= n10) go to 190
        l0 = l1 - n10
        go to 180
190     h1(l1) = h(i)
        cycle
200     h2(l0) = h2(l0) + h(i)
      end do

      !  rotate the new row of the observation matrix into triangle
      !  by givens transformations.
      if (n10 <= 0) go to 250
      !  rotation with the rows 1,2,...n10 of matrix a.
      do j = 1, n10
        piv = h1(1)
        if (piv /= 0.) go to 214
        do i = 1, kk
          h1(i) = h1(i + 1)
        end do

        h1(kk1) = 0._8
        cycle
        !  calculate the parameters of the givens transformation.
214     call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformation to the right hand side.
        call fprota(cos, sin, yi, z(j))
        !  transformations to the left hand side with respect to a2.
        do i = 1, kk
          call fprota(cos, sin, h2(i), a2(j, i))
        end do

        if (j == n10) go to 250
        i2 = min0(n10 - j, kk)
        !  transformations to the left hand side with respect to a1.
        do i = 1, i2
          i1 = i + 1
          call fprota(cos, sin, h1(i1), a1(j, i1))
          h1(i) = h1(i1)
        end do

        h1(i1) = 0.
      end do

      !  rotation with the rows n10+1,...n7 of matrix a.
250   do j = 1, kk
        ij = n10 + j
        if (ij <= 0) cycle
        piv = h2(j)
        if (piv == 0.) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a2(ij, j), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, z(ij))
        if (j == kk) go to 280
        j1 = j + 1
        !  transformations to left hand side.
        do i = j1, kk
          call fprota(cos, sin, h2(i), a2(ij, i))
        end do

      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
280   fp = fp + yi**2
      cycle
      !  rotation of the new row of the observation matrix into
      !  triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero
      !  at xi.
285   j = l5
      do i = 1, kk1
        j = j + 1
        piv = h(i)
        if (piv == 0.) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, z(j))
        if (i == kk1) exit
        i2 = 1
        i3 = i + 1
        !  transformations to left hand side.
        do i1 = i3, kk1
          i2 = i2 + 1
          call fprota(cos, sin, h(i1), a1(j, i2))
        end do

      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
      fp = fp + yi**2
    end do

    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
    !  backward substitution to obtain the b-spline coefficients c(j),j=1,.n
    call fpbacp(a1, a2, z, n7, kk, c, kk1, nest)
    !  calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
    do i = 1, k
      j = i + n7
      c(j) = c(i)
    end do

    if (iopt < 0) go to 660
    !  test whether the approximation sinf(x) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) go to 660
    !  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0.) go to 350
    !  if n=nmax, sinf(x) is an interpolating spline.
    if (n == nmax) go to 630
    !  increase the number of knots.
    !  if n=nest we cannot increase the number of knots because of the
    !  storage capacity limitation.
    if (n == nest) go to 620
    !  determine the number of knots nplus we are going to add.
    npl1 = nplus * 2
    rn = nplus
    if (fpold - fp > acc) npl1 = int(rn * fpms / (fpold - fp))
    nplus = min0(nplus * 2, max0(npl1, nplus / 2, 1))
    fpold = fp
    !  compute the sum(wi*(yi-s(xi))**2) for each knot interval
    !  t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0.
    i = 1
    l = k1
    do it = 1, m1
      if (x(it) < t(l)) go to 300
      new = 1
      l = l + 1
300   term = 0.
      l0 = l - k2
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      term = (w(it) * (term - y(it)))**2
      fpart = fpart + term
      if (new == 0) cycle
      if (l <= k2) then
        fpint(nrint) = term
        new = 0
      else
        store = term * half
        fpint(i) = fpart - store
        i = i + 1
        fpart = store
        new = 0
      end if
    end do

    fpint(nrint) = fpint(nrint) + fpart
    do l = 1, nplus
      !  add a new knot
      call fpknot(x, m, t, n, fpint, nrdata, nrint, nest, 1)
      !  if n=nmax we locate the knots as for interpolation.
      if (n == nmax) go to 5
      !  test whether we cannot further increase the number of knots.
      if (n == nest) cycle
    end do

    !  restart the computations with the new set of knots.
  end do

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 2: determination of the smoothing periodic spline sp(x).       c
  !  *************************************************************       c
  !  we have determined the number of knots and their position.          c
  !  we now compute the b-spline coefficients of the smoothing spline    c
  !  sp(x). the observation matrix a is extended by the rows of matrix   c
  !  b expressing that the kth derivative discontinuities of sp(x) at    c
  !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
  !  ponding weights of these additional rows are set to 1/sqrt(p).      c
  !  iteratively we then have to determine the value of p such that      c
  !  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c
  !  the least-squares constant function corresponds to p=0, and that    c
  !  the least-squares periodic spline corresponds to p=infinity. the    c
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
350 call fpdisc(t, n, k2, b, nest)
  !  initial value for p.
  p1 = 0.
  f1 = fp0 - s
  p3 = -one
  f3 = fpms
  n11 = n10 - 1
  n8 = n7 - 1
  p = 0.
  l = n7
  do i = 1, k
    j = k + 1 - i
    p = p + a2(l, j)
    l = l - 1
    if (l == 0) go to 356
  end do

  p = p + sum(a1(:n10, 1))
356 rn = n7
  p = rn / p
  ich1 = 0
  ich3 = 0
  !  iteration process to find the root of f(p) = s.
  do iter = 1, maxit
    !  form the matrix g  as the matrix a extended by the rows of matrix b.
    !  the rows of matrix b with weight 1/p are rotated into
    !  the triangularised observation matrix a.
    !  after triangularisation our n7 x n7 matrix g takes the form
    !            ! g1 '    !
    !        g = !    ' g2 !
    !            ! 0  '    !
    !  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
    !  matrix of bandwidth k+2. ( n11 = n7-k-1)
    pinv = one / p
    !  store matrix a into g

    c(:n7) = z(:n7)
    g1(:n7, k1) = a1(:n7, k1)
    g1(:n7, k2) = 0._8
    g2(:n7, 1) = 0._8
    g1(:n7, :k) = a1(:n7, :k)
    g2(:n7, 2:k + 1) = a2(:n7, :k)
    l = n10
    do j = 1, k1
      if (l <= 0) go to 375
      g2(l, 1) = a1(l, j)
      l = l - 1
    end do

375 do it = 1, n8
      !  fetch a new row of matrix b and store it in the arrays h1 (the part
      !  with respect to g1) and h2 (the part with respect to g2).
      yi = 0.
      h1(:k1) = 0._8
      h1(k2) = 0._8
      h2(:k1) = 0._8

      if (it > n11) go to 420
      l = it
      l0 = it
      do j = 1, k2
        if (l0 == n10) go to 400
        h1(j) = b(it, j) * pinv
        l0 = l0 + 1
      end do

      go to 470
400   l0 = 1
      do l1 = j, k2
        h2(l0) = b(it, l1) * pinv
        l0 = l0 + 1
      end do

      go to 470
420   l = 1
      i = it - n10
      do j = 1, k2
        i = i + 1
        l0 = i
430     l1 = l0 - k1
        if (l1 <= 0) go to 450
        if (l1 <= n11) go to 440
        l0 = l1 - n11
        go to 430
440     h1(l1) = b(it, j) * pinv
        cycle
450     h2(l0) = h2(l0) + b(it, j) * pinv
      end do

      if (n11 <= 0) go to 510
      !  rotate this row into triangle by givens transformations without
      !  square roots.
      !  rotation with the rows l,l+1,...n11.
470   do j = l, n11
        piv = h1(1)
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, g1(j, 1), cos, sin)
        !  transformation to right hand side.
        call fprota(cos, sin, yi, c(j))
        !  transformation to the left hand side with respect to g2.
        do i = 1, k1
          call fprota(cos, sin, h2(i), g2(j, i))
        end do

        if (j == n11) exit
        i2 = min0(n11 - j, k1)
        !  transformation to the left hand side with respect to g1.
        do i = 1, i2
          i1 = i + 1
          call fprota(cos, sin, h1(i1), g1(j, i1))
          h1(i) = h1(i1)
        end do

        h1(i1) = 0._8
      end do

      !  rotation with the rows n11+1,...n7
510   do j = 1, k1
        ij = n11 + j
        if (ij <= 0) cycle
        piv = h2(j)
        !  calculate the parameters of the givens transformation
        call fpgivs(piv, g2(ij, j), cos, sin)
        !  transformation to the right hand side.
        call fprota(cos, sin, yi, c(ij))
        if (j == k1) exit
        j1 = j + 1
        !  transformation to the left hand side.
        do i = j1, k1
          call fprota(cos, sin, h2(i), g2(ij, i))
        end do
      end do
    end do

    !  backward substitution to obtain the b-spline coefficients
    !  c(j),j=1,2,...n7 of sp(x).
    call fpbacp(g1, g2, c, n7, k1, c, k2, nest)
    !  calculate from condition (**) the b-spline coefficients c(n7+j),j=1,.
    do i = 1, k
      j = i + n7
      c(j) = c(i)
    end do

    !  computation of f(p).
    fp = 0.
    l = k1
    do it = 1, m1
      if (x(it) >= t(l)) l = l + 1
      l0 = l - k2
      term = 0.
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      fp = fp + (w(it) * (term - y(it)))**2
    end do

    !  test whether the approximation sp(x) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) go to 660
    !  test whether the maximal number of iterations is reached.
    if (iter == maxit) go to 600
    !  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3 /= 0) go to 580
    if ((f2 - f3) > acc) go to 575
    !  our initial choice of p is too large.
    p3 = p2
    f3 = f2
    p = p * con4
    if (p <= p1) p = p1 * con9 + p2 * con1
    cycle
575 if (f2 < 0.) ich3 = 1
580 if (ich1 /= 0) go to 590
    if ((f1 - f2) > acc) go to 585
    !  our initial choice of p is too small
    p1 = p2
    f1 = f2
    p = p / con4
    if (p3 < 0.) cycle
    if (p >= p3) p = p2 * con1 + p3 * con9
    cycle
585 if (f2 > 0.) ich1 = 1
    !  test whether the iteration process proceeds as theoretically
    !  expected.
590 if (f2 >= f1 .or. f2 <= f3) go to 610
    !  find the new value for p.
    p = fprati(p1, f1, p2, f2, p3, f3)

  end do

  !  error codes and messages.
600 ier = 3
  go to 660
610 ier = 2
  go to 660
620 ier = 1
  go to 660
630 ier = -1
  go to 660
640 ier = -2
  !  the least-squares constant function c1 is a solution of our problem.
  !  a constant function is a spline of degree k with all b-spline
  !  coefficients equal to that constant c1.
  do i = 1, k1
    rn = k1 - i
    t(i) = x(1) - rn * per
    c(i) = c1
    j = i + k1
    rn = i - 1
    t(j) = x(m) + rn * per
  end do

  n = nmin
  fp = fp0
  fpint(n) = fp0
  fpint(n - 1) = 0.
  nrdata(n) = 0
660 return
end subroutine fpperi
