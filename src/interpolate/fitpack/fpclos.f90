subroutine fpclos(iopt, idim, m, u, mx, x, w, k, s, nest, tol, maxit, k1, k2,&
  &n, t, nc, c, fp, fpint, z, a1, a2, b, g1, g2, q, nrdata, ier)
  !  ..
  !  ..scalar arguments..

  real(8) :: s, tol, fp
  integer :: iopt, idim, m, mx, k, nest, maxit, k1, k2, n, nc, ier
  !  ..array arguments..
  real(8) :: u(m), x(mx), w(m), t(nest), c(nc), fpint(nest), z(nc), a1(nest, k1),&
    & a2(nest, k), b(nest, k2), g1(nest, k2), g2(nest, k1), q(m, k1)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: acc, cos, d1, fac, fpart, fpms, fpold, fp0, f1, f2, f3, p, per, pinv, piv,&
    &p1, p2, p3, sin, store, term, ui, wi, rn, one, con1, con4, con9, half
  integer :: i, ich1, ich3, ij, ik, it, iter, i1, i2, i3, j, jj, jk, jper, j1, j2, kk,&
    &kk1, k3, l, l0, l1, l5, mm, m1, new, nk1, nk2, nmax, nmin, nplus, npl1,&
    &nrint, n10, n11, n7, n8
  !  ..local arrays..
  real(8) :: h(6), h1(7), h2(6), xi(10)
  !  ..function references..
  ! real(8) :: abs, fprati
  ! integer :: max0, min0
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
  !  given a set of knots we compute the least-squares closed curve      c
  !  sinf(u). if the sum f(p=inf) <= s we accept the choice of knots.    c
  !  if iopt=-1 sinf(u) is the requested curve                           c
  !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
  !    if fp <=s we will continue with the current set of knots.         c
  !    if fp > s we will increase the number of knots and compute the    c
  !       corresponding least-squares curve until finally fp<=s.         c
  !  the initial choice of knots depends on the value of s and iopt.     c
  !    if s=0 we have spline interpolation; in that case the number of   c
  !    knots equals nmax = m+2*k.                                        c
  !    if s > 0 and                                                      c
  !      iopt=0 we first compute the least-squares polynomial curve of   c
  !      degree k; n = nmin = 2*k+2. since s(u) must be periodic we      c
  !      find that s(u) reduces to a fixed point.                        c
  !      iopt=1 we start with the set of knots found at the last         c
  !      call of the routine, except for the case that s > fp0; then     c
  !      we compute directly the least-squares polynomial curve.         c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  m1 = m - 1
  kk = k
  kk1 = k1
  k3 = 3 * k + 1
  nmin = 2 * k1
  !  determine the length of the period of the splines.
  per = u(m) - u(1)
  if (iopt < 0) go to 50
  !  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  !  determine nmax, the number of knots for periodic spline interpolation
  nmax = m + 2 * k
  if (s > 0. .or. nmax == nmin) go to 30
  !  if s=0, s(u) is an interpolating curve.
  n = nmax
  !  test whether the required storage space exceeds the available one.
  if (n > nest) then
    ier = 1
    return
  end if

  !  find the position of the interior knots in case of interpolation.
5 if ((k / 2) * 2 == k) go to 20
  do i = 2, m1
    j = i + k
    t(j) = u(i)
  end do

  if (s > 0.) go to 50
  kk = k - 1
  kk1 = k
  if (kk > 0) go to 50
  t(1) = t(m) - per
  t(2) = u(1)
  t(m + 1) = u(m)
  t(m + 2) = t(3) + per
  jj = 0
  do i = 1, m1
    j = i
    do j1 = 1, idim
      jj = jj + 1
      c(j) = x(jj)
      j = j + n
    end do
  end do

  jj = 1
  j = m
  do j1 = 1, idim
    c(j) = c(jj)
    j = j + n
    jj = jj + n
  end do

  fp = 0.
  fpint(n) = fp0
  fpint(n - 1) = 0.
  nrdata(n) = 0
  ier = -1
  return

20 do i = 2, m1
    j = i + k
    t(j) = (u(i) + u(i - 1)) * half
  end do

  go to 50
  !  if s > 0 our initial choice depends on the value of iopt.
  !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
  !  polynomial curve. (i.e. a constant point).
  !  if iopt=1 and fp0>s we start computing the least-squares closed
  !  curve according the set of knots found at the last call of the
  !  routine.
30 if (iopt == 0) go to 35
  if (n == nmin) go to 35
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 50
  !  the case that s(u) is a fixed point is treated separetely.
  !  fp0 denotes the corresponding sum of squared residuals.
35 fp0 = 0.
  d1 = 0.
  z(:idim) = 0._8
  jj = 0
  do it = 1, m1
    wi = w(it)
    call fpgivs(wi, d1, cos, sin)
    do j = 1, idim
      jj = jj + 1
      fac = wi * x(jj)
      call fprota(cos, sin, fac, z(j))
      fp0 = fp0 + fac**2
    end do
  end do

  z(:idim) = z(:idim) / d1

  !  test whether that fixed point is a solution of our problem.
  fpms = fp0 - s
  if (fpms < acc .or. nmax == nmin) go to 640
  fpold = fp0
  !  test whether the required storage space exceeds the available one.
  if (n >= nest) then
    ier = 1
    return
  end if

  !  start computing the least-squares closed curve with one
  !  interior knot.
  nplus = 1
  n = nmin + 1
  mm = (m + 1) / 2
  t(k2) = u(mm)
  nrdata(1) = mm - 2
  nrdata(2) = m1 - mm
  !  main loop for the different sets of knots. m is a save upper
  !  bound for the number of trials.
50 do iter = 1, m
    !  find nrint, the number of knot intervals.
    nrint = n - nmin + 1
    !  find the position of the additional knots which are needed for
    !  the b-spline representation of s(u). if we take
    !      t(k+1) = u(1), t(n-k) = u(m)
    !      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
    !      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
    !  then s(u) will be a smooth closed curve if the b-spline
    !  coefficients satisfy the following conditions
    !      c((i-1)*n+n7+j) = c((i-1)*n+j), j=1,...k,i=1,2,...,idim (**)
    !  with n7=n-2*k-1.
    t(k1) = u(1)
    nk1 = n - k1
    nk2 = nk1 + 1
    t(nk2) = u(m)
    do j = 1, k
      i1 = nk2 + j
      i2 = nk2 - j
      j1 = k1 + j
      j2 = k1 - j
      t(i1) = t(j1) + per
      t(j2) = t(i2) - per
    end do

    !  compute the b-spline coefficients of the least-squares closed curve
    !  sinf(u). the observation matrix a is built up row by row while
    !  taking into account condition (**) and is reduced to triangular
    !  form by givens transformations .
    !  at the same time fp=f(p=inf) is computed.
    !  the n7 x n7 triangularised upper matrix a has the form
    !            ! a1 '    !
    !        a = !    ' a2 !
    !            ! 0  '    !
    !  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
    !  matrix of bandwidth k+1 ( n10 = n7-k).
    !  initialization.

    z(:nc) = 0._8
    a1(:nk1, :kk1) = 0._8

    n7 = nk1 - k
    n10 = n7 - kk
    jper = 0
    fp = 0.
    l = k1
    jj = 0
    item: do it = 1, m1
      !  fetch the current data point u(it),x(it)
      ui = u(it)
      wi = w(it)
      do j = 1, idim
        jj = jj + 1
        xi(j) = x(jj) * wi
      end do

      !  search for knot interval t(l) <= ui < t(l+1).
80    if (ui < t(l + 1)) go to 85
      l = l + 1
      go to 80
      !  evaluate the (k+1) non-zero b-splines at ui and store them in q.
85    call fpbspl(t, n, k, ui, l, h)
      do i = 1, k1
        q(it, i) = h(i)
        h(i) = h(i) * wi
      end do

      l5 = l - k1
      !  test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui
      if (l5 < n10) go to 285
      if (jper /= 0) go to 160
      !  initialize the matrix a2.
      a2(:n7, :kk) = 0._8
      jk = n10 + 1
      do i = 1, kk
        ik = jk
        do j = 1, kk1
          if (ik <= 0) cycle
          a2(ik, i) = a1(ik, j)
          ik = ik - 1
        end do

        jk = jk + 1
      end do

      jper = 1
      !  if one of the b-splines nj,k+1(u),j=n7+1,...nk1 is not zero at ui
      !  we take account of condition (**) for setting up the new row
      !  of the observation matrix a. this row is stored in the arrays h1
      !  (the part with respect to a1) and h2 (the part with
      !  respect to a2).
160   h1(:kk) = 0._8
      h2(:kk) = 0._8
      h1(kk1) = 0.
      j = l5 - n10
      do i = 1, kk1
        j = j + 1
        l0 = j
        do
180       l1 = l0 - kk
          if (l1 <= 0) go to 200
          if (l1 <= n10) go to 190
          l0 = l1 - n10
          go to 180
        end do

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
        h1(kk1) = 0.
        cycle ! go to 240
        !  calculate the parameters of the givens transformation.
214     call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformation to the right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

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
        j1 = ij
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

        if (j == kk) go to 280
        j1 = j + 1
        !  transformations to left hand side.
        do i = j1, kk
          call fprota(cos, sin, h2(i), a2(ij, i))
        end do

      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
280   do j2 = 1, idim
        fp = fp + xi(j2)**2
      end do

      cycle item
      !  rotation of the new row of the observation matrix into
      !  triangle in case the b-splines nj,k+1(u),j=n7+1,...n-k-1 are all zero
      !  at ui.
285   j = l5
      do i = 1, kk1
        j = j + 1
        piv = h(i)
        if (piv == 0.) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformations to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

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
      do j2 = 1, idim
        fp = fp + xi(j2)**2
      end do

    end do item

    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
    !  backward substitution to obtain the b-spline coefficients .
    j1 = 1
    do j2 = 1, idim
      call fpbacp(a1, a2, z(j1), n7, kk, c(j1), kk1, nest)
      j1 = j1 + n
    end do

    !  calculate from condition (**) the remaining coefficients.
    do i = 1, k
      j1 = i
      do j = 1, idim
        j2 = j1 + n7
        c(j2) = c(j1)
        j1 = j1 + n
      end do
    end do

    if (iopt < 0) return
    !  test whether the approximation sinf(u) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) return
    !  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0.) go to 350
    !  if n=nmax, sinf(u) is an interpolating curve.
    if (n == nmax) then
      ier = -1
      return
    end if

    !  increase the number of knots.
    !  if n=nest we cannot increase the number of knots because of the
    !  storage capacity limitation.
    if (n == nest) then
      ier = 1
      return
    end if

    !  determine the number of knots nplus we are going to add.
    npl1 = nplus * 2
    rn = nplus
    if (fpold - fp > acc) npl1 = int(rn * fpms / (fpold - fp))
    nplus = min(nplus * 2, max(npl1, nplus / 2, 1))
    fpold = fp
    !  compute the sum of squared residuals for each knot interval
    !  t(j+k) <= ui <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0.
    i = 1
    l = k1
    jj = 0
    do it = 1, m1
      if (u(it) < t(l)) go to 300
      new = 1
      l = l + 1
300   term = 0.
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
      if (l > k2) go to 315
      fpint(nrint) = term
      new = 0
      cycle
315   store = term * half
      fpint(i) = fpart - store
      i = i + 1
      fpart = store
      new = 0
    end do

    fpint(nrint) = fpint(nrint) + fpart
    knots: do l = 1, nplus
      !  add a new knot
      call fpknot(u, m, t, n, fpint, nrdata, nrint, nest, 1)
      !  if n=nmax we locate the knots as for interpolation
      if (n == nmax) go to 5
      !  test whether we cannot further increase the number of knots.
      if (n == nest) exit knots
    end do knots

    !  restart the computations with the new set of knots.
  end do

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 2: determination of the smoothing closed curve sp(u).          c
  !  **********************************************************          c
  !  we have determined the number of knots and their position.          c
  !  we now compute the b-spline coefficients of the smoothing curve     c
  !  sp(u). the observation matrix a is extended by the rows of matrix   c
  !  b expressing that the kth derivative discontinuities of sp(u) at    c
  !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
  !  ponding weights of these additional rows are set to 1/p.            c
  !  iteratively we then have to determine the value of p such that f(p),c
  !  the sum of squared residuals be = s. we already know that the least-c
  !  squares polynomial curve corresponds to p=0, and that the least-    c
  !  squares periodic spline curve corresponds to p=infinity. the        c
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
  ! do i = 1, n10
  !   p = p + a1(i, 1)
  ! end do

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
    c(:nc) = z(:nc)
    do i = 1, n7
      g1(i, k1) = a1(i, k1)
      g1(i, k2) = 0.
      g2(i, 1) = 0.
      do j = 1, k
        g1(i, j) = a1(i, j)
        g2(i, j + 1) = a2(i, j)
      end do
    end do

    l = n10
    do j = 1, k1
      if (l <= 0) exit
      g2(l, 1) = a1(l, j)
      l = l - 1
    end do

    do it = 1, n8
      !  fetch a new row of matrix b and store it in the arrays h1 (the part
      !  with respect to g1) and h2 (the part with respect to g2).
      xi(:idim) = 0._8
      h1(:k1) = 0._8
      h2(:k1) = 0._8
      h1(k2) = 0.
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
      !  rotate this row into triangle by givens transformations
      !  rotation with the rows l,l+1,...n11.
470   do j = l, n11
        piv = h1(1)
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, g1(j, 1), cos, sin)
        !  transformation to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), c(j1))
          j1 = j1 + n
        end do

        !  transformation to the left hand side with respect to g2.
        do i = 1, k1
          call fprota(cos, sin, h2(i), g2(j, i))
        end do

        if (j == n11) go to 510
        i2 = min0(n11 - j, k1)
        !  transformation to the left hand side with respect to g1.
        do i = 1, i2
          i1 = i + 1
          call fprota(cos, sin, h1(i1), g1(j, i1))
          h1(i) = h1(i1)
        end do

        h1(i1) = 0.
      end do

      !  rotation with the rows n11+1,...n7
510   do j = 1, k1
        ij = n11 + j
        if (ij <= 0) cycle
        piv = h2(j)
        !  calculate the parameters of the givens transformation
        call fpgivs(piv, g2(ij, j), cos, sin)
        !  transformation to the right hand side.
        j1 = ij
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), c(j1))
          j1 = j1 + n
        end do

        if (j == k1) exit
        j1 = j + 1
        !  transformation to the left hand side.
        do i = j1, k1
          call fprota(cos, sin, h2(i), g2(ij, i))
        end do

      end do

    end do

    !  backward substitution to obtain the b-spline coefficients
    j1 = 1
    do j2 = 1, idim
      call fpbacp(g1, g2, c(j1), n7, k1, c(j1), k2, nest)
      j1 = j1 + n
    end do

    !  calculate from condition (**) the remaining b-spline coefficients.
    do i = 1, k
      j1 = i
      do j = 1, idim
        j2 = j1 + n7
        c(j2) = c(j1)
        j1 = j1 + n
      end do
    end do

    !  computation of f(p).
    fp = 0.
    l = k1
    jj = 0
    do it = 1, m1
      if (u(it) < t(l)) go to 550
      l = l + 1
550   l0 = l - k2
      term = 0.
      do j2 = 1, idim
        fac = 0.
        j1 = l0
        do j = 1, k1
          j1 = j1 + 1
          fac = fac + c(j1) * q(it, j)
        end do

        jj = jj + 1
        term = term + (fac - x(jj))**2
        l0 = l0 + n
      end do

      fp = fp + term * w(it)**2
    end do

    !  test whether the approximation sp(u) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) return
    !  test whether the maximal number of iterations is reached.
    if (iter == maxit) then
      ier = 3
      return
    end if

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
    if (p3 < 0._8) cycle
    if (p >= p3) p = p2 * con1 + p3 * con9
    cycle
585 if (f2 > 0._8) ich1 = 1
    !  test whether the iteration process proceeds as theoretically
    !  expected.
590 if (f2 >= f1 .or. f2 <= f3) then
      ier = 2
      return
    end if

    !  find the new value for p.
    p = fprati(p1, f1, p2, f2, p3, f3)
  end do

  !  error codes and messages.
640 ier = -2
  !  the point (z(1),z(2),...,z(idim)) is a solution of our problem.
  !  a constant function is a spline of degree k with all b-spline
  !  coefficients equal to that constant.
  do i = 1, k1
    rn = k1 - i
    t(i) = u(1) - rn * per
    j = i + k1
    rn = i - 1
    t(j) = u(m) + rn * per
  end do

  n = nmin
  j1 = 0
  do j = 1, idim
    fac = z(j)
    j2 = j1
    do i = 1, k1
      j2 = j2 + 1
      c(j2) = fac
    end do
    j1 = j1 + n
  end do
  fp = fp0
  fpint(n) = fp0
  fpint(n - 1) = 0._8
  nrdata(n) = 0
end subroutine fpclos
