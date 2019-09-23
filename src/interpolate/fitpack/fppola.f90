subroutine fppola(iopt1, iopt2, iopt3, m, u, v, z, w, rad, s, nuest, nvest,&
  &eta, tol, maxit, ib1, ib3, nc, ncc, intest, nrest, nu, tu, nv, tv, c, fp, sup,&
  &fpint, coord, f, ff, row, cs, cosi, a, q, bu, bv, spu, spv, h, index, nummer,&
  &wrk, lwrk, ier)
!  ..scalar arguments..
  integer :: iopt1, iopt2, iopt3, m, nuest, nvest, maxit, ib1, ib3, nc, ncc,&
    &intest, nrest, nu, nv, lwrk, ier
  real(8) :: s, eta, tol, fp, sup
!  ..array arguments..
  integer :: index(nrest), nummer(m)
  real(8) :: u(m), v(m), z(m), w(m), tu(nuest), tv(nvest), c(nc), fpint(intest),&
    &coord(intest), f(ncc), ff(nc), row(nvest), cs(nvest), cosi(5, nvest),&
    &a(ncc, ib1), q(ncc, ib3), bu(nuest, 5), bv(nvest, 5), spu(m, 4), spv(m, 4),&
    &h(ib3), wrk(lwrk)
!  ..user supplied function..
  real(8) :: rad
!  ..local scalars..
  real(8) :: acc, arg, co, c1, c2, c3, c4, dmax, eps, fac, fac1, fac2, fpmax, fpms,&
    &f1, f2, f3, hui, huj, p, pi, pinv, piv, pi2, p1, p2, p3, r, ratio, si, sigma,&
    &sq, store, uu, u2, u3, wi, zi, rn, one, two, three, con1, con4, con9, half, ten
  integer :: i, iband, iband3, iband4, ich1, ich3, ii, il, in, ipar, ipar1, irot,&
    &iter, i1, i2, i3, j, jrot, j1, j2, k, l, la, lf, lh, ll, lu, lv, lwest, l1, l2,&
    &l3, l4, ncof, ncoff, nvv, nv4, nreg, nrint, nrr, nr1, nuu, nu4, num, num1,&
    &numin, nvmin, rank, iband1
!  ..local arrays..
  real(8) :: hu(4), hv(4)
!  ..function references..
  real(8) :: abs, atan, cos, fprati, sin, sqrt
  integer :: min0
!  ..subroutine references..
!    fporde,fpbspl,fpback,fpgivs,fprota,fprank,fpdisc,fprppo
!  ..
!  set constants
  one = 1
  two = 2
  three = 3
  ten = 10
  half = 0.5e0
  con1 = 0.1e0
  con9 = 0.9e0
  con4 = 0.4e-01
  pi = atan(one) * 4
  pi2 = pi + pi
  ipar = iopt2 * (iopt2 + 3) / 2
  ipar1 = ipar + 1
  eps = sqrt(eta)
  if (iopt1 < 0) go to 90
  numin = 9
  nvmin = 9 + iopt2 * (iopt2 + 1)
!  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  if (iopt1 == 0) go to 10
  if (s < sup) then
    if (nv < nvmin) go to 70
    go to 90
  endif
!  if iopt1 = 0 we begin by computing the weighted least-squares
!  polymomial of the form
!     s(u,v) = f(1)*(1-u**3)+f(2)*u**3+f(3)*(u**2-u**3)+f(4)*(u-u**3)
!  where f(4) = 0 if iopt2> 0 , f(3) = 0 if iopt2 > 1 and
!        f(2) = 0 if iopt3> 0.
!  the corresponding weighted sum of squared residuals gives the upper
!  bound sup for the smoothing factor s.
10 sup = 0.
  f(:4) = 0._8
  a(:4, :4) = 0._8
  do i = 1, m
    wi = w(i)
    zi = z(i) * wi
    uu = u(i)
    u2 = uu * uu
    u3 = uu * u2
    h(1) = (one - u3) * wi
    h(2) = u3 * wi
    h(3) = u2 * (one - uu) * wi
    h(4) = uu * (one - u2) * wi
    if (iopt3 /= 0) h(2) = 0.
    if (iopt2 > 1) h(3) = 0.
    if (iopt2 > 0) h(4) = 0.
    do j = 1, 4
      piv = h(j)
      if (piv == 0.) cycle
      call fpgivs(piv, a(j, 1), co, si)
      call fprota(co, si, zi, f(j))
      if (j == 4) cycle
      j1 = j + 1
      j2 = 1
      do l = j1, 4
        j2 = j2 + 1
        call fprota(co, si, h(l), a(j, j2))
      end do

    end do

    sup = sup + zi * zi
  end do

  if (a(4, 1) /= 0.) f(4) = f(4) / a(4, 1)
  if (a(3, 1) /= 0.) f(3) = (f(3) - a(3, 2) * f(4)) / a(3, 1)
  if (a(2, 1) /= 0.) f(2) = (f(2) - a(2, 2) * f(3) - a(2, 3) * f(4)) / a(2, 1)
  if (a(1, 1) /= 0.) then
    f(1) = (f(1) - a(1, 2) * f(2) - a(1, 3) * f(3) - a(1, 4) * f(4)) / a(1, 1)
  end if

!  find the b-spline representation of this least-squares polynomial
  c1 = f(1)
  c4 = f(2)
  c2 = f(4) / three + c1
  c3 = (f(3) + two * f(4)) / three + c1
  nu = 8
  nv = 8
  do i = 1, 4
    c(i) = c1
    c(i + 4) = c2
    c(i + 8) = c3
    c(i + 12) = c4
    tu(i) = 0.
    tu(i + 4) = one
    rn = 2 * i - 9
    tv(i) = rn * pi
    rn = 2 * i - 1
    tv(i + 4) = rn * pi
  end do

  fp = sup
!  test whether the least-squares polynomial is an acceptable solution
  fpms = sup - s
  if (fpms < acc) go to 960
!  test whether we cannot further increase the number of knots.
70 if (nuest < numin .or. nvest < nvmin) go to 950
!  find the initial set of interior knots of the spline in case iopt1=0.
  nu = numin
  nv = nvmin
  tu(5) = half
  nvv = nv - 8
  rn = nvv + 1
  fac = pi2 / rn
  do i = 1, nvv
    rn = i
    tv(i + 4) = rn * fac - pi
  end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 1 : computation of least-squares bicubic splines.              c
!  ******************************************************              c
!  if iopt1<0 we compute the least-squares bicubic spline according    c
!  to the given set of knots.                                          c
!  if iopt1>=0 we compute least-squares bicubic splines with in-       c
!  creasing numbers of knots until the corresponding sum f(p=inf)<=s.  c
!  the initial set of knots then depends on the value of iopt1         c
!    if iopt1=0 we start with one interior knot in the u-direction     c
!              (0.5) and 1+iopt2*(iopt2+1) in the v-direction.         c
!    if iopt1>0 we start with the set of knots found at the last       c
!              call of the routine.                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  main loop for the different sets of knots. m is a save upper bound
!  for the number of trials.
90 do iter = 1, m
!  find the position of the additional knots which are needed for the
!  b-spline representation of s(u,v).
    l1 = 4
    l2 = l1
    l3 = nv - 3
    l4 = l3
    tv(l2) = -pi
    tv(l3) = pi
    do i = 1, 3
      l1 = l1 + 1
      l2 = l2 - 1
      l3 = l3 + 1
      l4 = l4 - 1
      tv(l2) = tv(l4) - pi2
      tv(l3) = tv(l1) + pi2
    end do

    l = nu
    do i = 1, 4
      tu(i) = 0.
      tu(l) = one
      l = l - 1
    end do

!  find nrint, the total number of knot intervals and nreg, the number
!  of panels in which the approximation domain is subdivided by the
!  intersection of knots.
    nuu = nu - 7
    nvv = nv - 7
    nrr = nvv / 2
    nr1 = nrr + 1
    nrint = nuu + nvv
    nreg = nuu * nvv
!  arrange the data points according to the panel they belong to.
    call fporde(u, v, m, 3, 3, tu, nu, tv, nv, nummer, index, nreg)
    if (iopt2 == 0) go to 195
!  find the b-spline coefficients cosi of the cubic spline
!  approximations for cr(v)=rad(v)*cos(v) and sr(v) = rad(v)*sin(v)
!  if iopt2=1, and additionally also for cr(v)**2,sr(v)**2 and
    !  2*cr(v)*sr(v) if iopt2=2
    cosi(:nvv, :ipar) = 0._8
    a(:nvv, :nvv) = 0._8
!  the coefficients cosi are obtained from interpolation conditions
!  at the knots tv(i),i=4,5,...nv-4.
    do i = 1, nvv
      l2 = i + 3
      arg = tv(l2)
      call fpbspl(tv, nv, 3, arg, l2, hv)

      row(:nvv) = 0.

      ll = i
      do j = 1, 3
        if (ll > nvv) ll = 1
        row(ll) = row(ll) + hv(j)
        ll = ll + 1
      end do

      co = cos(arg)
      si = sin(arg)
      r = rad(arg)
      cs(1) = co * r
      cs(2) = si * r
      if (iopt2 == 1) go to 155
      cs(3) = cs(1) * cs(1)
      cs(4) = cs(2) * cs(2)
      cs(5) = cs(1) * cs(2)
155   do j = 1, nvv
        piv = row(j)
        if (piv == 0.) cycle
        call fpgivs(piv, a(j, 1), co, si)
        do l = 1, ipar
          call fprota(co, si, cs(l), cosi(l, j))
        end do

        if (j == nvv) exit
        j1 = j + 1
        j2 = 1
        do l = j1, nvv
          j2 = j2 + 1
          call fprota(co, si, row(l), a(j, j2))
        end do

      end do

    end do

    do l = 1, ipar

      cs(:nvv) = cosi(l, :nvv)

      call fpback(a, cs, nvv, nvv, cs, ncc)

      cosi(l, :nvv) = cs(:nvv)

    end do

!  find ncof, the dimension of the spline and ncoff, the number
!  of coefficients in the standard b-spline representation.
195 nu4 = nu - 4
    nv4 = nv - 4
    ncoff = nu4 * nv4
    ncof = ipar1 + nvv * (nu4 - 1 - iopt2 - iopt3)
!  find the bandwidth of the observation matrix a.
    iband = 4 * nvv
    if (nuu - iopt2 - iopt3 <= 1) iband = ncof
    iband1 = iband - 1
    !  initialize the observation matrix a.
    f(:ncof) = 0._8
    a(:ncof, :iband) = 0._8
!  initialize the sum of squared residuals.
    fp = 0.
    ratio = one + tu(6) / tu(5)
!  fetch the data points in the new order. main loop for the
!  different panels.
    do 380 num = 1, nreg
!  fix certain constants for the current panel; jrot records the column
!  number of the first non-zero element in a row of the observation
!  matrix according to a data point of the panel.
      num1 = num - 1
      lu = num1 / nvv
      l1 = lu + 4
      lv = num1 - lu * nvv + 1
      l2 = lv + 3
      jrot = 0
      if (lu > iopt2) jrot = ipar1 + (lu - iopt2 - 1) * nvv
      lu = lu + 1
!  test whether there are still data points in the current panel.
      in = index(num)
210   if (in == 0) go to 380
!  fetch a new data point.
      wi = w(in)
      zi = z(in) * wi
!  evaluate for the u-direction, the 4 non-zero b-splines at u(in)
      call fpbspl(tu, nu, 3, u(in), l1, hu)
!  evaluate for the v-direction, the 4 non-zero b-splines at v(in)
      call fpbspl(tv, nv, 3, v(in), l2, hv)
!  store the value of these b-splines in spu and spv resp.
      spu(in, :4) = hu(:4)
      spv(in, :4) = hv(:4)

!  initialize the new row of observation matrix.
      h(:iband) = 0._8

!  calculate the non-zero elements of the new row by making the cross
!  products of the non-zero b-splines in u- and v-direction and
!  by taking into account the conditions of the splines.
      row(:nvv) = 0._8
!  take into account the periodicity condition of the bicubic splines.
      ll = lv
      do i = 1, 4
        if (ll > nvv) ll = 1
        row(ll) = row(ll) + hv(i)
        ll = ll + 1
      end do

!  take into account the other conditions of the splines.
      if (iopt2 == 0 .or. lu > iopt2 + 1) go to 280
      do l = 1, ipar
        cs(l) = 0.
        do i = 1, nvv
          cs(l) = cs(l) + row(i) * cosi(l, i)
        end do
      end do

!  fill in the non-zero elements of the new row.
280   j1 = 0
      do 330 j = 1, 4
        jlu = j + lu
        huj = hu(j)
        if (jlu > iopt2 + 2) go to 320
        if ((jlu == 1) .or. (jlu == 2)) then
          go to 290
        elseif (jlu == 3) then
          go to 300
        elseif (jlu == 4) then
          go to 310
        end if
290     h(1) = huj
        j1 = 1
        go to 330
300     h(1) = h(1) + huj
        h(2) = huj * cs(1)
        h(3) = huj * cs(2)
        j1 = 3
        go to 330
310     h(1) = h(1) + huj
        h(2) = h(2) + huj * ratio * cs(1)
        h(3) = h(3) + huj * ratio * cs(2)
        h(4) = huj * cs(3)
        h(5) = huj * cs(4)
        h(6) = huj * cs(5)
        j1 = 6
        go to 330
320     if (jlu > nu4 .and. iopt3 /= 0) go to 330
        do i = 1, nvv
          j1 = j1 + 1
          h(j1) = row(i) * huj
        end do
330   end do

      h(:iband) = h(:iband) * wi

!  rotate the row into triangle by givens transformations.
      irot = jrot
      do 350 i = 1, iband
        irot = irot + 1
        piv = h(i)
        if (piv == 0.) go to 350
!  calculate the parameters of the givens transformation.
        call fpgivs(piv, a(irot, 1), co, si)
!  apply that transformation to the right hand side.
        call fprota(co, si, zi, f(irot))
        if (i == iband) go to 360
!  apply that transformation to the left hand side.
        i2 = 1
        i3 = i + 1
        do j = i3, iband
          i2 = i2 + 1
          call fprota(co, si, h(j), a(irot, i2))
        end do

350   end do

!  add the contribution of the row to the sum of squares of residual
!  right hand sides.
360   fp = fp + zi**2
!  find the number of the next data point in the panel.
      in = nummer(in)
      go to 210
380   continue
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
        if (a(i, 1) <= sigma) go to 410
      end do

!  backward substitution in case of full rank.
      call fpback(a, f, ncof, iband, c, ncc)
      rank = ncof
      q(:ncof, 1) = a(:ncof, 1) / dmax

      go to 430
!  in case of rank deficiency, find the minimum norm solution.
410   lwest = ncof * iband + ncof + iband
      if (lwrk < lwest) go to 925
      lf = 1
      lh = lf + ncof
      la = lh + iband
      ff(:ncof) = f(:ncof)
      q(:ncof, :iband) = a(:ncof, :iband)

      call fprank(q, ff, ncof, iband, ncc, sigma, c, sq, rank, wrk(la),&
        &wrk(lf), wrk(lh))

      q(:ncof, 1) = q(:ncof, 1) / dmax

!  add to the sum of squared residuals, the contribution of reducing
!  the rank.
      fp = fp + sq
!  find the coefficients in the standard b-spline representation of
!  the spline.
430   call fprppo(nu, nv, iopt2, iopt3, cosi, ratio, c, ff, ncoff)
!  test whether the least-squares spline is an acceptable solution.
      if (iopt1 < 0) then
        if (fp <= 0) go to 970
        go to 980
      endif
      fpms = fp - s
      if (abs(fpms) <= acc) then
        if (fp <= 0) go to 970
        go to 980
      endif
!  if f(p=inf) < s, accept the choice of knots.
      if (fpms < 0.) go to 580
!  test whether we cannot further increase the number of knots
      if (m < ncof) go to 935
!  search where to add a new knot.
!  find for each interval the sum of squared residuals fpint for the
!  data points having the coordinate belonging to that knot interval.
!  calculate also coord which is the same sum, weighted by the position
!  of the data points considered.
      do i = 1, nrint
        fpint(i) = 0.
        coord(i) = 0.
      end do

      do 490 num = 1, nreg
        num1 = num - 1
        lu = num1 / nvv
        l1 = lu + 1
        lv = num1 - lu * nvv
        l2 = lv + 1 + nuu
        jrot = lu * nv4 + lv
        in = index(num)
460     if (in == 0) go to 490
        store = 0.
        i1 = jrot
        do i = 1, 4
          hui = spu(in, i)
          j1 = i1
          do j = 1, 4
            j1 = j1 + 1
            store = store + hui * spv(in, j) * c(j1)
          end do

          i1 = i1 + nv4
        end do

        store = (w(in) * (z(in) - store))**2
        fpint(l1) = fpint(l1) + store
        coord(l1) = coord(l1) + store * u(in)
        fpint(l2) = fpint(l2) + store
        coord(l2) = coord(l2) + store * v(in)
        in = nummer(in)
        go to 460
490   end do

! bring together the information concerning knot panels which are
! symmetric with respect to the origin.
      do i = 1, nrr
        l1 = nuu + i
        l2 = l1 + nrr
        fpint(l1) = fpint(l1) + fpint(l2)
        coord(l1) = coord(l1) + coord(l2) - pi * fpint(l2)
      end do

!  find the interval for which fpint is maximal on the condition that
!  there still can be added a knot.
      l1 = 1
      l2 = nuu + nrr
      if (nuest < nu + 1) l1 = nuu + 1
      if (nvest < nv + 2) l2 = nuu
!  test whether we cannot further increase the number of knots.
      if (l1 > l2) go to 950
500   fpmax = 0.
      l = 0
      do i = l1, l2
        if (fpmax >= fpint(i)) cycle
        l = i
        fpmax = fpint(i)
      end do

      if (l == 0) go to 930
!  calculate the position of the new knot.
      arg = coord(l) / fpint(l)
!  test in what direction the new knot is going to be added.
      if (l > nuu) go to 530
!  addition in the u-direction
      l4 = l + 4
      fpint(l) = 0.
      fac1 = tu(l4) - arg
      fac2 = arg - tu(l4 - 1)
      if (fac1 > (ten * fac2) .or. fac2 > (ten * fac1)) go to 500
      j = nu
      do i = l4, nu
        tu(j + 1) = tu(j)
        j = j - 1
      end do

      tu(l4) = arg
      nu = nu + 1
      cycle
!  addition in the v-direction
530   l4 = l + 4 - nuu
      fpint(l) = 0.
      fac1 = tv(l4) - arg
      fac2 = arg - tv(l4 - 1)
      if (fac1 > (ten * fac2) .or. fac2 > (ten * fac1)) go to 500
      ll = nrr + 4
      j = ll
      do i = l4, ll
        tv(j + 1) = tv(j)
        j = j - 1
      end do

      tv(l4) = arg
      nv = nv + 2
      nrr = nrr + 1
      do i = 5, ll
        j = i + nrr
        tv(j) = tv(i) + pi
      end do

!  restart the computations with the new set of knots.
    end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 2: determination of the smoothing bicubic spline.               c
! ******************************************************               c
! we have determined the number of knots and their position. we now    c
! compute the coefficients of the smoothing spline sp(u,v).            c
! the observation matrix a is extended by the rows of a matrix, expres-c
! sing that sp(u,v) must be a constant function in the variable        c
! v and a cubic polynomial in the variable u. the corresponding        c
! weights of these additional rows are set to 1/(p). iteratively       c
! we than have to determine the value of p such that f(p) = sum((w(i)* c
! (z(i)-sp(u(i),v(i))))**2)  be = s.                                   c
! we already know that the least-squares polynomial corresponds to p=0,c
! and that the least-squares bicubic spline corresponds to p=infin.    c
! the iteration process makes use of rational interpolation. since f(p)c
! is a convex and strictly decreasing function of p, it can be approx- c
! imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c
! three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c
! f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c
! of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  evaluate the discontinuity jumps of the 3-th order derivative of
!  the b-splines at the knots tu(l),l=5,...,nu-4.
580 call fpdisc(tu, nu, 5, bu, nuest)
!  evaluate the discontinuity jumps of the 3-th order derivative of
!  the b-splines at the knots tv(l),l=5,...,nv-4.
    call fpdisc(tv, nv, 5, bv, nvest)
!  initial value for p.
    p1 = 0.
    f1 = sup - s
    p3 = -one
    f3 = fpms
    p = sum(a(:ncof, 1))
    rn = ncof
    p = rn / p
!  find the bandwidth of the extended observation matrix.
    iband4 = iband + ipar1
    if (iband4 > ncof) iband4 = ncof
    iband3 = iband4 - 1
    ich1 = 0
    ich3 = 0
    nuu = nu4 - iopt3 - 1
!  iteration process to find the root of f(p)=s.
    do 920 iter = 1, maxit
      pinv = one / p
!  store the triangularized observation matrix into q.

      ff(:ncof) = f(:ncof)
      ! JF: Here overwrite some components

      q(:ncof, :iband4) = 0._8
      q(:ncof, :iband) = a(:ncof, :iband)

!  extend the observation matrix with the rows of a matrix, expressing
!  that for u=constant sp(u,v) must be a constant function.
      do i = 5, nv4
        ii = i - 4
        row(:nvv) = 0._8

        ll = ii
        do l = 1, 5
          if (ll > nvv) ll = 1
          row(ll) = row(ll) + bv(ii, l)
          ll = ll + 1
        end do

        do j = 1, nuu
!  initialize the new row.

          h(:iband) = 0._8

!  fill in the non-zero elements of the row. jrot records the column
!  number of the first non-zero element in the row.
          if (j > iopt2) go to 665
          if (j == 2) go to 655
          do k = 1, 2
            cs(k) = 0.
            do l = 1, nvv
              cs(k) = cs(k) + cosi(k, l) * row(l)
            end do
          end do

          h(1) = cs(1)
          h(2) = cs(2)
          jrot = 2
          go to 675
655       do k = 3, 5
            cs(k) = 0.
            do l = 1, nvv
              cs(k) = cs(k) + cosi(k, l) * row(l)
            end do
          end do

          h(1:2) = cs(1:2) * ratio
          h(3:5) = cs(3:5)
          jrot = 2
          go to 675
665       h(:nvv) = row(:nvv)
          jrot = ipar1 + 1 + (j - iopt2 - 1) * nvv
675       h(:iband) = h(:iband) * pinv

          zi = 0.
!  rotate the new row into triangle by givens transformations.
          inner1: do irot = jrot, ncof
            piv = h(1)
            i2 = min0(iband1, ncof - irot)
            if (piv == 0.) then
              if (i2 <= 0) exit inner1
              go to 690
            endif
!  calculate the parameters of the givens transformation.
            call fpgivs(piv, q(irot, 1), co, si)
!  apply that givens transformation to the right hand side.
            call fprota(co, si, zi, ff(irot))
            if (i2 == 0) exit inner1
!  apply that givens transformation to the left hand side.
            do l = 1, i2
              l1 = l + 1
              call fprota(co, si, h(l1), q(irot, l1))
            end do

690         do l = 1, i2
              h(l) = h(l + 1)
            end do

            h(i2 + 1) = 0.
          end do inner1

        end do
      end do

!  extend the observation matrix with the rows of a matrix expressing
!  that for v=constant. sp(u,v) must be a cubic polynomial.
      do i = 5, nu4
        ii = i - 4
        do j = 1, nvv
          !  initialize the new row
          h(:iband4) = 0._8
!  fill in the non-zero elements of the row. jrot records the column
!  number of the first non-zero element in the row.
          j1 = 1
          do l = 1, 5
            il = ii + l - 1
            if (il == nu4 .and. iopt3 /= 0) cycle
            if (il > iopt2 + 1) go to 750
            if (il == 1) then
              go to 735
            else if (il == 2) then
              go to 740
            else if (il == 3) then
              go to 745
            end if

735         h(1) = bu(ii, l)
            j1 = j + 1
            cycle
740         h(1) = h(1) + bu(ii, l)
            h(2) = bu(ii, l) * cosi(1, j)
            h(3) = bu(ii, l) * cosi(2, j)
            j1 = j + 3
            cycle
745         h(1) = h(1) + bu(ii, l)
            h(2) = bu(ii, l) * cosi(1, j) * ratio
            h(3) = bu(ii, l) * cosi(2, j) * ratio
            h(4) = bu(ii, l) * cosi(3, j)
            h(5) = bu(ii, l) * cosi(4, j)
            h(6) = bu(ii, l) * cosi(5, j)
            j1 = j + 6
            cycle
750         h(j1) = bu(ii, l)
            j1 = j1 + nvv
          end do

          h(:iband4) = h(:iband4) * pinv

          zi = 0.
          jrot = 1
          if (ii > iopt2 + 1) jrot = ipar1 + (ii - iopt2 - 2) * nvv + j
!  rotate the new row into triangle by givens transformations.
          do irot = jrot, ncof
            piv = h(1)
            i2 = min0(iband3, ncof - irot)
            if (piv == 0.) then
              if (i2 <= 0) go to 810
              go to 780
            endif
!  calculate the parameters of the givens transformation.
            call fpgivs(piv, q(irot, 1), co, si)
!  apply that givens transformation to the right hand side.
            call fprota(co, si, zi, ff(irot))
            if (i2 == 0) go to 810
!  apply that givens transformation to the left hand side.
            do l = 1, i2
              l1 = l + 1
              call fprota(co, si, h(l1), q(irot, l1))
            end do

780         continue
            do l = 1, i2
              h(l) = h(l + 1)
            end do
            h(i2 + 1) = 0.
          end do

810       continue
        end do
      end do

!  find dmax, the maximum value for the diagonal elements in the
!  reduced triangle.
      dmax = 0.
      do 820 i = 1, ncof
        if (q(i, 1) <= dmax) go to 820
        dmax = q(i, 1)
820     continue
!  check whether the matrix is rank deficient.
        sigma = eps * dmax
        do 830 i = 1, ncof
          if (q(i, 1) <= sigma) go to 840
830       continue
!  backward substitution in case of full rank.
          call fpback(q, ff, ncof, iband4, c, ncc)
          rank = ncof
          go to 845
!  in case of rank deficiency, find the minimum norm solution.
840       lwest = ncof * iband4 + ncof + iband4
          if (lwrk < lwest) go to 925
          lf = 1
          lh = lf + ncof
          la = lh + iband4
          call fprank(q, ff, ncof, iband4, ncc, sigma, c, sq, rank, wrk(la), wrk(lf), wrk(lh))
845       q(:ncof, 1) = q(:ncof, 1) / dmax

!  find the coefficients in the standard b-spline representation of
!  the polar spline.
          call fprppo(nu, nv, iopt2, iopt3, cosi, ratio, c, ff, ncoff)
!  compute f(p).
          fp = 0.
          do 890 num = 1, nreg
            num1 = num - 1
            lu = num1 / nvv
            lv = num1 - lu * nvv
            jrot = lu * nv4 + lv
            in = index(num)
860         if (in == 0) go to 890
            store = 0.
            i1 = jrot
            do 880 i = 1, 4
              hui = spu(in, i)
              j1 = i1
              do 870 j = 1, 4
                j1 = j1 + 1
                store = store + hui * spv(in, j) * c(j1)
870             continue
                i1 = i1 + nv4
880             continue
                fp = fp + (w(in) * (z(in) - store))**2
                in = nummer(in)
                go to 860
890             continue
!  test whether the approximation sp(u,v) is an acceptable solution
                fpms = fp - s
                if (abs(fpms) <= acc) go to 980
!  test whether the maximum allowable number of iterations has been
!  reached.
                if (iter == maxit) go to 940
!  carry out one more step of the iteration process.
                p2 = p
                f2 = fpms
                if (ich3 /= 0) go to 900
                if ((f2 - f3) > acc) go to 895
!  our initial choice of p is too large.
                p3 = p2
                f3 = f2
                p = p * con4
                if (p <= p1) p = p1 * con9 + p2 * con1
                go to 920
895             if (f2 < 0.) ich3 = 1
900             if (ich1 /= 0) go to 910
                if ((f1 - f2) > acc) go to 905
!  our initial choice of p is too small
                p1 = p2
                f1 = f2
                p = p / con4
                if (p3 < 0.) go to 920
                if (p >= p3) p = p2 * con1 + p3 * con9
                go to 920
905             if (f2 > 0.) ich1 = 1
!  test whether the iteration process proceeds as theoretically
!  expected.
910             if (f2 >= f1 .or. f2 <= f3) go to 945
!  find the new value of p.
                p = fprati(p1, f1, p2, f2, p3, f3)
920             continue
!  error codes and messages.
925             ier = lwest
                go to 990
930             ier = 5
                go to 990
935             ier = 4
                go to 990
940             ier = 3
                go to 990
945             ier = 2
                go to 990
950             ier = 1
                go to 990
960             ier = -2
                go to 990
970             ier = -1
                fp = 0.
980             if (ncof /= rank) ier = -rank
990             return
              end

