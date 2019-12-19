subroutine fpgrdi(ifsu, ifsv, ifbu, ifbv, iback, u, mu, v, mv, z, mz, dz,&
  &iop0, iop1, tu, nu, tv, nv, p, c, nc, sq, fp, fpu, fpv, mm, mvnu, spu, spv,&
  &right, q, au, av1, av2, bu, bv, aa, bb, cc, cosi, nru, nrv)
  !  ..
  !  ..scalar arguments..
  real(8) :: p, sq, fp
  integer :: ifsu, ifsv, ifbu, ifbv, iback, mu, mv, mz, iop0, iop1, nu, nv, nc, mm, mvnu
  !  ..array arguments..
  real(8) :: u(mu), v(mv), z(mz), dz(3), tu(nu), tv(nv), c(nc), fpu(nu), fpv(nv)
  real(8) :: spu(mu, 4), spv(mv, 4), right(mm), q(mvnu), au(nu, 5), av1(nv, 6)
  real(8) :: av2(nv, 4), aa(2, mv), bb(2, nv), cc(nv), cosi(2, nv), bu(nu, 5), bv(nv, 5)
  integer :: nru(mu), nrv(mv)
  !  ..local scalars..
  real(8) :: arg, co, dz1, dz2, dz3, fac, fac0, pinv, piv, si, term, one, three, half
  integer :: i, ic, ii, ij, ik, iq, irot, it, iz, i0, i1, i2, i3, j, jk, jper
  integer :: j0, j1, k, k1, k2, l, l0, l1, l2, mvv, ncof, nrold, nroldu, nroldv, number
  integer :: numu, numu1, numv, numv1, nuu, nu4, nu7, nu8, nu9, nv11, nv4, nv7, nv8, n1
  !  ..local arrays..
  real(8) :: h(5), h1(5), h2(4)
  !  ..function references..
  integer :: min0
  real(8) :: cos, sin
  !  ..subroutine references..
  !    fpback,fpbspl,fpgivs,fpcyt1,fpcyt2,fpdisc,fpbacp,fprota
  !  ..
  !  let
  !               |   (spu)    |            |   (spv)    |
  !        (au) = | ---------- |     (av) = | ---------- |
  !               | (1/p) (bu) |            | (1/p) (bv) |
  !
  !                                | z  ' 0 |
  !                            q = | ------ |
  !                                | 0  ' 0 |
  !
  !  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
  !                coefficients.
  !       z      : the mu x mv matrix which contains the function values.
  !       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
  !                according to the least-squares problems in the u-,resp.
  !                v-direction.
  !       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
  !                containing the discontinuity jumps of the derivatives
  !                of the b-splines in the u-,resp.v-variable at the knots
  !  the b-spline coefficients of the smoothing spline are then calculated
  !  as the least-squares solution of the following over-determined linear
  !  system of equations
  !
  !    (1)  (av) c (au)' = q
  !
  !  subject to the constraints
  !
  !    (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
  !
  !    (3)  if iop0 = 0  c(1,j) = dz(1)
  !            iop0 = 1  c(1,j) = dz(1)
  !                      c(2,j) = dz(1)+(dz(2)*cosi(1,j)+dz(3)*cosi(2,j))*
  !                               tu(5)/3. = cc(j) , j=1,2,...nv-4
  !
  !    (4)  if iop1 = 1  c(nu-4,j) = 0, j=1,2,...,nv-4.
  !
  !  set constants
  one = 1
  three = 3
  half = 0.5_8
  !  initialization
  nu4 = nu - 4
  nu7 = nu - 7
  nu8 = nu - 8
  nu9 = nu - 9
  nv4 = nv - 4
  nv7 = nv - 7
  nv8 = nv - 8
  nv11 = nv - 11
  nuu = nu4 - iop0 - iop1 - 1
  if (p > 0.) pinv = one / p
  !  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv and iop0 and
  !  on the value of p whether the matrices (spu), (spv), (bu), (bv) and
  !  (cosi) still must be determined.
  if (ifsu /= 0) go to 30
  !  calculate the non-zero elements of the matrix (spu) which is the ob-
  !  servation matrix according to the least-squares spline approximation
  !  problem in the u-direction.
  l = 4
  l1 = 5
  number = 0
  do it = 1, mu
    arg = u(it)
10  if (arg < tu(l1) .or. l == nu4) go to 15
    l = l1
    l1 = l + 1
    number = number + 1
    go to 10
15  call fpbspl(tu, nu, 3, arg, l, h)
    spu(it, :4) = h(:4)
    nru(it) = number
  end do

  ifsu = 1
  !  calculate the non-zero elements of the matrix (spv) which is the ob-
  !  servation matrix according to the least-squares spline approximation
  !  problem in the v-direction.
30 if (ifsv /= 0) go to 85
  l = 4
  l1 = 5
  number = 0
  do it = 1, mv
    arg = v(it)
    do while (tv(l1) <= arg .and. l /= nv4)
      l = l1
      l1 = l + 1
      number = number + 1
    end do
    call fpbspl(tv, nv, 3, arg, l, h)
    spv(it, :4) = h(:4)
    nrv(it) = number
  end do

  ifsv = 1
  if (iop0 == 0) go to 85
  !  calculate the coefficients of the interpolating splines for cos(v)
  !  and sin(v).
  cosi(:, :nv4) = 0._8
  if (nv7 < 4) go to 85
  do i = 1, nv7
    l = i + 3
    arg = tv(l)
    call fpbspl(tv, nv, 3, arg, l, h)
    av1(i, :3) = h(:3)
    cosi(1, i) = cos(arg)
    cosi(2, i) = sin(arg)
  end do

  call fpcyt1(av1, nv7, nv)
  do j = 1, 2
    right(:nv7) = cosi(j, :nv7)
    call fpcyt2(av1, nv7, right, right, nv)
    cosi(j, 2:nv7 + 1) = right(:nv7)
    cosi(j, 1) = cosi(j, nv7 + 1)
    cosi(j, nv7 + 2) = cosi(j, 2)
    cosi(j, nv4) = cosi(j, 3)
  end do

85 if (p <= 0._8) go to 150
  !  calculate the non-zero elements of the matrix (bu).
  if (ifbu == 0 .and. nu8 /= 0) then
    call fpdisc(tu, nu, 5, bu, nu)
    ifbu = 1
  end if

  !  calculate the non-zero elements of the matrix (bv).
  if (ifbv == 0 .and. nv8 /= 0) then
    call fpdisc(tv, nv, 5, bv, nv)
    ifbv = 1
  end if

  !  substituting (2),(3) and (4) into (1), we obtain the overdetermined
  !  system
  !         (5)  (avv) (cr) (auu)' = (qq)
  !  from which the nuu*nv7 remaining coefficients
  !         c(i,j) , i=2+iop0,3+iop0,...,nu-4-iop1 ; j=1,2,...,nv-7 ,
  !  the elements of (cr), are then determined in the least-squares sense.
  !  simultaneously, we compute the resulting sum of squared residuals sq.
150 dz1 = dz(1)
  aa(1, 1:mv) = dz1
  if (nv8 /= 0 .and. p > 0._8) bb(1, :nv8) = 0._8
  mvv = mv
  if (iop0 == 0) go to 220
  fac = tu(5) / three
  dz2 = dz(2) * fac
  dz3 = dz(3) * fac
  do i = 1, nv4
    cc(i) = dz1 + dz2 * cosi(1, i) + dz3 * cosi(2, i)
  end do

  do i = 1, mv
    number = nrv(i)
    fac = 0.
    do j = 1, 4
      number = number + 1
      fac = fac + cc(number) * spv(i, j)
    end do

    aa(2, i) = fac
  end do

  if (nv8 /= 0 .and. p > 0.) then
    do i = 1, nv8
      number = i
      fac = 0.
      do j = 1, 5
        fac = fac + cc(number) * bv(i, j)
        number = number + 1
      end do

      bb(2, i) = fac * pinv
    end do

    mvv = mvv + nv8
  end if
  !  we first determine the matrices (auu) and (qq). then we reduce the
  !  matrix (auu) to upper triangular form (ru) using givens rotations.
  !  we apply the same transformations to the rows of matrix qq to obtain
  !  the (mv+nv8) x nuu matrix g.
  !  we store matrix (ru) into au and g into q.
220 l = mvv * nuu
  !  initialization.
  sq = 0.
  q(:l) = 0._8
  au(:nuu, :5) = 0._8
  l = 0
  nrold = 0
  n1 = nrold + 1
  do it = 1, mu
    number = nru(it)
    !  find the appropriate column of q.
250 right(:mvv) = 0._8
    if (nrold == number) go to 280
    if (p <= 0.) go to 410
    !  fetch a new row of matrix (bu).
    h(:5) = bu(n1, :5) * pinv
    i0 = 1
    i1 = 5
    go to 310
    !  fetch a new row of matrix (spu).
280 h(:4) = spu(it, :4)
    !  find the appropriate column of q.
    right(:mv) = z(l + 1:l + mv)
    ! l = l + mv
    i0 = 1
    i1 = 4
310 if (nu7 - number == iop1) i1 = i1 - 1
    j0 = n1
    !  take into account that we eliminate the constraints (3)
    do while (j0 - 1 <= iop0)
      fac0 = h(i0)
      right(:mv) = right(:mv) - fac0 * aa(j0, :mv)
      if (mv /= mvv) right(mv + 1:mv + nv8) = right(mv + 1:mv + nv8) - fac0 * bb(j0, 1:nv8)

      j0 = j0 + 1
      i0 = i0 + 1
    end do

    irot = nrold - iop0 - 1
    if (irot < 0) irot = 0
    !  rotate the new row of matrix (auu) into triangle.
    do i = i0, i1
      irot = irot + 1
      piv = h(i)
      if (piv == 0._8) cycle
      !  calculate the parameters of the givens transformation.
      call fpgivs(piv, au(irot, 1), co, si)
      !  apply that transformation to the rows of matrix (qq).
      iq = (irot - 1) * mvv
      do j = 1, mvv
        iq = iq + 1
        call fprota(co, si, right(j), q(iq))
      end do

      !  apply that transformation to the columns of (auu).
      if (i == i1) cycle
      i2 = 1
      i3 = i + 1
      do j = i3, i1
        i2 = i2 + 1
        call fprota(co, si, h(j), au(irot, i2))
      end do

    end do

    ! we update the sum of squared residuals
    do j = 1, mvv
      sq = sq + right(j)**2
    end do

    if (nrold == number) cycle
410 nrold = n1
    n1 = n1 + 1
    go to 250
  end do

  ncof = nuu * nv7
  !  initialization.
  c(:ncof) = 0._8
  av1(:nv4, :5) = 0._8
  av2(:nv4, :4) = 0._8
  nrold = 0
  do it = 1, mv
    number = nrv(it)
450 if (nrold == number) then
      !  fetch a new row of matrix (spv)
      h(5) = 0._8
      h(:4) = spv(it, :4)
      !  find the appropriate row of g.
      l = it
      do j = 1, nuu
        right(j) = q(l)
        l = l + mvv
      end do
    else
      if (p <= 0._8) go to 760
      !  fetch a new row of matrix (bv).
      n1 = nrold + 1
      h(:5) = bv(n1, :5) * pinv
      !  find the appropriate row of g.
      right(:nuu) = 0._8
      if (mv /= mvv) then
        l = mv + n1
        do j = 1, nuu
          right(j) = q(l)
          l = l + mvv
        end do
      end if
    end if

    !  test whether there are non-zero values in the new row of (avv)
    !  corresponding to the b-splines n(j,v),j=nv7+1,...,nv4.
    if (nrold < nv11) go to 710
    if (jper == 0) then
      !  initialize the matrix (av2).
      jk = nv11 + 1
      do i = 1, 4
        ik = jk
        do j = 1, 5
          if (ik <= 0) exit
          av2(ik, i) = av1(ik, j)
          ik = ik - 1
        end do
        jk = jk + 1
      end do

      jper = 1
    end if

    !  if one of the non-zero elements of the new row corresponds to one of
    !  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition
    !  (2) for setting up this row of (avv). the row is stored in h1( the
    !  part with respect to av1) and h2 (the part with respect to av2).
    h1(:5) = 0._8
    h2(:4) = 0._8
    j = nrold - nv11

    do i = 1, 5
      j = j + 1
      l0 = j
570   l1 = l0 - 4
      if (l1 <= 0) then
        h2(l0) = h2(l0) + h(i)
      else
        if (l1 <= nv11) go to 580
        l0 = l1 - nv11
        go to 570
580     h1(l1) = h(i)
      end if
    end do

    !  rotate the new row of (avv) into triangle.
    if (nv11 <= 0) go to 670
    !  rotations with the rows 1,2,...,nv11 of (avv).
    do j = 1, nv11
      piv = h1(1)
      i2 = min0(nv11 - j, 4)
      if (piv == 0.) go to 640
      !  calculate the parameters of the givens transformation.
      call fpgivs(piv, av1(j, 1), co, si)
      !  apply that transformation to the columns of matrix g.
      ic = j
      do i = 1, nuu
        call fprota(co, si, right(i), c(ic))
        ic = ic + nv7
      end do

      !  apply that transformation to the rows of (avv) with respect to av2.
      do i = 1, 4
        call fprota(co, si, h2(i), av2(j, i))
      end do

      !  apply that transformation to the rows of (avv) with respect to av1.
      if (i2 == 0) go to 670
      do i = 1, i2
        i1 = i + 1
        call fprota(co, si, h1(i1), av1(j, i1))
      end do

640   do i = 1, i2
        h1(i) = h1(i + 1)
      end do
      h1(i2 + 1) = 0.
    end do

    !  rotations with the rows nv11+1,...,nv7 of avv.
670 do j = 1, 4
      ij = nv11 + j
      if (ij <= 0) cycle
      piv = h2(j)
      if (piv == 0._8) cycle
      !  calculate the parameters of the givens transformation.
      call fpgivs(piv, av2(ij, j), co, si)
      !  apply that transformation to the columns of matrix g.
      ic = ij
      do i = 1, nuu
        call fprota(co, si, right(i), c(ic))
        ic = ic + nv7
      end do

      if (j == 4) cycle
      !  apply that transformation to the rows of (avv) with respect to av2.
      j1 = j + 1
      do i = j1, 4
        call fprota(co, si, h2(i), av2(ij, i))
      end do

    end do

    ! we update the sum of squared residuals
    do i = 1, nuu
      sq = sq + right(i)**2
    end do

    go to 750
    !  rotation into triangle of the new row of (avv), in case the elements
    !  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
710 irot = nrold
    do i = 1, 5
      irot = irot + 1
      piv = h(i)
      if (piv == 0._8) cycle
      !  calculate the parameters of the givens transformation.
      call fpgivs(piv, av1(irot, 1), co, si)
      !  apply that transformation to the columns of matrix g.
      ic = irot
      do j = 1, nuu
        call fprota(co, si, right(j), c(ic))
        ic = ic + nv7
      end do

      !  apply that transformation to the rows of (avv).
      if (i == 5) cycle
      i2 = 1
      i3 = i + 1
      do j = i3, 5
        i2 = i2 + 1
        call fprota(co, si, h(j), av1(irot, i2))
      end do
    end do

    ! we update the sum of squared residuals
    do i = 1, nuu
      sq = sq + right(i)**2
    end do

750 if (nrold == number) cycle
760 nrold = nrold + 1
    go to 450
  end do

  !  test whether the b-spline coefficients must be determined.
  if (iback /= 0) return
  !  backward substitution to obtain the b-spline coefficients as the
  !  solution of the linear system    (rv) (cr) (ru)' = h.
  !  first step: solve the system  (rv) (c1) = h.
  k = 1
  do i = 1, nuu
    call fpbacp(av1, av2, c(k), nv7, 4, c(k), 5, nv)
    k = k + nv7
  end do

  !  second step: solve the system  (cr) (ru)' = (c1).
  k = 0
  do j = 1, nv7
    k = k + 1
    l = k
    do i = 1, nuu
      right(i) = c(l)
      l = l + nv7
    end do

    call fpback(au, right, nuu, 5, right, nu)
    l = k
    do i = 1, nuu
      c(l) = right(i)
      l = l + nv7
    end do

  end do

  !  calculate from the conditions (2)-(3)-(4), the remaining b-spline
  !  coefficients.
  ncof = nu4 * nv4
  i = nv4
  j = 0
  q(:nv4) = dz1
  if (iop0 == 0) go to 815
  do l = 1, nv4
    i = i + 1
    q(i) = cc(l)
  end do

815 if (nuu /= 0) then
    do l = 1, nuu
      ii = i
      do k = 1, nv7
        i = i + 1
        j = j + 1
        q(i) = c(j)
      end do

      do k = 1, 3
        ii = ii + 1
        i = i + 1
        q(i) = q(ii)
      end do
    end do
  end if

  if (iop1 /= 0) q(i + 1:i + nv4) = 0._8
  c(:ncof) = q(:ncof)
  !  calculate the quantities
  !    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
  !    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
  !    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
  !                  tu(r+3) <= u(i) <= tu(r+4)
  !    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
  !                  tv(r+3) <= v(j) <= tv(r+4)
  fp = 0._8
  fpu(:nu) = 0._8
  fpv(:nv) = 0._8
  iz = 0
  nroldu = 0
  !  main loop for the different grid points.
  do i1 = 1, mu
    numu = nru(i1)
    numu1 = numu + 1
    nroldv = 0
    do i2 = 1, mv
      numv = nrv(i2)
      numv1 = numv + 1
      iz = iz + 1
      !  evaluate s(u,v) at the current grid point by making the sum of the
      !  cross products of the non-zero b-splines at (u,v), multiplied with
      !  the appropriate b-spline coefficients.
      term = 0.
      k1 = numu * nv4 + numv
      do l1 = 1, 4
        k2 = k1
        fac = spu(i1, l1)
        do l2 = 1, 4
          k2 = k2 + 1
          term = term + fac * spv(i2, l2) * c(k2)
        end do

        k1 = k1 + nv4
      end do

      !  calculate the squared residual at the current grid point.
      term = (z(iz) - term)**2
      !  adjust the different parameters.
      fp = fp + term
      fpu(numu1) = fpu(numu1) + term
      fpv(numv1) = fpv(numv1) + term
      fac = term * half
      if (numv /= nroldv) then
        fpv(numv1) = fpv(numv1) - fac
        fpv(numv) = fpv(numv) + fac
      end if
      nroldv = numv
      if (numu == nroldu) cycle
      fpu(numu1) = fpu(numu1) - fac
      fpu(numu) = fpu(numu) + fac
    end do

    nroldu = numu
  end do

end subroutine fpgrdi
