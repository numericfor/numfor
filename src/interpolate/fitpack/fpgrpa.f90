subroutine fpgrpa(ifsu, ifsv, ifbu, ifbv, idim, ipar, u, mu, v, mv, z, mz,&
  &tu, nu, tv, nv, p, c, nc, fp, fpu, fpv, mm, mvnu, spu, spv, right, q, au, au1,&
  &av, av1, bu, bv, nru, nrv)
  !  ..
  !  ..scalar arguments..
  real(8) :: p, fp
  integer :: ifsu, ifsv, ifbu, ifbv, idim, mu, mv, mz, nu, nv, nc, mm, mvnu
  !  ..array arguments..
  real(8) :: u(mu), v(mv), z(mz * idim), tu(nu), tv(nv), c(nc * idim), fpu(nu)
  real(8) :: fpv(nv), spu(mu, 4), spv(mv, 4), right(mm * idim), q(mvnu), au(nu, 5)
  real(8) :: au1(nu, 4), av(nv, 5), av1(nv, 4), bu(nu, 5), bv(nv, 5)
  integer :: ipar(2), nru(mu), nrv(mv)
  !  ..local scalars..
  real(8) :: arg, fac, term, one, half, value
  integer :: i, id, ii, it, iz, i1, i2, j, jz, k, k1, k2, l, l1, l2, mvv, k0, muu
  integer :: ncof, nroldu, nroldv, number, nmd, numu, numu1, numv, numv1, nuu, nvv
  integer :: nu4, nu7, nu8, nv4, nv7, nv8
  !  ..local arrays..
  real(8) :: h(5)
  !  ..subroutine references..
  !    fpback,fpbspl,fpdisc,fpbacp,fptrnp,fptrpe
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
  !    (2)  c(nu-3+i,j) = c(i,j), i=1,2,3 ; j=1,2,...,nv-4
  !            if(ipar(1) /= 0)
  !
  !    (3)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
  !            if(ipar(2) /= 0)
  !
  !  set constants
  one = 1
  half = 0.5
  !  initialization
  nu4 = nu - 4
  nu7 = nu - 7
  nu8 = nu - 8
  nv4 = nv - 4
  nv7 = nv - 7
  nv8 = nv - 8
  muu = mu
  if (ipar(1) /= 0) muu = mu - 1
  mvv = mv
  if (ipar(2) /= 0) mvv = mv - 1
  !  it depends on the value of the flags ifsu,ifsv,ifbu  and ibvand
  !  on the value of p whether the matrices (spu), (spv), (bu) and (bv)
  !  still must be determined.
  if (ifsu == 0) then
    !  calculate the non-zero elements of the matrix (spu) which is the ob-
    !  servation matrix according to the least-squares spline approximation
    !  problem in the u-direction.
    l = 4
    l1 = 5
    number = 0
    do it = 1, muu
      arg = u(it)
      do while (tu(l1) < arg .and. l /= nu4)
        l = l1
        l1 = l + 1
        number = number + 1
      end do
      call fpbspl(tu, nu, 3, arg, l, h)
      spu(it, :4) = h(:4)
      nru(it) = number
    end do

    ifsu = 1
  end if

  !  calculate the non-zero elements of the matrix (spv) which is the ob-
  !  servation matrix according to the least-squares spline approximation
  !  problem in the v-direction.
  if (ifsv /= 0) go to 100
  l = 4
  l1 = 5
  number = 0
  do it = 1, mvv
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
100 if (p <= 0._8) go to 150
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

  !  substituting (2)  and (3) into (1), we obtain the overdetermined
  !  system
  !         (4)  (avv) (cr) (auu)' = (qq)
  !  from which the nuu*nvv remaining coefficients
  !         c(i,j) , i=1,...,nu-4-3*ipar(1) ; j=1,...,nv-4-3*ipar(2) ,
  !  the elements of (cr), are then determined in the least-squares sense.
  !  we first determine the matrices (auu) and (qq). then we reduce the
  !  matrix (auu) to upper triangular form (ru) using givens rotations.
  !  we apply the same transformations to the rows of matrix qq to obtain
  !  the (mv) x nuu matrix g.
  !  we store matrix (ru) into au (and au1 if ipar(1)=1) and g into q.
150 if (ipar(1) == 0) then
    nuu = nu4
    call fptrnp(mu, mv, idim, nu, nru, spu, p, bu, z, au, q, right)
  else
    nuu = nu7
    call fptrpe(mu, mv, idim, nu, nru, spu, p, bu, z, au, au1, q, right)
  end if
  !  we determine the matrix (avv) and then we reduce this matrix to
  !  upper triangular form (rv) using givens rotations.
  !  we apply the same transformations to the columns of matrix
  !  g to obtain the (nvv) x (nuu) matrix h.
  !  we store matrix (rv) into av (and av1 if ipar(2)=1) and h into c.
  if (ipar(2) == 0) then
    nvv = nv4
    call fptrnp(mv, nuu, idim, nv, nrv, spv, p, bv, q, av, c, right)
  else
    nvv = nv7
    call fptrpe(mv, nuu, idim, nv, nrv, spv, p, bv, q, av, av1, c, right)
  end if

  !  backward substitution to obtain the b-spline coefficients as the
  !  solution of the linear system    (rv) (cr) (ru)' = h.
  !  first step: solve the system  (rv) (c1) = h.
  ncof = nuu * nvv
  k = 1
  if (ipar(2) == 0) then
    do ii = 1, idim
      do i = 1, nuu
        call fpback(av, c(k), nvv, 5, c(k), nv)
        k = k + nvv
      end do
    end do
  else
    do ii = 1, idim
      do i = 1, nuu
        call fpbacp(av, av1, c(k), nvv, 4, c(k), 5, nv)
        k = k + nvv
      end do
    end do
  end if

  !  second step: solve the system  (cr) (ru)' = (c1).
  if (ipar(1) == 0) then
    do ii = 1, idim
      k = (ii - 1) * ncof
      do j = 1, nvv
        k = k + 1
        l = k
        do i = 1, nuu
          right(i) = c(l)
          l = l + nvv
        end do

        call fpback(au, right, nuu, 5, right, nu)
        l = k
        do i = 1, nuu
          c(l) = right(i)
          l = l + nvv
        end do

      end do
    end do

  else
    do ii = 1, idim
      k = (ii - 1) * ncof
      do j = 1, nvv
        k = k + 1
        l = k
        do i = 1, nuu
          right(i) = c(l)
          l = l + nvv
        end do

        call fpbacp(au, au1, right, nuu, 4, right, 5, nu)
        l = k
        do i = 1, nuu
          c(l) = right(i)
          l = l + nvv
        end do

      end do
    end do

  end if

  !  calculate from the conditions (2)-(3), the remaining b-spline
  !  coefficients.
  if (ipar(2) == 0) go to 600
  i = 0
  j = 0
  do id = 1, idim

    do l = 1, nuu
      ii = i
      do k = 1, nvv
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
  end do

  ncof = nv4 * nuu
  nmd = ncof * idim
  c(:nmd) = q(:nmd)
600 if (ipar(1) == 0) go to 700
  i = 0
  j = 0
  n33 = 3 * nv4
  do id = 1, idim
    ii = i
    do k = 1, ncof
      i = i + 1
      j = j + 1
      q(i) = c(j)
    end do

    do k = 1, n33
      ii = ii + 1
      i = i + 1
      q(i) = q(ii)
    end do
  end do

  ncof = nv4 * nu4
  nmd = ncof * idim
  c(:nmd) = q(:nmd)
  !  calculate the quantities
  !    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
  !    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
  !    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
  !                  tu(r+3) <= u(i) <= tu(r+4)
  !    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
  !                  tv(r+3) <= v(j) <= tv(r+4)
700 fp = 0._8
  fpu(:nu) = 0._8
  fpv(:nv) = 0._8
  nroldu = 0
  !  main loop for the different grid points.
  do i1 = 1, muu
    numu = nru(i1)
    numu1 = numu + 1
    nroldv = 0
    iz = (i1 - 1) * mv
    do i2 = 1, mvv
      numv = nrv(i2)
      numv1 = numv + 1
      iz = iz + 1
      !  evaluate s(u,v) at the current grid point by making the sum of the
      !  cross products of the non-zero b-splines at (u,v), multiplied with
      !  the appropriate b-spline coefficients.
      term = 0.
      k0 = numu * nv4 + numv
      jz = iz
      do id = 1, idim
        k1 = k0
        value = 0.
        do l1 = 1, 4
          k2 = k1
          fac = spu(i1, l1)
          do l2 = 1, 4
            k2 = k2 + 1
            value = value + fac * spv(i2, l2) * c(k2)
          end do

          k1 = k1 + nv4
        end do

        !  calculate the squared residual at the current grid point.
        term = term + (z(jz) - value)**2
        jz = jz + mz
        k0 = k0 + ncof
      end do

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

end subroutine fpgrpa
