subroutine fpgrre(ifsx, ifsy, ifbx, ifby, x, mx, y, my, z, mz, kx, ky, tx, nx,&
  &ty, ny, p, c, nc, fp, fpx, fpy, mm, mynx, kx1, kx2, ky1, ky2, spx, spy, right, q,&
  &ax, ay, bx, by, nrx, nry)
  !  ..
  !  ..scalar arguments..
  real(8) :: p, fp
  integer :: ifsx, ifsy, ifbx, ifby, mx, my, mz, kx, ky, nx, ny, nc, mm, mynx,&
    &kx1, kx2, ky1, ky2
  !  ..array arguments..
  real(8) :: x(mx), y(my), z(mz), tx(nx), ty(ny), c(nc), spx(mx, kx1), spy(my, ky1)
  real(8) :: right(mm), q(mynx), ax(nx, kx2), bx(nx, kx2), ay(ny, ky2), by(ny, ky2), &
    &fpx(nx), fpy(ny)
  integer :: nrx(mx), nry(my)
  !  ..local scalars..
  real(8) :: arg, cos, fac, pinv, piv, sin, term, one, half
  integer :: i, ibandx, ibandy, ic, iq, irot, it, iz, i1, i2, i3, j, k, k1, k2, l
  integer :: l1, l2, ncof, nk1x, nk1y, nrold, nroldx, nroldy, number, numx, numx1, numy, numy1, n1
  !  ..local arrays..
  real(8) :: h(7)
  !  ..subroutine references..
  !    fpback,fpbspl,fpgivs,fpdisc,fprota
  !  ..
  !  the b-spline coefficients of the smoothing spline are calculated as
  !  the least-squares solution of the over-determined linear system of
  !  equations  (ay) c (ax)' = q       where
  !
  !               |   (spx)    |            |   (spy)    |
  !        (ax) = | ---------- |     (ay) = | ---------- |
  !               | (1/p) (bx) |            | (1/p) (by) |
  !
  !                                | z  ' 0 |
  !                            q = | ------ |
  !                                | 0  ' 0 |
  !
  !  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
  !                b-spline coefficients.
  !       z      : the my x mx matrix which contains the function values.
  !       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
  !                matrices according to the least-squares problems in
  !                the x- and y-direction.
  !       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
  !                matrices which contain the discontinuity jumps of the
  !                derivatives of the b-splines in the x- and y-direction.
  one = 1._8
  half = 0.5_8
  nk1x = nx - kx1
  nk1y = ny - ky1
  if (p > 0._8) pinv = one / p
  !  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
  !  the value of p whether the matrices (spx),(spy),(bx) and (by) still
  !  must be determined.
  if (ifsx /= 0) go to 50
  !  calculate the non-zero elements of the matrix (spx) which is the
  !  observation matrix according to the least-squares spline approximat-
  !  ion problem in the x-direction.
  l = kx1
  l1 = kx2
  number = 0
  do it = 1, mx
    arg = x(it)
    do while (arg >= tx(l1) .and. l /= nk1x)
      l = l1
      l1 = l + 1
      number = number + 1
    end do

    call fpbspl(tx, nx, kx, arg, l, h)

    spx(it, :kx1) = h(:kx1)
    nrx(it) = number
  end do

  ifsx = 1
50 if (ifsy /= 0) go to 100
  !  calculate the non-zero elements of the matrix (spy) which is the
  !  observation matrix according to the least-squares spline approximat-
  !  ion problem in the y-direction.
  l = ky1
  l1 = ky2
  number = 0
  do it = 1, my
    arg = y(it)
    do while (arg >= ty(l1) .and. l /= nk1y)
      l = l1
      l1 = l + 1
      number = number + 1
    end do

    call fpbspl(ty, ny, ky, arg, l, h)
    spy(it, :ky1) = h(:ky1)
    nry(it) = number
  end do

  ifsy = 1
100 if (p <= 0._8) go to 120
  !  calculate the non-zero elements of the matrix (bx).
  if (ifbx /= 0 .or. nx == 2 * kx1) go to 110
  call fpdisc(tx, nx, kx2, bx, nx)
  ifbx = 1
  !  calculate the non-zero elements of the matrix (by).
110 if (ifby /= 0 .or. ny == 2 * ky1) go to 120
  call fpdisc(ty, ny, ky2, by, ny)
  ifby = 1
  !  reduce the matrix (ax) to upper triangular form (rx) using givens
  !  rotations. apply the same transformations to the rows of matrix q
  !  to obtain the my x (nx-kx-1) matrix g.
  !  store matrix (rx) into (ax) and g into q.
120 l = my * nk1x
  !  initialization.
  q(:l) = 0._8

  ax(:nk1x, kx2) = 0._8

  l = 0
  nrold = 0
  !  ibandx denotes the bandwidth of the matrices (ax) and (rx).
  ibandx = kx1
  do it = 1, mx
    number = nrx(it)
150 if (nrold == number) go to 180
    if (p <= 0._8) go to 260
    ibandx = kx2
    !  fetch a new row of matrix (bx).
    n1 = nrold + 1
    h(:kx2) = bx(n1, :kx2) * pinv
    !  find the appropriate column of q.
    right(:my) = 0.
    irot = nrold
    go to 210
    !  fetch a new row of matrix (spx).
180 h(ibandx) = 0.
    h(:kx1) = spx(it, :kx1)
    !  find the appropriate column of q.
    do j = 1, my
      l = l + 1
      right(j) = z(l)
    end do

    irot = number
    !  rotate the new row of matrix (ax) into triangle.
210 do i = 1, ibandx
      irot = irot + 1
      piv = h(i)
      if (piv == 0._8) cycle
      !  calculate the parameters of the givens transformation.
      call fpgivs(piv, ax(irot, 1), cos, sin)
      !  apply that transformation to the rows of matrix q.
      iq = (irot - 1) * my
      do j = 1, my
        iq = iq + 1
        call fprota(cos, sin, right(j), q(iq))
      end do

      !  apply that transformation to the columns of (ax).
      if (i == ibandx) go to 250
      i2 = 1
      i3 = i + 1
      do j = i3, ibandx
        i2 = i2 + 1
        call fprota(cos, sin, h(j), ax(irot, i2))
      end do

    end do

250 if (nrold == number) cycle
260 nrold = nrold + 1
    go to 150
  end do

  !  reduce the matrix (ay) to upper triangular form (ry) using givens
  !  rotations. apply the same transformations to the columns of matrix g
  !  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
  !  store matrix (ry) into (ay) and h into c.
  ncof = nk1x * nk1y
  !  initialization.
  c(:ncof) = 0._8
  ay(:nk1y, :ky2) = 0._8
  nrold = 0
  !  ibandy denotes the bandwidth of the matrices (ay) and (ry).
  ibandy = ky1
  do it = 1, my
    number = nry(it)
300 if (nrold == number) go to 330
    if (p <= 0.) go to 410
    ibandy = ky2
    !  fetch a new row of matrix (by).
    n1 = nrold + 1
    h(:ky2) = by(n1, :ky2) * pinv
    !  find the appropriate row of g.
    right(:nk1x) = 0._8
    irot = nrold
    go to 360
    !  fetch a new row of matrix (spy)
330 h(ibandy) = 0._8
    h(:ky1) = spy(it, :ky1)
    !  find the appropriate row of g.
    l = it
    do j = 1, nk1x
      right(j) = q(l)
      l = l + my
    end do

    irot = number
    !  rotate the new row of matrix (ay) into triangle.
360 do i = 1, ibandy
      irot = irot + 1
      piv = h(i)
      if (piv == 0.) cycle
      !  calculate the parameters of the givens transformation.
      call fpgivs(piv, ay(irot, 1), cos, sin)
      !  apply that transformation to the columns of matrix g.
      ic = irot
      do j = 1, nk1x
        call fprota(cos, sin, right(j), c(ic))
        ic = ic + nk1y
      end do

      !  apply that transformation to the columns of matrix (ay).
      if (i == ibandy) exit
      i2 = 1
      i3 = i + 1
      do j = i3, ibandy
        i2 = i2 + 1
        call fprota(cos, sin, h(j), ay(irot, i2))
      end do

    end do

    if (nrold == number) cycle
410 nrold = nrold + 1
    go to 300
  end do

  !  backward substitution to obtain the b-spline coefficients as the
  !  solution of the linear system    (ry) c (rx)' = h.
  !  first step: solve the system  (ry) (c1) = h.
  k = 1
  do i = 1, nk1x
    call fpback(ay, c(k), nk1y, ibandy, c(k), ny)
    k = k + nk1y
  end do

  !  second step: solve the system  c (rx)' = (c1).
  k = 0
  do j = 1, nk1y
    k = k + 1
    l = k
    do i = 1, nk1x
      right(i) = c(l)
      l = l + nk1y
    end do

    call fpback(ax, right, nk1x, ibandx, right, nx)
    l = k
    do i = 1, nk1x
      c(l) = right(i)
      l = l + nk1y
    end do
  end do

  !  calculate the quantities
  !    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
  !    fp = sumi=1,mx(sumj=1,my(res(i,j)))
  !    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
  !                  tx(r+kx) <= x(i) <= tx(r+kx+1)
  !    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
  !                  ty(r+ky) <= y(j) <= ty(r+ky+1)
  fp = 0._8
  fpx(:nx) = 0._8
  fpy(:ny) = 0._8
  nk1y = ny - ky1
  iz = 0
  nroldx = 0
  !  main loop for the different grid points.
  do i1 = 1, mx
    numx = nrx(i1)
    numx1 = numx + 1
    nroldy = 0
    do i2 = 1, my
      numy = nry(i2)
      numy1 = numy + 1
      iz = iz + 1
      !  evaluate s(x,y) at the current grid point by making the sum of the
      !  cross products of the non-zero b-splines at (x,y), multiplied with
      !  the appropriate b-spline coefficients.
      term = 0.
      k1 = numx * nk1y + numy
      do l1 = 1, kx1
        k2 = k1
        fac = spx(i1, l1)
        do l2 = 1, ky1
          k2 = k2 + 1
          term = term + fac * spy(i2, l2) * c(k2)
        end do

        k1 = k1 + nk1y
      end do

      !  calculate the squared residual at the current grid point.
      term = (z(iz) - term)**2
      !  adjust the different parameters.
      fp = fp + term
      fpx(numx1) = fpx(numx1) + term
      fpy(numy1) = fpy(numy1) + term
      fac = term * half
      if (numy == nroldy) go to 530
      fpy(numy1) = fpy(numy1) - fac
      fpy(numy) = fpy(numy) + fac
530   nroldy = numy
      if (numx == nroldx) cycle
      fpx(numx1) = fpx(numx1) - fac
      fpx(numx) = fpx(numx) + fac
    end do

    nroldx = numx
  end do

end subroutine fpgrre

