subroutine fpbisp(tx, nx, ty, ny, c, kx, ky, x, mx, y, my, z, wx, wy, lx, ly)
  !  ..scalar arguments..
  integer :: nx, ny, kx, ky, mx, my
  !  ..array arguments..
  integer :: lx(mx), ly(my)
  real(8) :: tx(nx), ty(ny), c((nx - kx - 1) * (ny - ky - 1)), x(mx), y(my), z(mx * my), &
    &wx(mx, kx + 1), wy(my, ky + 1)
  !  ..local scalars..
  integer :: kx1, ky1, l, l1, l2, m, nkx1, nky1
  real(8) :: arg, sp, tb, te
  !  ..local arrays..
  real(8) :: h(6)
  !  ..subroutine references..
  !    fpbspl
  !  ..
  kx1 = kx + 1
  nkx1 = nx - kx1
  tb = tx(kx1)
  te = tx(nkx1 + 1)
  l = kx1
  l1 = l + 1
  do i = 1, mx
    arg = x(i)
    if (arg < tb) arg = tb
    if (arg > te) arg = te
10  if (arg < tx(l1) .or. l == nkx1) go to 20
    l = l1
    l1 = l + 1
    go to 10
20  call fpbspl(tx, nx, kx, arg, l, h)
    lx(i) = l - kx1
    do j = 1, kx1
      wx(i, j) = h(j)
    end do

  end do

  ky1 = ky + 1
  nky1 = ny - ky1
  tb = ty(ky1)
  te = ty(nky1 + 1)
  l = ky1
  l1 = l + 1
  do i = 1, my
    arg = y(i)
    if (arg < tb) arg = tb
    if (arg > te) arg = te
50  if (arg < ty(l1) .or. l == nky1) go to 60
    l = l1
    l1 = l + 1
    go to 50
60  call fpbspl(ty, ny, ky, arg, l, h)
    ly(i) = l - ky1
    do j = 1, ky1
      wy(i, j) = h(j)
    end do

  end do

  m = 0
  do i = 1, mx
    l = lx(i) * nky1
    do i1 = 1, kx1
      h(i1) = wx(i, i1)
    end do

    do j = 1, my
      l1 = l + ly(j)
      sp = 0.
      do i1 = 1, kx1
        l2 = l1
        do j1 = 1, ky1
          l2 = l2 + 1
          sp = sp + c(l2) * h(i1) * wy(j, j1)
        end do

        l1 = l1 + nky1
      end do

      m = m + 1
      z(m) = sp
    end do

  end do

end subroutine fpbisp

