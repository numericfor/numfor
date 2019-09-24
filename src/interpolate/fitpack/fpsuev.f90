subroutine fpsuev(idim, tu, nu, tv, nv, c, u, mu, v, mv, f, wu, wv, lu, lv)
  !  ..scalar arguments..
  integer :: idim, nu, nv, mu, mv
  !  ..array arguments..
  integer :: lu(mu), lv(mv)
  real(8) :: tu(nu), tv(nv), c((nu - 4) * (nv - 4) * idim), u(mu), v(mv),&
    &f(mu * mv * idim), wu(mu, 4), wv(mv, 4)
  !  ..local scalars..
  integer :: i, i1, j, j1, k, l, l1, l2, l3, m, nuv, nu4, nv4
  real(8) :: arg, sp, tb, te
  !  ..local arrays..
  real(8) :: h(4)
  !  ..subroutine references..
  !    fpbspl
  !  ..
  nu4 = nu - 4
  tb = tu(4)
  te = tu(nu4 + 1)
  l = 4
  l1 = l + 1
  do i = 1, mu
    arg = u(i)
    if (arg < tb) arg = tb
    if (arg > te) arg = te
10  if (arg < tu(l1) .or. l == nu4) go to 20
    l = l1
    l1 = l + 1
    go to 10
20  call fpbspl(tu, nu, 3, arg, l, h)
    lu(i) = l - 4
    wu(i, :4) = h(:4)
  end do

  nv4 = nv - 4
  tb = tv(4)
  te = tv(nv4 + 1)
  l = 4
  l1 = l + 1
  do i = 1, mv
    arg = v(i)
    if (arg < tb) arg = tb
    if (arg > te) arg = te
50  if (arg < tv(l1) .or. l == nv4) go to 60
    l = l1
    l1 = l + 1
    go to 50
60  call fpbspl(tv, nv, 3, arg, l, h)
    lv(i) = l - 4
    wv(i, :4) = h(:4)
  end do

  m = 0
  nuv = nu4 * nv4
  do k = 1, idim
    l3 = (k - 1) * nuv
    do i = 1, mu
      l = lu(i) * nv4 + l3

      h(:4) = wu(i, :4)

      do j = 1, mv
        l1 = l + lv(j)
        sp = 0.
        do i1 = 1, 4
          l2 = l1
          do j1 = 1, 4
            l2 = l2 + 1
            sp = sp + c(l2) * h(i1) * wv(j, j1)
          end do

          l1 = l1 + nv4
        end do

        m = m + 1
        f(m) = sp
      end do
    end do
  end do

end subroutine fpsuev
