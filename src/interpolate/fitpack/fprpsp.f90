subroutine fprpsp(nt, np, co, si, c, f, ncoff)
  !  given the coefficients of a spherical spline function, subroutine
  !  fprpsp calculates the coefficients in the standard b-spline re-
  !  presentation of this bicubic spline.
  !  ..
  !  ..scalar arguments
  integer :: nt, np, ncoff
  !  ..array arguments
  real(8) :: co(np), si(np), c(ncoff), f(ncoff)
  !  ..local scalars
  real(8) :: cn, c1, c2, c3
  integer :: i, ii, j, k, l, ncof, npp, np4, nt4
  !  ..
  nt4 = nt - 4
  np4 = np - 4
  npp = np4 - 3
  ncof = 6 + npp * (nt4 - 4)
  c1 = c(1)
  cn = c(ncof)
  j = ncoff
  do i = 1, np4
    f(i) = c1
    f(j) = cn
    j = j - 1
  end do

  i = np4
  j = 1
  do l = 3, nt4
    ii = i
    if (l == 3 .or. l == nt4) go to 30
    do k = 1, npp
      i = i + 1
      j = j + 1
      f(i) = c(j)
    end do

    go to 50
30  if (l == nt4) c1 = cn
    c2 = c(j + 1)
    c3 = c(j + 2)
    j = j + 2
    do k = 1, npp
      i = i + 1
      f(i) = c1 + c2 * co(k) + c3 * si(k)
    end do

50  do k = 1, 3
      ii = ii + 1
      i = i + 1
      f(i) = f(ii)
    end do

  end do

  c(:ncoff) = f(:ncoff)

end subroutine fprpsp
