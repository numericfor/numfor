subroutine fpsysy(a, n, g)
  ! subroutine fpsysy solves a linear n x n symmetric system
  !    (a) * (b) = (g)
  ! on input, vector g contains the right hand side ; on output it will
  ! contain the solution (b).
  !  ..
  !  ..scalar arguments..
  integer :: n
  !  ..array arguments..
  real(8) :: a(6, 6), g(6)
  !  ..local scalars..
  real(8) :: fac
  integer :: i, i1, j, k
  !  ..
  g(1) = g(1) / a(1, 1)
  if (n == 1) return
  !  decomposition of the symmetric matrix (a) = (l) * (d) *(l)'
  !  with (l) a unit lower triangular matrix and (d) a diagonal
  !  matrix

  a(2:n, 1) = a(2:n, 1) / a(1, 1)
  do i = 2, n
    i1 = i - 1
    do k = i, n
      fac = a(k, i)
      do j = 1, i1
        fac = fac - a(j, j) * a(k, j) * a(i, j)
      end do
      a(k, i) = fac
      if (k > i) a(k, i) = fac / a(i, i)
    end do
  end do

  !  solve the system (l)*(d)*(l)'*(b) = (g).
  !  first step : solve (l)*(d)*(c) = (g).
  do i = 2, n
    i1 = i - 1
    fac = g(i)
    do j = 1, i1
      fac = fac - g(j) * a(j, j) * a(i, j)
    end do

    g(i) = fac / a(i, i)
  end do

  !  second step : solve (l)'*(b) = (c)
  i = n
  do j = 2, n
    i1 = i
    i = i - 1
    fac = g(i)
    do k = i1, n
      fac = fac - g(k) * a(k, i)
    end do

    g(i) = fac
  end do

end subroutine fpsysy
