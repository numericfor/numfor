program test
  use utils
  implicit none

  integer, parameter :: N = 100000
  real(dp), dimension(N) :: a, b, c
  integer :: i, j
  type(timer) T1

  ! Test timer
  call T1%start

  a = 1._dp
  c = 0
  do j = 1, 2
    do i = 1, N
      b(i) = (-1)**i * a(i)
    end do
    do i = 1, N
      c = c + b(i)
    end do
    c = -c
  end do

  call T1%finish()

  call T1%show()

end program test

