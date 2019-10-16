program ex_timer
  USE numfor, only: dp, timer
  implicit none

  integer, parameter :: N = 100000
  real(dp), dimension(N) :: a, b, c
  integer :: i, j
  type(timer) T1

  call T1%start()               ! Start the timer
  ! ---- Do your operations ....
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
  ! -----------------------------
  call T1%stop()                ! Stop the timer
  call T1%show()                ! Show the results
end program ex_timer

