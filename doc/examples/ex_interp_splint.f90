program ex_spl_int
  USE numfor, only: dp, Zero, M_PI, linspace
  USE numfor, only: UnivSpline, splrep, splint
  implicit none
  integer, parameter :: N = 6
  real(dp), dimension(N) :: x
  real(dp), dimension(N) :: y
  real(dp) :: yI

  type(UnivSpline) :: tck

  x = linspace(Zero, M_PI, N)
  y = sin(x)
  !< [integrate]
  call splrep(x, y, tck=tck, s=0._dp)
  call splint(tck, Zero, M_PI / 2._dp, yI) ! yI = 1.000630799770
  !< [integrate]
  print "(A, f14.12)", 'Integral entre 0 y π/2: ', yI
  ! Integral entre 0 y π/2: 1.000630799770
  call splint(tck, Zero, M_PI, yI)
  print "(A, f14.12)", 'Integral entre 0 y π: ', yI
  ! Integral entre 0 y π: 2.001261599540
end program ex_spl_int

