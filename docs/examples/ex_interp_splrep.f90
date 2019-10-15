program ex_splprep
  USE numfor, only: dp, Zero, M_PI, str, linspace, save_array
  USE numfor, only: UnivSpline, splrep, splev
  implicit none

  integer, parameter :: N = 6
  integer, parameter :: Nnew = 29
  real(dp), dimension(N) :: x
  real(dp), dimension(N) :: y
  real(dp), dimension(Nnew) :: xnew
  real(dp), dimension(Nnew) :: ynew
  character(len=:), allocatable :: header
  real(dp) :: s
  type(UnivSpline) :: tck

  ! Create data
  s = 0._dp
  x = linspace(Zero, M_PI, N)
  y = sin(x)
  xnew = linspace(Zero, M_PI, Nnew)

  !< [Using it]
  ! After setting data in x
  ! interpolate: (step one) Create interpolation and parameters u
  call splrep(x, y, tck=tck, s=s)
  ! interpolate: (step two) Apply to points u
  call splev(xnew, tck, ynew)
  !< [Using it]

  ! save to file
  header = "x                   y"
  call save_array([xnew, ynew], 2, fmt="f17.14", header=header)
end program ex_splprep
! Prints:
!
! 0.00000000000000  0.00000000000000
! 0.11219973762821  0.11401345780438
! 0.22439947525641  0.22522911715938
! 0.33659921288462  0.33269883705930
! 0.44879895051283  0.43547447649845
! 0.56099868814103  0.53260789447113
!  ...
!  ...
! 2.46839422782055  0.62315094997166
! 2.58059396544876  0.53260789447113
! 2.69279370307697  0.43547447649845
! 2.80499344070517  0.33269883705930
! 2.91719317833338  0.22522911715938
! 3.02939291596159  0.11401345780438
! 3.14159265358979  0.00000000000000
!
