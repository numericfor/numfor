!> test
program test
  USE basic, only: dp, timer
  USE strings, only: center, Str
  USE grids, only: linspace, arange
  USE sorting, only: searchsorted, sort
  implicit none
  integer, parameter :: Ndat = 1000 * 10000
  real(dp), dimension(Ndat) :: a
  real(dp), dimension(Ndat) :: b
  real(dp) :: elem
  type(timer) :: T1, T2, T
  integer :: n, i
  real(dp), dimension(:), allocatable :: s
  ! integer, dimension(:), allocatable :: ns

  ! print *, real(huge(1), kind=dp) / Ndat
  ! print *, Str(huge(1))
  call random_number(elem)
  elem = elem * 10

  print "(A)", repeat('*', 65)
  print "(A)", center(' Sort '//Str(Ndat)//' elements ', 65, '*')
  print "(A)", center(' Average times for '//Str(3)//' runs ', 65, '*')

  T1 = 0
  T2 = 0
  do i = 1, 3

    call random_number(a)
    a = a * 10_dp - 5._dp
    b = a

    ! Sorting in ascending order
    call T%start()
    call sort(a)
    call T%stop()
    T1 = T1 + T

    ! Sorting in descending order
    call T%start()
    call sort(b, reverse=.True.)
    call T%stop()
    T2 = T2 + T

  end do

  T1 = T1 / 3
  T2 = T2 / 3
  print "(A)", repeat('*', 65)
  print "(A)", center(' Ascending algorithm ', 65, '-')
  call T1%show()
  print "(A)", center(' Descending algorithm ', 65, '-')
  call T2%show()
  print "(A)", repeat('*', 65)

  ! Sort in ascending order
  call sort(a)
  n = searchsorted(a, elem)
  s = linspace(0.1_dp, 9.9_dp, 11) - 5._dp

  print "(A)", center(' Test searchsorted ', 65, '*')
  print "(A)", repeat('*', 65)
  do i = 1, size(s)
    n = searchsorted(a, s(i))
    print "(A,f0.7,A,I8.8,(2(A,f0.7)))", "Put ", s(i), " in n=", n, " between ", a(n), ' and ', a(n + 1)
  end do

end program test
