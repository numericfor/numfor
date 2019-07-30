!> test
program test
  USE basic, only: dp, timer
  USE strings, only: center
  USE grids, only: linspace, arange
  USE sorting, only: searchsorted, sort
  implicit none
  integer, parameter :: Ndat = 12 * 100000
  real(dp), dimension(Ndat) :: a
  real(dp), dimension(Ndat) :: b
  real(dp) :: elem
  type(timer) :: T1, T2, T
  integer :: n, i
  real(dp), dimension(:), allocatable :: s
  ! integer, dimension(:), allocatable :: ns

  print *, real(huge(1), kind=dp) / Ndat
  call random_number(elem)
  elem = elem * 10

  print *, 'elem', elem
  T1 = 0
  T2 = 0
  do i = 1, 3

    call random_number(a)
    a = a * 10_dp - 5._dp
    b = a
    ! print "(10(g0.4,1x))", a(:20)

    ! print "(A)", 'Sorting in descending order'
    call T%start()
    call sort(a)
    call T%stop()
    T1 = T1 + T

    call T%start()
    call sort(b, reverse=.True.)
    call T%stop()
    T2 = T2 + T

  end do

  print "(10(g0.8, 1x))", a(:7)
  print "(10(g0.8, 1x))", a(Ndat - 6:)
  print "(10(g0.8, 1x))", b(:7)
  print "(10(g0.8, 1x))", b(Ndat - 6:)

  print "(A)", center('Total times', 65, '-')
  print "(A)", repeat('*', 60)
  print "(A)", center('Descending algorithm', 65, '-')
  print "(4(I4.4,1x))", T1%date_elapsed(4:7)
  call T1%show()
  print "(A)", center('Ascending algorithm', 65, '-')
  print "(4(I4.4,1x))", T2%date_elapsed(4:7)
  call T2%show()
  print "(A)", repeat('*', 60)

  T1 = T1 / 3
  T2 = T2 / 3

  print "(A)", center('Average times', 65, '-')
  print "(A)", repeat('*', 60)
  print "(A)", center('Descending algorithm', 65, '-')
  print "(4(I4.4,1x))", T1%date_elapsed(4:7)
  call T1%show()
  print "(A)", center('Ascending algorithm', 65, '-')
  print "(4(I4.4,1x))", T2%date_elapsed(4:7)
  call T2%show()
  print "(A)", repeat('*', 60)

  ! call T1%show()
  ! call T2%show()

  call sort(a)
  ! print "(10(g0.4,1x))", a(:20)

! print "(A)", 'Sorted in ascending order'
! print "(10(g0.4,1x))", a

  n = searchsorted(a, elem)
  print *, ''
  ! print "(A,f7.3,1x,A,I8.8)", 'put: ', elem, 'in position after n = ', n

  s = linspace(1._dp, 9.3_dp, 8)

! print "(A)", "Array:"
! print "(6(f8.4,1x))", a

  do i = 1, size(s)
    n = searchsorted(a, s(i))
    print "(A,f9.7,A,I8.8,(2(A,f9.7)))", "Put ", s(i), " in n=", n, " between ", a(n), ' and ', a(n + 1)
  end do

! ns = arange(0, 13, 2)

! print "(A,7(I2.2,1x))", "Array: ", ns
! print "(A)", "  Find   in"
! do i = 1, size(s)
!   print "(f8.4,1x,I2.2)", s(i), searchsorted(ns, s(i))
! end do

end program test
