program ex_rand_int
  use numfor, only: dp, i8, str, random_seed, random_int
  implicit none
  integer :: i
  integer(i8) :: rn
  integer, parameter :: Ndim = 10**7
  integer, dimension(6) :: hist

  hist = 0
  call random_seed()     ! Uses a random seed
  do i = 1, Ndim
    rn = random_int(6_i8) + 1
    hist(rn) = hist(rn) + 1
  end do
  print "(A)", 'hist = '//str(hist)

  hist = 0
  do i = 1, Ndim
    rn = random_int(1_i8, 6_i8)
    hist(rn) = hist(rn) + 1
  end do
  print "(A)", 'hist = '//str(hist)
  ! Outputs:

end program ex_rand_int
