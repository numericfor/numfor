program test_rand_gauss
  use numfor, only: dp, str, random_seed, random_normal
  implicit none
  integer :: i
  real(dp) :: rn
  integer, parameter :: Ndim = 10**0
  ! integer, dimension(6) :: hist

  ! hist = 0
  call random_seed()     ! Uses a random seed
  do i = 1, Ndim
    rn = random_gauss()
    ! hist(rn) = hist(rn) + 1
  end do
  print "(A)", 'hist = '//str(hist)

  ! Outputs:

end program test_rand_gauss
