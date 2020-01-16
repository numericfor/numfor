program test_rand_gauss
  use numfor, only: dp, Zero, str, random_seed, random_normal
  implicit none
  real(dp), dimension(2, 3) :: samples
  call random_seed(0)           ! Uses default seed
  call random_normal(loc=1._dp, scale=2._dp, x=samples)
  print "(A)", str(samples)
  ! Outputs:
  ! [[ 2.579691898234, -0.3742516980564]
  !  [ 1.1897226266753, 1.4022523097265]
  !  [ 0.4435270050626, 0.8899184474152]]
end program test_rand_gauss
