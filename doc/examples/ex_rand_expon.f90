program test_rand_expon
  use numfor, only: dp, Zero, str, random_seed, random_exponential
  implicit none
  real(dp), dimension(2, 3) :: samples
  call random_seed(0)           ! Uses default seed
  call random_exponential(2._dp, samples)
  print "(A)", str(samples)
  ! Outputs:
  ! [[ 3.0912457577868, 0.5766454637618]
  !  [ 2.4803832429473, 5.8624300413216]
  !  [ 0.0389183312912, 1.0380588489303]]
end program test_rand_expon
