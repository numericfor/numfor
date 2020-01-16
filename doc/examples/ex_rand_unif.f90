program ex_rand_unif
  use numfor, only: dp, str, random_seed, random_sample, random_uniform
  implicit none
  !< [usesample]
  real(dp), dimension(2, 3) :: r1
  call random_seed(0)     ! Uses the default seed
  ! 2d-array of random numbers
  call random_sample(r1)
  print "(A)", str(r1)
  ! Outputs:
  ! [[ 0.7868209548678, 0.250480340688]
  !  [ 0.7106712289787, 0.946667800961]
  !  [ 0.0192710581958, 0.4049021448162]]
  !< [usesample]
  !< [useuniform]
  call random_uniform(-5._dp, 5._dp, r1)
  print "(A)", str(r1)
  ! Outputs:
  ! [[ -2.4868218207196, -4.7728756137207]
  !  [ 0.2064315257349, -1.5532969392081]
  !  [ -2.2580439639714, 0.6103210017639]]
  !< [useuniform]
end program ex_rand_unif
