program ex_random
  use numfor, only: i8, random_seed, random_real
  implicit none
  call random_seed(0_i8)     ! This uses the default seed
  print *, random_real()
  ! Outputs:
  !    0.78682095486780212
end program ex_random

