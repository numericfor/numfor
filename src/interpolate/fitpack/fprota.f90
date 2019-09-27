subroutine fprota(cos, sin, a, b)
  !  subroutine fprota applies a givens rotation to a and b.
  !  ..
  !  ..scalar arguments..
  real(8) :: cos, sin, a, b
  ! ..local scalars..
  real(8) :: stor1, stor2
  !  ..
  stor1 = a
  stor2 = b
  b = cos * stor2 + sin * stor1
  a = cos * stor1 - sin * stor2
end subroutine fprota
