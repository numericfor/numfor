program ex_rand_int
  use numfor, only: dp, i8, str, random_seed, random_int
  implicit none
  integer :: i
  integer(i8) :: rn
  integer(i8) :: ymin, ymax, range

  ymin = ishft(-2_i8**62, 1)
  ymax = 2_i8**62 + 1
  range = ymax / 2 - ymin / 2
  print *, ymin, ymax, ishft(ymax, 1)
  print *, ymax / 10, ymin / 10, ymax / 10 - ymin / 10

  call random_seed()     ! Uses a random seed
  do i = 1, 10
    rn = random_int(6_i8) + 1
    print *, rn
  end do
  print *, 'size', bit_size(rn)
  ! Outputs:

end program ex_rand_int
