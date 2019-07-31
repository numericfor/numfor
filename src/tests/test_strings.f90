program test
  use basic
  use strings
  implicit none

  character(len=:), allocatable :: st1, st2
  integer :: i1
  logical :: l1
  complex(dp) x, y, z
  x = (1._dp, +1.e-50_dp)
  y = (1._dp, -1.e-60_dp)
  z = x / y
  print *, z

  st1 = "NumFor - Library for Simple Numerical computing."

  i1 = count_sub(st1, "Num")
  if (i1 /= 2) print *, '-- Failed count (2) --'
  i1 = count_sub(st1, "num")
  if (i1 /= 0) stop '-- Failed count (0) --'

  print *, 'Original: ', st1
  st2 = replace(st1, "Num", "Sci")
  print *, 'Replacing (Num -> Sci): ', st2
  st2 = replace(st2, "Sci", "Num")
  if (st1 /= st2) stop 'Failed replace'
  print *, 'UPPER:', upper(st1)
  print *, 'LOWER:', lower(st1)

  st1 = "-One"
  print *, "MYSTRING = ", st1
  st2 = zfill(st1, 10)
  print *, '| ST1 | ZFILL(st1, 10) | = |', st1, '|', st2, '|'

  st1 = "  -One"
  st2 = zfill(st1, 10)
  print *, '| ST1 | ZFILL(st1, 10) | = |', st1, '|', st2, '|'

  st1 = "Three"
  !print *,"MYSTRING = ",st1
  st2 = center(st1, 2)
  st2 = center(st1, 5, fillchar="x")

  st2 = center(st1, 7, fillchar="x")
  print *, "CENTER = ", st2
  st2 = center(st1, 20, fillchar="x")
  print *, "CENTER = ", ":", st2, ":"
  st2 = center(st1, 20)
  print *, "CENTER = ", ":", st2, ":"

  st1 = "NumFor - Library for Simple Numerical computing. - NumFor"

  st2 = reverse(st1)
  IF (st1 /= reverse(st2)) stop "Reverse failed"
  print *, center(' Original ', 20, '-')
  print *, st1
  print *, center(' Reversed ', 20, '-')
  print *, st2

  l1 = issub(st1, 'for') ! TRUE
  IF (.not. l1) stop 'failed issub 1'

  l1 = issub(st1, "Vision") ! FALSE
  IF (l1) stop 'failed issub 2'

  l1 = issub(st1, "FOR") ! False
  IF (l1) stop 'failed issub 3'

  l1 = issub(upper(st1), "FOR") ! True
  IF (.not. l1) stop 'failed issub 4'

  st2 = '  hola    '
  print *, '|', st2, '|'
  print *, '|', adjustl(st2), '|'

  IF (swapcase(swapcase(st1)) .ne. st1) stop 'swapcase failed'

  ! test strip
  print *, 'Original:', st1, ' --- Chars to remove:', "Num"
  print *, 'Stripped:', strip(st1, "Num")

  ! test strip
  print *, 'Original:', st1, ' --- Chars to remove:', "NumFor"
  print *, 'Stripped:', strip(st1, "NumFor")

  print *, 'I2str: - int:', 123, '- string: |'//str(123)//'|'
  print *, 'r2str:  - real:', 123.1219000_sp, '- string: |'//str(123.1219000_sp)//'|'
  print *, 'dp2str: - real:', 123.1221_dp, '- string: |'//str(123.1221_dp)//'|'
  print *, 'dp2str: - real:', -123.121_dp, '- string: |'//str(-123.121_dp)//'|'
  print *, 'dp2str: - real:', -123.1221e-12_dp, '- string: |'//str(-123.1221e-12_dp)//'|'
  print *, 'dp2str: - real:', M_PI, '- string: |'//str(M_PI)//'|'
  print *, '|', '  12312 ', '|'
  print *, 'rjust: |', rjust('  12312 ', 20), '|'
  print *, 'ljust: |', ljust('  12312 ', 20), '|'

  st1 = achar(9)//'tab-indented with 4 trailing spaces    '
  print *, '|'//st1//'|'//rstrip(st1)//'|'
  print *, '|'//st1//'|'//lstrip(st1)//'|'
  print *, '|'//st1//'|'//strip(st1)//'|'
  print *, '|'//st1//'|'//trim(adjustl(st1))//'|'
end program test
