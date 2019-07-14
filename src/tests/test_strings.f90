program test
  use strings
  implicit none

  character(len=:), allocatable :: c1, c2
  integer :: i1
  logical :: l1

  c1 = "NumFor - Library for Simple Numerical computing."

  i1 = count_sub(c1, "Num")
  if (i1 /= 2) print *, '-- Failed count (2) --'
  i1 = count_sub(c1, "num")
  if (i1 /= 0) stop '-- Failed count (0) --'

  print *, 'Original: ', c1
  c2 = replace(c1, "Num", "Sci")
  print *, 'Replacing (Num -> Sci): ', c2
  c2 = replace(c2, "Sci", "Num")
  if (c1 /= c2) stop 'Failed replace'
  print *, 'UPPER:', upper(c1)
  print *, 'LOWER:', lower(c1)

  c1 = "-One"
  print *, "MYSTRING = ", c1
  c2 = zfill(c1, 10)
  print *, '| C1 | ZFILL(c1, 10) | = |', c1, '|', c2, '|'

  c1 = "  -One"
  c2 = zfill(c1, 10)
  print *, '| C1 | ZFILL(c1, 10) | = |', c1, '|', c2, '|'

  c1 = "Three"
  !print *,"MYSTRING = ",c1
  c2 = center(c1, 2)
  c2 = center(c1, 5, fillchar="x")

  c2 = center(c1, 7, fillchar="x")
  print *, "CENTER = ", c2
  c2 = center(c1, 20, fillchar="x")
  print *, "CENTER = ", ":", c2, ":"
  c2 = center(c1, 20)
  print *, "CENTER = ", ":", c2, ":"

  c1 = "NumFor - Library for Simple Numerical computing. - NumFor"

  c2 = reverse(c1)
  IF (c1 /= reverse(c2)) stop "Reverse failed"

  l1 = issub(c1, 'for') ! TRUE
  IF (.not. l1) stop 'failed issub 1'

  l1 = issub(c1, "Vision") ! FALSE
  IF (l1) stop 'failed issub 2'

  l1 = issub(c1, "FOR") ! False
  IF (l1) stop 'failed issub 3'

  l1 = issub(upper(c1), "FOR") ! True
  IF (.not. l1) stop 'failed issub 4'

  c2 = '  hola    '
  print *, '|', c2, '|'
  print *, '|', adjustl(c2), '|'

  IF (swapcase(swapcase(c1)) .ne. c1) stop 'swapcase failed'

  ! test strip
  print *, 'Original:', c1, ' --- Chars to remove:', "Num"
  print *, 'Stripped:', strip(c1, "Num")

  ! test strip
  print *, 'Original:', c1, ' --- Chars to remove:', "NumFor"
  print *, 'Stripped:', strip(c1, "NumFor")

  print *, 'I2str: - int:', 123, '- string: |'//str(123)//'|'
  print *, 'r2str:  - real:', 123.121000_4, '- string: |'//str(123.121000_4)//'|'
  print *, 'dp2str: - real:', 123.1221_8, '- string: |'//str(123.1221_8)//'|'
  print *, 'dp2str: - real:', -123.121_8, '- string: |'//str(-123.121_8)//'|'
  print *, 'dp2str: - real:', -123.1221e-12_8, '- string: |'//str(-123.1221e-12_8)//'|'
  print *, '|', '  12312 ', '|'
  print *, 'rjust: |', rjust('  12312 ', 20), '|'
  print *, 'ljust: |', ljust('  12312 ', 20), '|'

  c1 = achar(9)//'tab-indented with 4 trailing spaces    '
  print *, '|'//c1//'|'//rstrip(c1)//'|'
  print *, '|'//c1//'|'//lstrip(c1)//'|'
  print *, '|'//c1//'|'//strip(c1)//'|'
  print *, '|'//c1//'|'//trim(adjustl(c1))//'|'
end program test
