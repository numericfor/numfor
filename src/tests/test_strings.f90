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
  !print *,"MYSTRING = ",c1
  !print *,"MYSTRING = ",c1
  c2 = upper(c1)
  !print *,"UPPER = ",upper(c1)
  c2 = lower(c1)
  !print *,"LOWER = ",lower(c1)

  ! c1="Ironman,Thor,Thanos,Black Panther,Winter Soldier"
  ! !print *,"MYSTRING = ",c1
  ! c2=split(c1,",",3)
  ! !print *,"SPLIT = ",c2
  ! c2=split(c1,",",4)
  ! !print *,"SPLIT = ",c2
  ! c2=split(c1,",",0)
  ! !print *,"SPLIT = ",c2
  ! c2=split(c1,",",3)
  ! !print *,"SPLIT = ",c2

  ! c1=",Ironman,Thor,,Thanos,Black Panther,Thanos,Thanos,Thor,Thor,,Winter Soldier,,,"
  ! !print *,"MYSTRING = ",c1
  ! c2=uniq(c1,",")
  ! !print *,"UNIQ = ",c2

  c1 = "Thor"
  print *, "MYSTRING = ", c1
  c2 = zfill(c1, 10)
  print *, "ZFILL = ", c2
  c2 = zfill(c1, 2)
  print *, "ZFILL = ", c2

  c1 = "Thor"
  !print *,"MYSTRING = ",c1
  c2 = center(c1, 2)
  !print *,"CENTER = ",c2
  c2 = center(c1, 5, fillchar="x")
  !print *,"CENTER = ",c2
  c2 = center(c1, 7, fillchar="x")
  print *, "CENTER = ", c2
  c2 = center(c1, 20, fillchar="x")
  print *, "CENTER = ", ":", c2, ":"
  c2 = center(c1, 20)
  print *, "CENTER = ", ":", c2, ":"

  c1 = "Ironman,Thor,Thanos,Black Panther,Winter Soldier"
  !print *,"MYSTRING = ",c1
  c2 = reverse(c1)
  !print *,"REVERSE = ",":",c2,":"

  c1 = "Ironman,Thor,Thanos,Black Panther,Winter Soldier"
  l1 = issub(c1, "Thor") ! TRUE
  ! print *, '1.- l1', l1
  if (.not. l1) stop 'failed issub 1'
  l1 = issub(c1, "Vision") ! FALSE
  ! print *, '2.- l1', l1
  if (l1) stop 'failed issub 2'
  l1 = issub(c1, "THOR") ! False
  ! print *, '3.- l1', l1
  if (l1) stop 'failed issub 3'
  l1 = issub(upper(c1), "THOR") ! True
  ! print *, '4.- l1', l1
  if (.not. l1) stop 'failed issub 4'

  c1 = "ironman,Thor,Thanos,Black Panther,Winter Soldier"
  if (swapcase(swapcase(c1)) .ne. c1) stop 'swapcase failed'

  ! test strip
  print *, 'Original:', c1, ' --- Chars to remove:', "ire"
  print *, 'Stripped:', strip(c1, "ire")

  print *, 'I2str: - int:', 123, '- string: |'//str(123)//'|'
  print *, 'r2str:  - real:', 123.121000_4, '- string: |'//str(123.121000_4)//'|'
  print *, 'dp2str: - real:', 123.1221_8, '- string: |'//str(123.1221_8)//'|'
  print *, 'dp2str: - real:', -123.121_8, '- string: |'//str(-123.121_8)//'|'
  print *, 'dp2str: - real:', -123.1221e-12_8, '- string: |'//str(-123.1221e-12_8)//'|'
end program test
