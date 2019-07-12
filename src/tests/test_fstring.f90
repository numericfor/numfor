program test
  USE fString

  implicit none

  type(fStr) :: mystr1
  type(fStr) :: mystr2
  type(fStr) :: mystr3
  character(len=:), allocatable :: c1
  character(len=10) :: iomsg = ' '
  integer :: i1
  ! logical :: l1

  c1 = "NumFor - Library for Simple Numerical computing."
  mystr1 = c1
  mystr2 = fStr(123123.e6)
  mystr3 = fStr("hola")
  print *, mystr1
  print *, mystr2
  print *, mystr3
  write (6, *) mystr3//mystr2
  write (6, *) mystr3//'chau'
  write (6, '(A)') mystr3       ! Esto no funciona
  call mystr3%writef(6, '(A)', [2], i1, iomsg)
  print *, mystr1%startswith('Num')
  print *, mystr1%startswith('Nume')
  print *, reverse(mystr1)
  print *, mystr1%reverse()
  print *, mystr1 == mystr2
  print *, 3 * fStr("hola") * 3
  print *, len(mystr1)

  print *, mystr1%replace('Num', 'Sci')
  ! call writef(mystr2)

  print *, 'i appears ', mystr1%count('i'), 'times'
end program test
