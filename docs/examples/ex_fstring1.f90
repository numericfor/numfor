program example_fstring1
  USE fString, only: fStr, len
  implicit none
  type(fStr) :: dir
  type(fStr) :: rfname
  type(fStr) :: output
  type(fStr) :: title
  real(8) :: param1
  real(8) :: param2
  integer :: param3
  integer, parameter :: Nmax = 5
  integer :: cnt
  integer :: width

  print *, 'Enter directory name'
  ! We can just set the output folder it
  ! dir = 'out'
  ! Or we could read the path from standard input.
  read (*, *) dir
  dir = dir//'/'                ! Add separator

  width = 40 + len(dir)         ! Used for centering
  title = ' These will the files created: '

  ! Set some parameters
  param1 = 1.3_8
  param2 = 0.07_8
  param3 = 42

  ! rfname holds the root of the filename
  rfname = "out"//fStr(param1)//'cm'//fStr(param3)//'eV_'
  print *, ''
  print *, dir%lower()
  print *, title%center(width, '=')
  do cnt = 1, Nmax
    ! Create full filename
    ! Here we have to use trim because dir is of fixed length
    output = dir%lower()//rfname//fStr(cnt)//'.dat'
    print *, output%center(width)
  end do
end program example_fstring1
