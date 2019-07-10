program example_strings1
  USE strings, only: str, lower, center
  implicit none
  character(len=1024) :: dir
  character(len=:), allocatable :: rfname
  character(len=:), allocatable :: output
  character(len=:), allocatable :: title
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
  read (*, '(A)') dir

  dir = trim(dir)//'/'          ! Add separator

  width = 40 + len(trim(dir))   ! Used for centering
  title = ' These will the files created: '

  ! Set some parameters
  param1 = 1.3_8
  param2 = 0.07_8
  param3 = 42

  ! rfname holds the root of the filename
  rfname = "out"//str(param1)//'cm'//str(param3)//'eV_'
  print *, ''
  print '(A)', center(title, width, '=')
  do cnt = 1, Nmax
    ! Create full filename
    ! Here we have to use trim because dir is of fixed length
    output = lower(trim(dir))//rfname//str(cnt)//'.dat'
    print '(A)', center(output, width)
  end do
end program example_strings1
