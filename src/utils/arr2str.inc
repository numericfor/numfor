  ! This should be common
  character(len=:), allocatable :: Sout !< String created
  !! Examples:
  !!
  character(len=:), allocatable :: S_
  integer :: n, i, j, k
  integer :: lol                    ! length of line
  integer, parameter :: maxlol = 70 ! maximum length of line
  n = size(vec)
  allocate (character(n * 15) :: Sout)
  j = 2
  lol = 1
  Sout = "[ "
  do i = 1, n
    S_ = str(vec(i))
    k = len_trim(S_) + 2
    lol = lol + k
    if (lol >= maxlol) then
      Sout = Sout(:j)//nl//"  "
      lol = 0
      j = j + 3
    end if
    Sout = Sout(:j)//strip(S_)//', '
    j = j + k
  end do
  Sout = Sout(:j - 2)//"]"
