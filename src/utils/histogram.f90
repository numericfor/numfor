module histogram

  type, public :: histogram
    integer                       :: n = 0
    real(8), dimension(:), pointer  :: range
    real(8), dimension(:), pointer  :: bin
  end type histogram

end module histogram
