!> @file oopstring.f90 file provides an object-oriented type for working with string
!! @date "2019-08-28 13:15:30"

!> Module defining the `fStr` object and its methods. Documentation: @ref modfstring
module fString
  USE basic, only: dp
  USE strings, only: str,&
    & st_upper => upper, st_lower => lower, st_swapcase => swapcase, &
    & st_startswith => startswith, st_endswith => endswith, &
    & st_lstrip => lstrip, st_rstrip => rstrip, st_strip => strip, &
    & st_reverse => reverse, st_count => count_sub, &
    & st_center => center, st_find => find, st_replace => replace
  implicit none

  private
  public fStr, len, str

  !> This type defines a string class, with its methods
  type :: fStr
    character(len=:), allocatable :: val !< This holds the string of characters

  contains

    procedure, private, pass(S) :: init
    procedure, private, pass(S) :: equal
    procedure, private, pass(S) :: not_equal
    procedure :: readf
    procedure :: writef
    procedure, private, pass(S) :: left_multiply_int
    procedure, private :: right_multiply_int
    procedure, private :: concat_ss
    procedure, private :: concat_sc
    procedure, private, pass(S) :: concat_cs

    procedure :: upper
    procedure :: lower
    procedure :: swapcase

    procedure :: reverse
    procedure :: endswith
    procedure :: startswith

    procedure :: lstrip
    procedure :: rstrip
    procedure :: strip

    procedure :: count => count_sub

    procedure :: center
    procedure :: find
    procedure :: replace

    ! Generic interfaces. Operator Overloading
    generic :: Assignment(=) => init
    generic :: Operator(==) => equal
    generic :: Operator(/=) => not_equal
    generic :: Operator(*) => left_multiply_int, right_multiply_int
    generic :: Operator(//) => concat_ss, concat_cs, concat_sc
    generic :: read (formatted) => readf
    generic :: write (formatted) => writef

  end type fStr

  interface str
    module procedure :: char_str
  end interface str

  interface len
    module procedure :: len_str
  end interface len

  interface fStr
    module procedure set_anyvalue
  end interface fStr

contains

  !> char converts a fStr to a standard allocatable character
  function char_str(S) result(output)
    implicit none
    class(fStr), intent(IN) :: S
    character(len=:), allocatable :: output
    output = S%val
  end function char_str

  subroutine init(S, value)
    class(fStr), intent(Out) :: S
    character(len=*), intent(In) :: value
    S%val = str(value)
  end subroutine init

  !> @copydoc strings::upper()
  function upper(S) result(output)
    implicit none
    class(fStr), intent(IN) :: S
    type(fStr) :: output
    output%val = st_upper(S%val)
  end function upper

  !> @copydoc strings::lower()
  function lower(S) result(output)
    implicit none
    class(fStr), intent(IN) :: S
    type(fStr) :: output
    output%val = st_lower(S%val)
  end function lower

  !> @copydoc strings::swapcase()
  function swapcase(S) result(output)
    implicit none
    class(fStr), intent(IN) :: S
    type(fStr) :: output
    output%val = st_swapcase(S%val)
  end function swapcase

  !> @copydoc strings::startswith()
  function startswith(S, prefix) result(y)
    implicit none
    logical :: y
    class(fStr), intent(IN) :: S
    character(len=*), intent(IN) :: prefix
    y = st_startswith(S%val, prefix)
  end function startswith

  !> @copydoc strings::endswith()
  function endswith(S, suffix) result(y)
    implicit none
    logical :: y
    class(fStr), intent(IN) :: S
    character(len=*), intent(IN) :: suffix
    y = st_endswith(S%val, suffix)
  end function endswith

  !> @copydoc strings::reverse()
  function reverse(S) result(output)
    class(fStr), intent(in) :: S
    type(fStr) :: output
    output%val = st_reverse(S%val)
  end function reverse

  !> @copydoc strings::lstrip()
  function lstrip(S, chars) result(output)
    class(fStr), intent(in) :: S
    character(len=*), optional, intent(in) :: chars
    type(fStr) :: output
    output%val = st_lstrip(S%val, chars)
  end function lstrip

  !> @copydoc strings::rstrip()
  function rstrip(S, chars) result(output)
    class(fStr), intent(in) :: S
    character(len=*), optional, intent(in) :: chars
    type(fStr) :: output
    output%val = st_rstrip(S%val, chars)
  end function rstrip

  !> @copydoc strings::strip()
  function strip(S, chars) result(output)
    class(fStr), intent(in) :: S
    character(len=*), optional, intent(in) :: chars
    type(fStr) :: output
    output%val = st_strip(S%val, chars)
  end function strip

  !> @copydoc strings::count_sub()
  function count_sub(S, sub, start, end) result(output)
    class(fStr), intent(in) :: S
    character(len=*), intent(in) :: sub
    integer, intent(in), optional :: start
    integer, intent(in), optional :: end
    integer :: output
    output = st_count(S%val, sub, start, end)
  end function count_sub

  !> @copydoc strings::center()
  function center(S, width, fillchar) result(output)
    class(fStr), intent(in) :: S
    integer, intent(in) :: width
    character(len=1), optional, intent(in) :: fillchar
    type(fStr) :: output
    output%val = st_center(S%val, width, fillchar)
  end function center

  !> @copydoc strings::find()
  function find(S, sub, start, end) result(output)
    class(fStr), intent(IN) :: S
    integer :: output
    character(len=*), intent(IN) :: sub
    integer, optional, intent(IN) :: start
    integer, optional, intent(IN) :: end
    output = st_find(S%val, sub, start, end)
  end function find

  !> @copydoc strings::replace()
  function replace(S, old, new, count) result(output)
    class(fStr), intent(IN) :: S
    type(fStr) :: output
    character(len=*), intent(IN) :: old
    character(len=*), intent(IN) :: new
    integer, optional, intent(IN) :: count
    output%val = st_replace(S%val, old, new, count)
  end function replace

  function set_anyvalue(value) result(S)
    type(fStr) :: S
    class(*), intent(in) :: value
    select type (value)
    type is (fStr)
      S%val = value%val
    type is (character(len=*))
      S%val = value
    type is (integer)
      S%val = str(value)
    type is (real)
      S%val = str(value)
    type is (real(dp))
      S%val = str(value)
    end select
  end function set_anyvalue

  subroutine writef(S, unit, iotype, v_list, iostat, iomsg)
    class(fStr), intent(in)     :: S   !< Object to write.
    integer, intent(in)         :: unit      !< Internal unit to write to.
    character(*), intent(in)    :: iotype    !< LISTDIRECTED or DTxxx
    integer, intent(in)         :: v_list(:) !< parameters from fmt spec.
    integer, intent(out)        :: iostat    !< non zero on error, etc.
    character(*), intent(inout) :: iomsg     !< define if iostat non zero.
    integer         :: vlist
    ! These are dummy statements, to avoid "unused-dummy-argument" warning/errors
    vlist = iachar(iotype)
    vlist = size(v_list)
    if (allocated(S%val)) then
      write (unit, "(a)", IOSTAT=iostat, IOMSG=iomsg) S%val
    else
      write (unit, "(a)", IOSTAT=iostat, IOMSG=iomsg) ' '
    endif
  end subroutine

  subroutine readf(S, unit, iotype, v_list, iostat, iomsg)
    class(fStr), intent(inout)  :: S
    integer, intent(in)         :: unit
    character(*), intent(in)    :: iotype
    integer, intent(in)         :: v_list(:)
    integer, intent(out)        :: iostat
    character(*), intent(inout) :: iomsg
    character(len=8 * 1024)     :: temp ! long lines
    integer         :: vlist
    ! These are dummy statements, to avoid "unused-dummy-argument" warning/errors
    vlist = iachar(iotype)
    vlist = size(v_list)

    read (unit, fmt='(a)', iostat=iostat, iomsg=iomsg) temp
    S%val = trim(adjustl(temp))
  end subroutine readf

  !> Calculates the length of the string, like the intrinsic function
  pure function len_str(S) result(output)
    implicit none
    class(fStr), intent(in) :: S !< Original string
    integer :: output            !< Length of the string
    output = len(S%val)
  end function len_str

  !> equal: comparison between two strings
  function equal(S, other) result(output)
    logical :: output           !< .True. if equal, .False. if not
    class(fStr), intent(IN) :: S !< Left-side of comparison
    class(fStr), intent(IN) :: other !< Right-side of comparison
    output = (S%val == other%val)
  end function equal

  !> not_equal: comparison between two strings
  function not_equal(S, other) result(output)
    logical :: output           !< .False. if equal, .True. if not
    class(fStr), intent(IN) :: S !< Left-side of comparison
    class(fStr), intent(IN) :: other !< Right-side of comparison
    output = (S%val /= other%val)
  end function not_equal

  !> Multiply by an integer
  function left_multiply_int(n, S) result(output)
    implicit none
    class(fStr), intent(IN) :: S      !< Original string
    integer, intent(IN) :: n          !< number of repetitions
    type(fStr) :: output              !< String multiplied n times
    output%val = repeat(S%val, n)
  end function left_multiply_int

  !> @copydoc left_multiply_int
  function right_multiply_int(S, nr) result(output)
    implicit none
    class(fStr), intent(IN) :: S
    integer, intent(IN) :: nr
    type(fStr) :: output
    output%val = repeat(S%val, nr)
  end function right_multiply_int

  function concat_sc(S, value) result(output)
    class(fStr), intent(in) :: S          !< First string
    character(len=*), intent(in) :: value !< Second string
    type(fStr) :: output                  !< Joined string
    output%val = S%val//value
  end function concat_sc

  function concat_cs(value, S) result(output)
    class(fStr), intent(in) :: S          !< First string
    character(len=*), intent(in) :: value !< Second string
    type(fStr) :: output                  !< Joined string
    output%val = value//S%val
  end function concat_cs

  function concat_ss(S, value) result(output)
    class(fStr), intent(in) :: S     !< first string
    class(fStr), intent(in) :: value !< Second string
    type(fStr) :: output             !< Joined string
    output%val = S%val//value%val
  end function concat_ss

end module fString
