!> @file strings.f90 provides routines for common string manipulation

!> This module defines functions to manipulate strings of characters.
!! Documentation: @ref modstrings
module strings

  character(*), Parameter :: ascii_lowercase = 'abcdefghijklmnopqrstuvwxyz'
  character(*), Parameter :: ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  Public
  Private :: c2str, i2str, r2str, dp2str

  !> `str()` converts a number (integer or real) to a string
  !!
  !! Use:
  !! ```
  !! USE strings, only: str
  !! print *, 'I2str: - int:', 123, '- string: |'//str(123)//'|'
  !! print *, 'r2str:  - real:', 123.121000_4, '- string: |'//str(123.12100_4)//'|'
  !! print *, 'dp2str: - real:', 123.1221_8, '- string: |'//str(123.1221_8)//'|'
  !! print *, 'dp2str: - real:', -123.121_8, '- string: |'//str(-123.121_8)//'|'
  !! print *, 'dp2str: - real:', -123.1221e-12_8, '- string: |'//str(-123.1221e-12_8)//'|'
  !! ```
  !! Gives:
  !! ```
  !! !    I2str: - int:         123 - string: |123|
  !! !    r2str:  - real:   123.121002     - string: |123.121002|
  !! !    dp2str: - real:   123.12210000000000      - string: |123.1221|
  !! !    dp2str: - real:  -123.12100000000000      - string: |-123.121|
  !! !    dp2str: - real:  -1.2312210000000001E-010 - string: |-0.12312210000000001E-009|
  !! !
  !! ```
  interface str
    module procedure :: c2str, i2str, r2str, dp2str
  end interface str

  !> Repeat string `S` a number of times
  !!
  !! Example:
  !! ```
  !! ! Define variables ...
  !! c1 = "hola"
  !! print *, 3 * c1
  !! ! produces: "holaholahola"```
  !! print *, c1 * 3
  !! ! produces: "holaholahola"```
  !!
  interface multiply
    module procedure :: left_multiply_int, right_multiply_int
  end interface multiply

contains

  !> @brief Returns a copy of the string converted to uppercase.
  function upper(S) result(Sout)
    implicit none
    character(len=*), intent(in) :: S     !< Original string
    character(len=:), allocatable :: Sout !< String converted to uppercase
    integer :: i
    integer :: n
    Sout = S
    do i = 1, len_trim(S)
      n = INDEX(ascii_lowercase, Sout(i:i))
      IF (n /= 0) Sout(i:i) = ascii_uppercase(n:n)
    end do
  end function upper

  !> @brief Returns a copy of the string converted to lowercase.
  function lower(S) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    character(len=:), allocatable :: Sout !< String converted to uppercase
    integer :: i
    integer :: n
    Sout = S
    do i = 1, len_trim(S)
      n = index(ascii_uppercase, Sout(i:i))
      IF (n /= 0) Sout(i:i) = ascii_lowercase(n:n)
    end do
  end function lower

  !> @brief Return str with case of letters swapped
  function swapcase(S) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    character(len=:), allocatable :: Sout !< String with cases swapped
    integer :: i, nl, nu
    Sout = S
    do i = 1, len(S)
      nu = index(ascii_uppercase, Sout(i:i))
      if (nu /= 0) then         ! if is uppercase
        Sout(i:i) = ascii_lowercase(nu:nu)
      else
        nl = index(ascii_lowercase, Sout(i:i)) ! is lower
        IF (nl /= 0) Sout(i:i) = ascii_uppercase(nl:nl)
      end if
    end do
  end function swapcase

  !> @brief Reverse a string.
  function reverse(S) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    character(len=:), allocatable :: Sout !< Reversed string
    integer :: i, n
    Sout = S
    n = len(S)
    do i = 1, len(S)
      Sout(i:i) = S(n:n)
      n = n - 1
    end do
  end function reverse

  !> Return True if S starts with the specified prefix, False otherwise.
  !!
  !! @note that differs from python method in that does not accept a tuple
  function endswith(S, suffix) result(y)
    implicit none
    character(len=*), intent(IN) :: S      !< Original string
    character(len=*), intent(IN) :: suffix !< substring to test
    logical :: y                !< .True. if `S` ends with `sub`.
    y = (S(len(S) - len(suffix):) == suffix)
  end function endswith

  !> Return True if S starts with the specified prefix, False otherwise.
  !!
  !! @note that differs from python method in that does not accept a tuple as prefix
  function startswith(S, prefix) result(y)
    implicit none
    character(len=*), intent(IN) :: S      !< Original string
    character(len=*), intent(IN) :: prefix !< substring to test
    logical :: y                !< True if `S` starts with `prefix`
    y = (S(:len(prefix)) == prefix)
  end function startswith

  !> This function returns a copy of the string with leading chars removed
  !!
  function lstrip(S, chars) result(Sout)
    implicit none
    character(len=*), intent(IN) :: S               !< Original string
    character(len=*), optional, intent(IN) :: chars !< chars to remove from S
    character(len=:), allocatable :: Sout                        !< String with chars removed
    integer :: n
    character(len=:), allocatable :: ch_

    ch_ = ' '
    if (Present(chars)) ch_ = chars

    n = verify(S, ch_)
    Sout = S(n:)
  end function lstrip

  !> This function returns a copy of the string with trailing chars removed
  !!
  function rstrip(S, chars) result(Sout)
    implicit none
    character(len=*), intent(IN) :: S               !< Original string
    character(len=*), optional, intent(IN) :: chars !< chars to remove from S
    character(len=:), allocatable :: Sout !< String with chars removed
    integer :: n
    character(len=:), allocatable :: ch_

    ch_ = ' '
    if (Present(chars)) ch_ = chars

    n = verify(S, chars, back=.True.)
    Sout = S(:n)
  end function rstrip

  !> This function returns a copy of the string with leading and trailing chars removed
  !!
  function strip(S, chars) result(Sout)
    implicit none
    character(len=*), intent(IN) :: S               !< Original string
    character(len=*), optional, intent(IN) :: chars !< chars to remove from S
    character(len=:), allocatable :: Sout !< String with chars removed
    character(len=:), allocatable :: ch_

    ch_ = ' '
    if (Present(chars)) ch_ = chars

    Sout = lstrip(rstrip(S, ch_), ch_)
  end function strip

  !> Return the number of occurrences of substring sub in string S[start:end].
  function count_sub(S, sub, start, end) result(cnt)
    implicit none
    character(len=*), intent(in) :: S   !< Original string
    character(len=*), intent(in) :: sub !< substring to count
    integer, intent(in), optional :: start !< initial position to consider
    integer, intent(in), optional :: end !< final position to consider
    integer :: cnt              !< Number of occurrences of sub in S
    integer :: i
    integer :: start_, end_
    cnt = 0
    start_ = 1
    end_ = len(trim(S))

    if (Present(start)) then
      IF (start >= 1) start_ = start
    end if
    if (Present(end)) then
      IF (end <= end_) end_ = end
    end if
    IF (start_ > end_) return

    i = start_
    do
      i = index(S(start_:end_), sub)
      IF (i == 0) return
      start_ = start_ + i
      cnt = cnt + 1
    end do
  end function count_sub

  !> Multiply a string by an integer
  function left_multiply_int(n, S) result(Sout)
    implicit none
    character(len=*), intent(in) :: S     !< Original string
    integer, intent(in) :: n              !< number of repetitions
    character(len=:), allocatable :: Sout !< String multiplied n times
    integer :: i
    Sout = ''
    do i = 1, n
      Sout = Sout//S
    end do
  end function left_multiply_int

  !> Multiply a string by an integer
  function right_multiply_int(S, nr) result(Sout)
    implicit none
    character(len=*), intent(in) :: S     !< Original string
    integer, intent(in) :: nr             !< number of repetitions
    character(len=:), allocatable :: Sout !< String multiplied n times
    Sout = left_multiply_int(nr, S)
  end function right_multiply_int

  !> @brief Returns .True. if sub is present in S, .False. otherwise
  function issub(S, sub) result(is_in)
    implicit none
    character(len=*), intent(in) :: S   !< Original string
    character(len=*), intent(in) :: sub !< substring to find
    logical :: is_in                    !< True if `sub` is present in `S`

    is_in = (count_sub(S, sub) > 0)
  end function issub

  !> Return a left-justified string of length width.
  !! @note The string is never truncated
  function ljust(S, width, fillchar) result(Sout)
    implicit none
    character(len=*), intent(in) :: S     !< Original string
    integer, intent(in) :: width          !< width of padded string
    character(len=1), optional, intent(in) :: fillchar !< Char used for filling
    character(len=:), allocatable :: Sout !< 0-padded string
    character(len=1) :: ch_
    integer :: lS
    ch_ = ' '; IF (present(fillchar)) ch_ = fillchar
    lS = len_trim(S)
    if (width <= lS) then
      Sout = S
    else
      Sout = repeat(ch_, width - lS)//S
    endif
  end function ljust

  !> Return a right-justified string of length width.
  !! @note The string is never truncated
  function rjust(S, width, fillchar) result(Sout)
    implicit none
    character(len=*), intent(in) :: S     !< Original string
    integer, intent(in) :: width          !< width of padded string
    character(len=1), optional, intent(in) :: fillchar !< Char used for filling
    character(len=:), allocatable :: Sout !< 0-padded string
    character(len=1) :: ch_
    integer :: lS
    ch_ = ' '; IF (present(fillchar)) ch_ = fillchar
    lS = len_trim(S)
    if (width <= lS) then
      Sout = S
    else
      Sout = S//repeat(ch_, width - lS)
    endif
  end function rjust

  !> Pad a string with zeroes ("0") to specified width. If width is <= input
  !!        string width, then the original string is returned.
  !! @note The string is never truncated
  function zfill(S, width) result(Sout)
    implicit none
    character(len=*), intent(in) :: S     !< Original string
    integer, intent(in) :: width          !< width of padded string
    character(len=:), allocatable :: Sout !< 0-padded string
    integer :: lS
    lS = len_trim(S)
    if (width <= lS) then
      Sout = S
    else
      if (startswith(S, '-')) then
        Sout = "-"//repeat("0", width - lS)//S(2:)
      else
        Sout = repeat("0", width - lS)//S
      end if

    endif
  end function zfill

  !> Center a string to a specified width.  The default character to fill in the
  !!        centered string is a blank character.
  function center(S, width, fillchar) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    integer, intent(in) :: width      !< Total width of centered string
    character(len=1), optional, intent(in) :: fillchar !< Padding character (default to space)
    character(len=:), allocatable :: Sout !< A centered string of length width.
    character(len=1) :: ch_
    integer :: Lfill

    ch_ = " "
    IF (present(fillchar)) ch_ = fillchar

    if (width > len_trim(S)) then
      Lfill = width - len_trim(S)
      if (mod(Lfill, 2) == 0) then
        Sout = repeat(ch_, Lfill / 2)//trim(S)//repeat(ch_, Lfill / 2)
      else
        Sout = repeat(ch_, (Lfill / 2) + 1)//trim(S)//repeat(ch_, Lfill / 2)
      endif
    endif
  end function center

  !> Return the lowest index in S where substring sub is found
  function find(S, sub, start, end) result(nindex)
    implicit none
    integer :: nindex !< position where found
    character(len=*), intent(IN) :: S      !< original string
    character(len=*), intent(IN) :: sub    !< substring to find
    integer, optional, intent(IN) :: start !< initial position to search
    integer, optional, intent(IN) :: end   !< final position to search

    integer :: start_, end_
    nindex = 0
    start_ = 1
    end_ = len(S)

    if (Present(start)) then
      IF (start >= 1) start_ = start
    end if
    if (Present(end)) then
      IF (end <= end_) end_ = end
    end if
    IF (start_ > end_) return

    nindex = index(S(start_:end_), sub)
    IF (nindex /= 0) nindex = nindex + (start_ - 1)
  end function find

  !> Return a copy with all occurrences of substring old replaced by new
  function replace(S, old, new, count) result(Sout)
    implicit none
    character(len=*), intent(in) :: S   !< original string
    character(len=*), intent(in) :: old !< substring to replace
    character(len=*), intent(in) :: new !< substring to substitute from old
    integer, optional, intent(in) :: count !< Maximum number of occurrences to replace
    character(len=:), allocatable :: Sout !< New string created

    integer :: i, c, nold
    integer :: count_

    Sout = S
    nold = len(old)

    count_ = -1; IF (Present(count)) count_ = count

    c = 0
    do
      IF (c == count_) exit
      i = index(Sout, old)
      IF (i == 0) exit
      Sout = Sout(:i - 1)//new//Sout(i + nold:)
      c = c + 1
    end do
  end function replace

  ! Casts a character string into an string
  function c2str(Sin) result(Sout)
    implicit none
    character(len=*), intent(IN) :: Sin !< Original string
    character(len=:), allocatable :: Sout !< String converted
    Sout = Sin
  end function c2str

  ! Casts a default integer into an string
  function i2str(nin) result(Sout)
    implicit none
    integer, intent(IN) :: nin !< Integer to convert
    character(len=:), allocatable :: Sout !< String converted
    character(len=10) :: S_
    write (S_, '(i8)') nin
    Sout = strip(S_)
  end function i2str

  ! Casts a real(8) into an string
  function dp2str(rin) result(Sout)
    implicit none
    real(8), intent(IN) :: rin !< number to convert
    character(len=:), allocatable :: Sout !< String converted
    character(len=25) :: S_
    if ((rin > 1.e-6) .and. (rin < 1.e6)) then
      write (S_, '(f18.8)') rin
    else
      write (S_, '(es21.8)') rin
    end if

    Sout = lower(rstrip(lstrip(S_), '0 '))
  end function dp2str

! Casts a real(4) into an string
  function r2str(rin) result(Sout)
    implicit none
    real(4), intent(IN) :: rin !< number to convert
    character(len=:), allocatable :: Sout !< String converted
    character(len=21) :: S_
    write (S_, '(g0)') rin
    Sout = lower(rstrip(S_, '0 '))
  end function r2str

end module strings
