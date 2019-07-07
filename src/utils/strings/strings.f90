!> This module defines functions to manipulate strings of characters.
!!
!! It provides functions similar to those present in Python
module strings

  character(*), Parameter :: ascii_lowercase = 'abcdefghijklmnopqrstuvwxyz'
  character(*), Parameter :: ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  Public
contains

  !> @brief Return a copy of the string converted to uppercase.
  function str_upper(S) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    character(len=:), allocatable :: Sout !< String converted to uppercase
    integer :: i
    integer :: n
    Sout = S
    do i = 1, len_trim(S)
      n = INDEX(ascii_lowercase, Sout(i:i)) ! Find the location of letter
      IF (n /= 0) Sout(i:i) = ascii_uppercase(n:n) ! If current substring is a lower case letter, make it upper case
    end do
  end function str_upper

  !> @brief Return a copy of the string converted to lowercase.
  function str_lower(S) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    character(len=:), allocatable :: Sout !< String converted to uppercase
    integer :: i
    integer :: n
    Sout = S
    do i = 1, len_trim(S)
      n = index(ascii_uppercase, Sout(i:i)) ! Find the location of letter
      IF (n /= 0) Sout(i:i) = ascii_lowercase(n:n) ! If current substring is a upper case letter, make it lower case
    end do
  end function str_lower

  !> @brief Return str with case of letters swapped
  function str_swapcase(S) result(Sout)
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
  end function str_swapcase

  !> @brief Reverse a string.
  function str_reverse(S) result(Sout)
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
  end function str_reverse

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
  !! @note that differs from python method in that does not accept a tuple
  function startswith(S, prefix) result(y)
    implicit none
    character(len=*), intent(IN) :: S      !< Original string
    character(len=*), intent(IN) :: prefix !< substring to test
    logical :: y                !< .True. if `S` starts with `sub`
    y = (S(:len(prefix)) == prefix)
  end function startswith

  !> This function returns a copy of the string with leading chars removed
  !!
  function str_lstrip(S, chars) result(Sout)
    implicit none
    character(len=*), intent(IN) :: S               !< Original string
    character(len=*), optional, intent(IN) :: chars !< chars to remove from S
    character(len=:), allocatable :: Sout                        !< String with chars removed

    integer :: n
    n = verify(S, chars)
    Sout = S(n:)
  end function str_lstrip

  !> This function returns a copy of the string with trailing chars removed
  !!
  function str_rstrip(S, chars) result(Sout)
    implicit none
    character(len=*), intent(IN) :: S               !< Original string
    character(len=*), optional, intent(IN) :: chars !< chars to remove from S
    character(len=:), allocatable :: Sout !< String with chars removed
    integer :: n

    n = verify(S, chars, back=.True.)
    Sout = S(:n)
  end function str_rstrip

  !> This function returns a copy of the string with leading and trailing chars removed
  !!
  function str_strip(S, chars) result(Sout)
    implicit none
    character(len=*), intent(IN) :: S               !< Original string
    character(len=*), optional, intent(IN) :: chars !< chars to remove from S
    character(len=:), allocatable :: Sout !< String with chars removed
    Sout = str_lstrip(str_rstrip(S, chars), chars)
  end function str_strip

  !> Return the number of occurrences of substring sub in string S[start:end].
  function str_count(S, sub, start, end) result(count)
    implicit none
    character(len=*), intent(in) :: S !< Original string to which count occurrences
    character(len=*), intent(in) :: sub !< substring to count
    integer, intent(in), optional :: start !< initial position to consider
    integer, intent(in), optional :: end   !< final position to consider
    integer :: count                       !< Number of occurrences of sub in S
    integer :: i
    integer :: start_, end_
    count = 0
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
      count = count + 1
    end do
  end function str_count

  !> @brief Return .True. if sub is present in string, .False. otherwise
  function str_issub(S, sub) result(is_in)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    character(len=*), intent(in) :: sub !< substring to find
    logical :: is_in

    is_in = (str_count(S, sub) > 0)
  end function str_issub

  !> Pad a string with zeroes ("0") to specified width. If width is <= input
  !!        string width, then the original string is returned.
  !! @note The string is never truncated
  function str_zfill(S, width) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    integer, intent(in) :: width        !< width of padded string
    character(len=:), allocatable :: Sout
    integer :: lS
    lS = len_trim(S)
    if (width <= lS) then
      Sout = S
    else
      Sout = repeat("0", width - lS)//S
    endif
  end function str_zfill

  !> Center a string to a specified width.  The default character to fill in the
  !>        centered string is a blank character.
  function str_center(S, width, fillchar) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< Original string
    integer, intent(in) :: width      !< Total width of centered string
    character(len=1), optional, intent(in) :: fillchar !< Padding character (default to space)
    character(len=:), allocatable :: Sout !< A centered string of length width.
    character(len=1) :: fchar_
    integer :: Lfill

    fchar_ = " "
    IF (present(fillchar)) fchar_ = fillchar

    ! Sout = trim(S)

    if (width > len_trim(S)) then
      Lfill = width - len_trim(S)
      if (mod(Lfill, 2) == 0) then
        Sout = repeat(fchar_, Lfill / 2)//trim(S)//repeat(fchar_, Lfill / 2)
      else
        Sout = repeat(fchar_, (Lfill / 2) + 1)//trim(S)//repeat(fchar_, Lfill / 2)
      endif
    endif
  end function str_center

  !> Return the lowest index in S where substring sub is found
  function str_find(S, sub, start, end) result(nindex)
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
  end function str_find

  !> Return a copy with all occurrences of substring old replaced by new
  function str_replace(S, old, new, count) result(Sout)
    implicit none
    character(len=*), intent(in) :: S !< original string
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
      i = index(Sout, old)
      IF (i == 0) exit
      Sout = Sout(:i - 1)//new//Sout(i + nold:)
      c = c + 1
      IF (c == count_) exit
    end do

  end function str_replace

  ! !> @brief Split string based on a character delimiter and return string given by the
  ! !>        column number.
  ! !> @param[in] str - string to work on
  ! !> @param[in] delim -
  ! !> @param[in] col - delimited column to return
  ! !> @return string
  ! function str_split(S, sep, maxsplit) result(Sout)
  !   implicit none
  !   character(len=*), intent(in) :: S             !< original string
  !   character(len=1), optional, intent(in) :: sep !< character delimiter
  !   integer, optional, intent(in) :: maxsplit     !< Maximum number of splits to do
  !   character(len=:), allocatable :: Sout
  !   character(len=:), allocatable :: ctemp
  !   character(len=:), allocatable :: cwork
  !   integer :: i, cnt, lastpos
  !   cnt = 0
  !   lastpos = 0
  !   Sout = ""
  !   if (maxsplit .le. 0 .or. maxsplit .gt. (str_count(S, sep) + 1) .or. str_count(S, sep) .eq. 0) then
  !     Sout = S
  !     return
  !   endif
  !   ctemp = S//sep
  !   do i = 1, len_trim(ctemp)
  !     if (ctemp(i:i) .eq. sep) then
  !       cnt = cnt + 1
  !       if (cnt .eq. maxsplit) then
  !         cwork = ctemp(lastpos + 1:i - 1)
  !         exit
  !       endif
  !       lastpos = i
  !     endif
  !   end do
  !   Sout = cwork
  ! end function str_split

  ! ! -------------------------------------------------------------------------------------
  ! ! Function: str_uniq
  ! !> @brief Removed duplicative entries from a \b delimited string.
  ! !> @param[in] str - string to work on
  ! !> @param[in] delim - character delimiter
  ! !> @return modified string
  ! ! -------------------------------------------------------------------------------------
  ! function str_uniq(str, delim) result(strout)
  !   implicit none
  !   character(len=*), intent(in) :: str
  !   character(len=1), intent(in) :: delim
  !   character(len=:), allocatable :: strout
  !   character(len=:), allocatable :: ctemp
  !   character(len=:), allocatable :: col
  !   integer :: n, nn, ncols, ndelims, nuniq, strlen
  !   logical(kind=1) :: strfound
  !   ncols = 0
  !   ndelims = 0
  !   nuniq = 0
  !   strlen = 0
  !   ctemp = trim(adjustl(str))
  !   strlen = len(ctemp)
  !   ndelims = str_count(ctemp, delim)
  !   ncols = ndelims + 1
  !   strout = ""
  !   do n = 1, ncols
  !     col = str_split(ctemp, delim, n)
  !     if (n .eq. 1) then
  !       strout = col//delim
  !       nuniq = len(strout)
  !     endif
  !     strfound = .false.
  !     do nn = 1, nuniq
  !       if (col .eq. strout(nn:(nn + len(col)) - 1)) then
  !         strfound = .true.
  !         exit
  !       endif
  !     end do
  !     if (.not. strfound) then
  !       strout = strout//col//delim
  !       nuniq = len(strout)
  !     endif
  !   end do
  !   if (strout(1:1) .eq. delim) strout = strout(2:len(strout))
  !   strlen = len(strout)
  !   if (strout(strlen:strlen) .eq. delim) strout = strout(1:len(strout) - 1)
  ! end function str_uniq

end module strings
