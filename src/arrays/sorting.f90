!> sort provides a framework for searching elements and sorting arrays
!!
module sorting
  USE basic, only: dp, Small
  integer, parameter :: minval_bisection = 100 !< Minimum value of elements for using bisection
  integer, parameter :: max_size_for_insertion_sort = 20 !! max size for using insertion sort.

  !> searchsorted: Find index where an element should be inserted in an array to maintain order.
  !!
  !! Find the index into an ascending sorted array `x` such that, if `elem` was inserted
  !! after the index, the order of `x` would be preserved.
  !!
  !! @note Bisection is used to find the required insertion point if number of elements is higher of
  !! a certain threshold
  !!
  !! @note If `elem` is outside the limits of `x` then:
  !!   - If below the vector then (first-1) is returned.
  !!   - If above the vector then "last index" is returned.
  !!
  !! Examples:
  !!
  !! With real numbers:
  !! ```
  !! ...
  !! USE arrays, only: searchsorted
  !! ...
  !! real(dp), dimension(6) :: a = [1._dp, 3._dp, 5._dp,9._dp, 32._dp, 124._dp]
  !! real(dp) :: elem
  !! integer :: n
  !! elem = 6.34_dp
  !! n = searchsorted(a, elem)
  !! ```
  !! will give `n = 3`,
  !! while
  !! ```
  !! n = searchsorted(a, 0.5_dp)
  !! ```
  !! will give `n = 0`
  !! and
  !! ```
  !! n = searchsorted(a, 1000._dp)
  !! ```
  !! will give `n = 6`
  !!
  interface searchsorted
    module procedure :: searchsorted_dp, searchsorted_i, searchsorted_idp, searchsorted_dpi
  end interface searchsorted

  interface swap
    module procedure :: swap_dp, swap_i
  end interface swap

  private
  public searchsorted, sort

contains

  !> Routine to sort a vector
  !!
  !! Sorting routine, slightly modified from the original sources:
  !! @ref https://github.com/jacobwilliams/fortran-search-and-sort.git
  !!
  !! It uses a simple insertion point algorithm for small arrays (less than 20 elements)
  !! and quicksort algorithm for larger arrays.
  !!
  !! Examples:
  !!```
  !!  real(dp), dimension(4) :: a = [ 1._dp, 0.5_dp, -2._dp, 9._dp]
  !!  call sort(a)
  !!  ! gives: -2 0.5 1 9
  !!  call sort(a, reverse=.TRUE.)
  !!  ! gives: 9 1 0.5 -2
  !!```
  !!
  !! @note Implementation detail: Sorting in descending order is currently implemented by sorting
  !! -vec and later multiplying by (-1) again.  From benchmark appears that there is not an
  !! important efficiency impact. Further study would be welcome.
  subroutine sort(vec, reverse)
    implicit none
    real(dp), dimension(:), intent(inout) :: vec !< Vector to sort
    logical, optional, intent(IN) :: reverse     !< .TRUE. to sort in descending order
    logical :: reverse_
    reverse_ = .FALSE.
    IF (present(reverse)) reverse_ = reverse

    IF (reverse_) vec = -vec
    call quicksort(1, size(vec))
    IF (reverse_) vec = -vec
  contains

    include 'sorting.inc'
  end subroutine sort

  ! Routines used by sorting routine
  subroutine swap_dp(v1, v2)
    implicit none
    real(dp), intent(INOUT) :: v1 !<
    real(dp), intent(INOUT) :: v2 !<
    real(dp) :: tmp
    tmp = v1
    v1 = v2
    v2 = tmp
  end subroutine swap_dp

  subroutine swap_i(v1, v2)
    implicit none
    integer, intent(INOUT) :: v1 !<
    integer, intent(INOUT) :: v2 !<
    integer :: tmp
    tmp = v1
    v1 = v2
    v2 = tmp
  end subroutine swap_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!             searchsorted implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Search in an array of real(dp) elements. @see searchsorted documentation
  pure function searchsorted_dp(x, elem) result(n)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< Array sorted in ascending order
    real(dp), intent(IN) :: elem            !< element to insert
    integer :: n                            !< index of closest edge to the left of elem

    real(dp), parameter :: Delta = Small ! Make comparison between real numbers more robust
    include "bisection.inc"
  end function searchsorted_dp

  !> Search in an array of integer elements. @see searchsorted documentation
  pure function searchsorted_i(x, elem) result(n)
    implicit none
    integer, dimension(:), intent(IN) :: x !< Array sorted in ascending order
    integer, intent(IN) :: elem            !< element to insert
    integer :: n                           !< index of closest edge to the left of elem

    integer, parameter :: Delta = 0 ! Comparing integers is exact
    include "bisection.inc"
  end function searchsorted_i

  pure function searchsorted_idp(x, elem) result(n)
    implicit none
    integer, dimension(:), intent(IN) :: x !< Array sorted in ascending order
    real(dp), intent(IN) :: elem           !< element to insert
    integer :: n                           !< index of closest edge to the left of elem

    real(dp), parameter :: Delta = Small ! Comparing integers to reals is not exact
    include "bisection.inc"
  end function searchsorted_idp

  pure function searchsorted_dpi(x, elem) result(n)
    implicit none
    real(dp), dimension(:), intent(IN) :: x !< Array sorted in ascending order
    integer, intent(IN) :: elem            !< element to insert
    integer :: n                           !< index of closest edge to the left of elem

    real(dp), parameter :: Delta = Small ! Comparing real to integer is not exact
    include "bisection.inc"
  end function searchsorted_dpi

end module sorting
