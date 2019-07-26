!> sort provides a framework for searching elements and sorting arrays
!!
module sort
  USE basic, only: dp, Small
  integer, parameter :: minval_bisection = 100 !< Minimum value of elements for using bisection
  !! Otherwise, will use simple sequential lookup

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
    module procedure :: searchsorted_dp, searchsorted_i
  end interface searchsorted

  private
  public searchsorted

contains

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

end module sort
