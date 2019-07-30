!> The module `histograms` provides infrastructure for the calculation of histograms. i.e: a count of frequency
!!
!!
module histograms

  USE basic, only: dp, Zero, Small, print_msg
  USE grids, only: linspace, mean, std
  USE sorting, only: searchsorted

  !> type histogram holds the data from an histogram
  type, public :: histog
    integer :: n = 0            !< Total number of values parsed
    real(dp), dimension(:), allocatable :: hist !< Values of the histogram
    real(dp), dimension(:), allocatable :: bin_edges !< Array of size (len(hist)+1)
  contains

    procedure :: clean => clean_histog
  end type histog

  public :: histogram
  private

contains
  !> clean-up
  subroutine clean_histog(h)
    implicit none
    class(histog), intent(INOUT) :: h !< histogram to clean
    IF (allocated(h%bin_edges)) deallocate (h%bin_edges)
    IF (allocated(h%hist)) deallocate (h%hist)
    h%n = 0
  end subroutine clean_histog

  !> Computes the histogram of an array of data
  !!
  !! @note It counts the number of occurrences in the interval [bins(i), bins(i+1) )
  !! which is closed on the left and open on the right (does not include the upper limit)
  !! except for the last one which includes both limits
  !! @note When neither `bins` nor `Nbins` are present the routine will calculate automatically
  !! the number of bins from the input data.
  !!
  !! Examples:
  !! --------
  !! ```
  !! h = histogram(linspace(0.5_dp, 5.5_dp, 6), Nbins=5)
  !! print "(A, 6(1x,g0.2))", 'Bin Edges:', h%bin_edges
  !! print "(A, 5(1x,g0.2))", 'Histogram:', h%hist
  !! ```
  !! produces:
  !! ```
  !!   Bin Edges: 0.50 1.5 2.5 3.5 4.5 5.5
  !!   Histogram: 1.0 1.0 1.0 1.0 1.0
  !! ```
  !!
  !! ```
  !! h = histogram([1._dp, 2._dp, 1._dp], bins=[0._dp, 1._dp, 2._dp, 3._dp])
  !! print "(A, 4(1x,g0.2))", 'Bin Edges:', h%bin_edges
  !! print "(A, 3(1x,g0.2))", 'Histogram:', h%hist
  !! ```
  !! produces:
  !! ```
  !!   Bin Edges: 0.0 1.0 2.0 3.0
  !!   Histogram: 0.0 2.0 1.0
  !! ```
  !!
  function histogram(a, Nbins, bins, range, weights, density) result(h)
    implicit none
    type(histog) :: h !< The histogram construct.
    real(dp), dimension(:), target, intent(IN) :: a !< Input data
    integer, optional, intent(IN) :: Nbins !< Number of equal-width bins to use
    real(dp), dimension(:), optional, intent(IN) :: bins  !< A monotonically increasing array of bin edges,
    !! including the rightmost edge, allowing for non-uniform bin widths.
    real(dp), dimension(2), optional, intent(IN) :: range !< The lower and upper range of the bins.
    !! If not provided uses `min` and `max`
    real(dp), dimension(size(a)), optional, intent(IN) :: weights !< Array of weights,
    !! of the same size as `a`.  Each value in `a` only contributes its associated weight towards the
    !! bin count (instead of 1)
    logical, optional, intent(IN) :: density !< If ``True``, the result is the probability *density*
    !! function at the bin, normalized such that the *integral* over the range is 1.
    real(dp), dimension(:), pointer :: p
    logical :: d_
    integer :: Nbins_
    integer :: N, nn
    integer :: i, j, k
    integer :: status
    integer, parameter :: BLOCK = 65536
    real(dp):: weight
    real(dp) :: first_edge, last_edge
    character(len=200):: msg

    d_ = .FALSE.; IF (present(density)) d_ = density
    weight = 1.0_dp
    h%n = size(a)
    call get_bins(a, Nbins, bins, range, h)
    Nbins_ = size(h%bin_edges) - 1

    ! Allocate memory for hist component
    allocate (h%hist(Nbins_), stat=status, errmsg=msg)
    IF (status /= 0) call print_msg('Error allocating hist:'//trim(msg), 'histogram')

    first_edge = h%bin_edges(1); last_edge = h%bin_edges(Nbins_ + 1)
    h%hist = Zero

    nn = 0
    ! The use of a BLOCK does not appear to be very important for speed. May be for memory use.
    blocks: do concurrent(k=1:size(a):BLOCK)
      ! blocks: do k = 1, size(a), BLOCK
      p => a(k:)
      N = min(h%n - k + 1, BLOCK) ! Number of elements in block

      do concurrent(j=1:N)
        i = searchsorted(h%bin_edges, p(j))
        if ((i > 0) .and. (i <= Nbins_)) then
          IF (present(weights)) weight = weights(j)
          h%hist(i) = h%hist(i) + weight
          nn = nn + 1
        end if
      end do
    end do blocks

    h%n = nn
    IF (d_) h%hist = h%hist / (h%bin_edges(2:) - h%bin_edges(1:Nbins_)) / sum(h%hist)

  end function histogram

  ! TODO: We need to test which automatic method works best
  ! Discussion. In order to make the last bin a closed-interval on both sides:
  ! [bins(N-2), bins(N-1)] which included both limits
  ! as opposed to [bins(N-2), bins(N-1)], which is closed on the left and open on the right
  ! we add a small quantity to the right edge (2*epsilon)
  ! where epsilon is the smallest difference between two real(dp) numbers
  subroutine get_bins(a, Nbins, bins, range, h)
    implicit none

    real(dp), dimension(:), intent(IN) :: a !<
    integer, optional, intent(IN) :: Nbins !<
    real(dp), dimension(:), optional, intent(IN) :: bins !<
    real(dp), dimension(2), optional, intent(IN) :: range !<
    type(histog), intent(INOUT) :: h
    real(dp) :: first_edge, last_edge, width
    integer :: Nbins_
    integer :: status
    character(len=200) :: msg

    if (present(range)) then     ! Clip to range if present
      IF (range(1) > range(2)) call print_msg('Max must be larger than min', 'histogram')
      first_edge = range(1)
      last_edge = range(2)
    else                        ! Otherwise use extreme values
      first_edge = minval(a)
      last_edge = maxval(a)
    end if

    ! Then, we correct the limits if they are equal
    if (first_edge == last_edge) then
      first_edge = first_edge - 0.5_dp
      last_edge = last_edge + 0.5_dp
    end if

    ! Finally we add a correction in order to compare real numbers
    first_edge = first_edge - Small
    last_edge = last_edge + Small

    if (Present(bins)) then
      ! Check that they are ordered
      IF (any((bins(:size(bins) - 1) > bins(2:)))) &
        & call print_msg('Bins must increase monotonically', 'histogram')

      Nbins_ = size(bins)
      h%bin_edges = bins
      ! We extend the last one in epsilon to make it inclusive
      h%bin_edges(Nbins_) = h%bin_edges(Nbins_) + 1.5_dp * epsilon(1._dp)
      return

    else if (present(Nbins)) then
      Nbins_ = Nbins + 1
    else
      ! We try to guess the "optimal bin edges"
      ! TODO: Here we should check wich method is best
      ! width = width_sturges(a)
      ! width = width_rice(a)
      width = width_doane(a)
      Nbins_ = ceiling(abs((last_edge - first_edge)) / width)
    end if

    allocate (h%bin_edges(Nbins_), stat=status, errmsg=msg)
    IF (status /= 0) call print_msg('Error allocation bin_edges: '//msg, 'histogram')
    h%bin_edges = linspace(first_edge, last_edge, Nbins_, endpoint=.True.)
    h%bin_edges(Nbins_) = h%bin_edges(Nbins_) + 1.5_dp * Small
  end subroutine get_bins

  !> width_sturges
  !!
  function width_sturges(x) result(w)
    implicit none
    real(8) :: w !<
    real(8), dimension(:), intent(IN) :: x !< Array to analyze
    w = (maxval(x) - minval(x)) * log10(2._dp) / log10(size(x) + 1._dp)
  end function width_sturges!> width_sturges

  !>
  !!
  function width_rice(x) result(w)
    implicit none
    real(8) :: w !<
    real(8), dimension(:), intent(IN) :: x !< Array to analyze
    w = (maxval(x) - minval(x)) / (2._dp * size(x)**(1._dp / 3._dp))
  end function width_rice

  !> Doane's histogram bin estimator.
  !!
  !! Improved version of Sturges' formula which works better for non-normal data.
  !! @see stats.stackexchange.com/questions/55134/doanes-formula-for-histogram-binning
  function width_doane(x) result(w)
    implicit none
    real(dp) :: w !< width
    real(dp), dimension(:), intent(IN) :: x !< Array to analyze
    real(dp), dimension(size(x)) :: temp
    real(dp) :: g1
    real(dp) :: sigma
    real(dp) :: log10_2
    integer :: N
    w = (maxval(x) - minval(x))
    N = size(x)
    IF (N <= 2) return
    log10_2 = log10(2._dp)
    sigma = std(x)
    if (sigma > Zero) then
      temp = x - mean(x)
      temp = (sigma / temp)**3
      g1 = abs(mean(temp)) / sqrt(6._dp * (N - 2) / ((N + 1._dp) * (N + 3._dp)))
      w = w * log10_2 / (log10_2 + log10(real(N, kind=dp)) + log10(1._dp + g1))
    end if
  end function width_doane

end module histograms
