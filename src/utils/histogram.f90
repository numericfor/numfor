module histogram

  !> type histogram holds the data from an histogram
  type, public :: histogram
    integer :: n = 0            !< Total number of values parsed
    real(dp), dimension(:), pointer :: hist !< Values of the histogram
    real(dp), dimension(2, :), pointer :: bin_edges !< Array (2xlen(hist)) the edges
  end type histogram

  !> histogram
  !!
  function histogram(a, Nbins, bins, range, weights, density) result(h)
    implicit none
    type(histogram) :: h !<
    real(dp), intent(IN) :: a !<
    integer, optional, intent(IN) :: Nbins !<
    real(dp), dimension(*), optional, intent(IN) :: bins !<
    real(dp), dimension(*), optional, intent(IN) :: bins !<
    real(dp), dimension(*), optional, intent(IN) :: weights !<
    logical, optional, intent(IN) :: density !<
    ! real(dp), pointer
    real(dp) :: first_edge, last_edge

    if (Present(bins)) then
      if (size(shape(bins)) == 1) then
        ! Si tiene rango 1:
        IF (bins(:size(bins) - 1) > bins(1:)) &
          & call print_msg('bins must increase monotonically', 'histogram')
        allocate (h%bin_edges(2, size(bins) - 1))
        h%bin_edges(1) = bins(:size(bins) - 1)
        h%bin_edges(2) = bins(1:)
      else if (Present(Nbins)) then
        mina = min(a); maxa = max(a)
        h%bin_edges[1] = linspace(mina, maxa, Nbins, endpoint=.False.)
        h%bin_edges[2, :Nbins - 1] = h%bin_edges[1, 1:]
        h%bin_edges[2, Nbins] = maxa
      else
        h%bin_edges = get_auto_bins()
      end if

    end function histogram
  end module histogram
