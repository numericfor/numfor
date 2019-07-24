program test
  use basic, only: dp, Zero, timer
  use strings, only: center
  use grids, only: linspace, geomspace, mean, std
  USE histograms
  implicit none

  real(dp), dimension(:), allocatable :: b
  type(histog) :: hist
  integer, parameter ::   Nbins = 50
  real(dp), dimension(Nbins) :: bins
  integer :: Ndat

  integer :: u
  integer :: j
  type(timer) :: T

  Ndat = 1000 * 100000

  print *, Ndat < huge(1), huge(1) / Ndat

  allocate (b(Ndat))
  b = geomspace(1.e-30_dp, 1.0_dp, num=Ndat)
  b = b + linspace(0._dp, 1.0_dp, num=Ndat)

  bins = linspace(minval(b), maxval(b), Nbins)

  call T%start()
  ! hist = histogram(b, Nbins=Nbins, range=[0.15_dp, 1.85_dp], density=.True.)
  hist = histogram(b, bins=bins, density=.True.)
  call T%stop()

  Ndat = size(hist%hist)

  ! Print the bins
  ! print '(4(g0.5,1x))', hist%bin_edges

  open (newunit=u, file='data/histog1.dat')
  write (u, '(2(g0.5,1x))') (0.5_dp * (hist%bin_edges(j) + hist%bin_edges(j + 1)), hist%hist(j), j=1, Ndat)
  close (u)
  print '(A)', repeat('-', 65)
  print '(A)', center(' Test histogram 1 ', 65, '*')
  print '(A)', 'Result in file: data/histog1.dat'

  print '(A)', repeat('-', 65)
  print '(A)', center(' Statistics: ', 65, '*')
  print '(A)', center(repeat('-', len('Statistics:')), 65, ' ')
  print *, 'Number of points:', size(b)
  print *, 'Number of bins:', size(hist%hist)
  print *, 'Number of terms in range:', hist%n
  print *, 'Mean value:', mean(hist%hist)
  print *, 'Standard deviation:', std(hist%hist)
  print *, 'Integral:', sum(hist%hist * (hist%bin_edges(2:) - hist%bin_edges(1:Ndat - 1)))
  print '(A)', center(' Timing: ', 35, '')
  call T%show()

end program test

