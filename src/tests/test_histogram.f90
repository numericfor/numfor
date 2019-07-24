program test
  use basic, only: dp, Zero, Small, timer
  use strings, only: center
  use grids, only: linspace, geomspace, mean, std
  USE histograms
  implicit none

  real(dp), dimension(:), allocatable :: b
  type(histog) :: hist, h
  integer, parameter ::   Nbins = 51
  real(dp), dimension(Nbins) :: bins
  integer :: Ndat

  integer :: u
  integer :: j
  type(timer) :: T

  Ndat = 10 * 10000

  allocate (b(Ndat))
  b = linspace(0._dp, 1.0_dp, num=Ndat)
  b = b + geomspace(1.e-30_dp, 1.0_dp, num=Ndat)

  ! print '(A)', center('data', 30, '-')
  ! print '(5(g0.3,1x))', b

  bins = linspace(minval(b), maxval(b), Nbins)
  print "(A,2(g0,1x))", 'min,max=', minval(b), maxval(b)
  ! print '(A)', center('bins', 30, '-')
  ! print '(6(g0.3,1x))', bins

  call T%start()
  ! hist = histogram(b, Nbins=Nbins, range=[0.15_dp, 1.85_dp], density=.True.)
  hist = histogram(b, bins=bins, density=.False.)
  call T%stop()

  Ndat = size(hist%hist)
  print *, "Menor?", minval(b) == hist%bin_edges(1), minval(b) - hist%bin_edges(1)
  print *, "Mayor?", maxval(b) == hist%bin_edges(Nbins), maxval(b) - hist%bin_edges(Nbins)

  h = histogram(linspace(0.5_dp, 5.5_dp, 6), Nbins=5)
  print "(6(g0.2,1x))", linspace(0.5_dp, 5.5_dp, 6)
  print "(6(g0.2,1x))", h%bin_edges
  print "(5(g0.2,1x))", h%hist
  ! Print the bins
  ! print '(4(g0.5,1x))', hist%bin_edges

  open (newunit=u, file='data/histog1.dat')
  write (u, '(3(g0.5,1x))') (hist%bin_edges(j), hist%bin_edges(j + 1), hist%hist(j), j=1, Ndat)
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
  print *, 'Mean value:', mean(b)
  print *, 'Standard deviation:', std(b)
  print *, 'Integral:', sum(hist%hist * (hist%bin_edges(2:) - hist%bin_edges(1:Ndat - 1)))
  print '(A)', center(' Timing: ', 35, '')
  call T%show()

end program test

