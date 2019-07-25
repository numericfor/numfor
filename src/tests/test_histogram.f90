program test
  use basic, only: dp, Zero, Small, timer
  use strings, only: center
  use grids, only: linspace, geomspace, mean, std
  USE histograms
  implicit none

  real(dp), dimension(:), allocatable :: b
  real(dp), dimension(:), allocatable :: w
  type(histog) :: hist, h
  integer, parameter ::   Nbins = 99
  real(dp), dimension(Nbins) :: bins
  integer :: Ndat

  integer :: u
  integer :: j
  type(timer) :: T

  Ndat = 10000 * 10000
  print *, 'huge/Ndat', huge(1) / Ndat

  allocate (b(Ndat))
  b = linspace(0._dp, 1.0_dp, num=Ndat)
  b = b + geomspace(1.e-30_dp, 1.0_dp, num=Ndat)

  w = 5._dp / (1.0_dp + b)        ! Weights
  ! print '(A)', center('data', 30, '-')
  ! print '(5(g0.3,1x))', b

  bins = linspace(minval(b), maxval(b), Nbins)

  call T%start()
  ! hist = histogram(b, Nbins=Nbins, range=[0.15_dp, 1.85_dp], density=.True.)
  hist = histogram(b, bins=bins, weights=w, density=.True.)
  call T%stop()

  print '(A)', repeat('-', 65)
  print "(A)", center("Example for numpy documentation", 65, '-')
  h = histogram(linspace(0.5_dp, 5.5_dp, 6), Nbins=5)
  print "(A, 6(1x,g0.2))", 'Values', linspace(0.5_dp, 5.5_dp, 6)
  print "(A, 6(1x,g0.2))", 'Bin Edges:', h%bin_edges
  print "(A, 5(1x,g0.2))", 'Histogram:', h%hist

  print '(A)', repeat('-', 65)
  print "(A)", center(" Example for numpy documentation ", 65, '-')
  h = histogram([1._dp, 2._dp, 1._dp], bins=[0._dp, 1._dp, 2._dp, 3._dp])
  print "(A, 3(1x,g0.2))", 'Values:', [1._dp, 2._dp, 1._dp]
  print "(A, 4(1x,g0.2))", 'Bin Edges:', h%bin_edges
  print "(A, 3(1x,g0.2))", 'Histogram:', h%hist

  open (newunit=u, file='data/histog1.dat')
  write (u, '(3(g0.5,1x))') (hist%bin_edges(j), hist%bin_edges(j + 1), hist%hist(j), j=1, size(hist%hist))
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
  print *, 'Integral:', sum(hist%hist * (hist%bin_edges(2:) - hist%bin_edges(1:size(hist%hist) - 1)))
  print '(A)', center(' Timing: ', 35, '')
  call T%show()

end program test

