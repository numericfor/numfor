program test
  use basic, only: dp, Zero, timer
  use strings, only: center
  use grids, only: linspace, geomspace, mean, std
  USE histograms
  implicit none

  real(dp), dimension(:), allocatable :: b
  type(histog) :: hist
  real(dp), dimension(50) :: bins
  integer :: Ndat
  integer :: u
  integer :: j
  type(timer) :: T

  bins = linspace(Zero, 2.5_dp, 50)

  Ndat = 9000000
  allocate (b(Ndat))
  b = geomspace(1.e-30_dp, 1.0_dp, num=Ndat)
  b = b + linspace(0._dp, 1.0_dp, num=Ndat)

  call T%start()
  hist = histogram(b, Nbins=50, range=[0.15_dp, 1.85_dp], density=.True.)
  call T%stop()

  Ndat = size(hist%bin_edges(1, :))

  open (newunit=u, file='data/histog1.dat')
  write (u, '(2(g0.5,1x))') (0.5_dp * (hist%bin_edges(1, j) + hist%bin_edges(2, j)), hist%hist(j), j=1, Ndat)
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
  print *, 'Integral:', sum(hist%hist) * (hist%bin_edges(2, 1) - hist%bin_edges(1, 1))
  print '(A)', center(' Timing: ', 35, '')
  call T%show()

end program test

