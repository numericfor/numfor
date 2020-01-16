program test_rand_expon
  use numfor, only: dp, Zero, str, random_seed, random_exponential, ran_exponential_pdf, histog, histogram,&
    & save_array, mean, timer
  implicit none
  integer :: i
  integer, parameter :: Ndim = 10**6
  real(dp), dimension(Ndim) :: samples
  real(dp), dimension(80) :: x, y
  type(histog) :: h
  type(timer) :: T
  real(dp), dimension(1) :: times

  times = 0._dp
  call random_seed()     ! Uses a random seed
  do i = 1, size(times)
    call T%start()
    call random_exponential(mu=1._dp, x=samples)
    call T%stop()
    times(i) = T%elapsed
  end do
  print "(A)", str(size(times))//" loops. "//str(1000 * mean(times))//" ms for each loop"
  h = histogram(samples, Nbins=size(x), range=[-7._dp, 10._dp], density=.True.)

  x = (h%bin_edges(:size(h%bin_edges) - 1) + h%bin_edges(2:)) / 2
  y = ran_exponential_pdf(1._dp, x)

  call save_array([x, y, h%hist], 3, 'data/rand_exponential_hist.dat')
  ! Outputs:

end program test_rand_expon
