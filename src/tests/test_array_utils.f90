!> test
program test
  USE basic
  USE grids
  USE array_utils
  implicit none
  real(dp), dimension(21) :: x, y
  real(dp), dimension(2, 21) :: z
  character(len=:), allocatable :: filename

  x = linspace(0._dp, 10._dp, 21, endpoint=.True.)
  y = -x**2 / 10._dp
  z(1, :) = x
  z(2, :) = y
  call save_array(x)  ! One array in one column to stdout
  filename = "data/output1.dat"
  call save_array(x, fname=filename) ! One array in one column to file
  filename = "data/output2.dat"
  call save_array([x, y], 2, filename) ! Two arrays in two columns to file
  filename = "data/output3.dat"
  call save_array([x, y], 2, filename, fmt="f0.8") ! Changing output format
  filename = "data/output4.dat"
  ! 2D arrays
  call save_array(z, fname=filename, fmt="(f0.8,1x,es15.8)") ! Different formats for each column
  filename = "data/output5.dat"
  call save_array(z, fname=filename, delimiter=',', fmt="f0.8") ! Save with a comma as delimiter
  filename = "data/output6.dat"
  call save_array(z, fname=filename, delimiter=',', header="A nice header", fmt="g0.8") ! Add header
  filename = "data/output7.dat"
  call save_array(z, fname=filename, delimiter=',', header="A nice header"&
    &, footer="And a nice footer", fmt="g0.8") ! Add footer
  filename = "data/output8.dat"
  call save_array(z, fname=filename, delimiter=',', comment="!", header="A nice header"&
    &, footer="And a nice footer", fmt="g0.8") !  Change comment char

end program test
