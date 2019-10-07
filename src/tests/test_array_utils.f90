!> test
program test
  USE basic
  USE grids
  USE array_utils
  implicit none
  real(dp), dimension(21) :: x, y
  character(len=:), allocatable :: filename

  x = linspace(0._dp, 10._dp, 21, endpoint=.True.)
  y = -x**2 / 10._dp
  call save_array(x)  ! One array in one column to stdout
  filename = "data/output1.dat"
  call save_array(x, fname=filename) ! One array in one column to file
  filename = "data/output2.dat"
  call save_array([x, y], 2, filename) ! Two arrays in two columns to file
  filename = "data/output3.dat"
  call save_array([x, y], 2, filename, fmt="f0.8") ! Changing output format
  filename = "data/output4.dat"
  call save_array([x, y], 2, filename, fmt="(f0.8,1x,es15.8)") ! Different formats for each column

end program test
