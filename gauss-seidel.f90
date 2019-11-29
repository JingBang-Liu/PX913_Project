program main

use kinds
use domain_tools

implicit none

integer, parameter :: nx=10, ny=10
real(REAL64), dimension(:), allocatable :: x, y
real(REAL64), dimension(2) :: xrange=(/-1,1/) 
real(REAL64), dimension(2) :: yrange=(/-1,1/)

call create_axis(x, nx, xrange, 1)
call create_axis(y, ny, yrange, 1)

print*, x
print*, y

end program


