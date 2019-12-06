module gauss_seidel

use kinds

implicit none

contains

subroutine charge(rho, nx, ny, problem)

  real(real64), dimension(:,:) :: rho 
  integer :: i,j, nx, ny
  character(len=20) :: problem
  
  print*, problem  
  
  select case (problem)
    case('null')
      rho=0

    case('single')
      do i=1, nx
        do j=1, ny
          rho(i,j)=EXP(-(i/0.1)**2-(j/0.1)**2)
        end do
      end do

    case('double')
      do i=1, nx
        do j=1, ny
          rho(i,j)=EXP(-((i+0.25)/0.1)**2-((j+0.25)/0.1)**2)+EXP(-((i-0.75)/0.2)**2-((j-0.75)/0.2)**2)
        end do
      end do
  end select
end subroutine

subroutine potential(phi, nx, ny, rho, dx, dy)

  real(real64), dimension(:,:) :: phi, rho
  real(real64) :: e_tot, d_sum, ratio, dx, dy
  integer :: nx, ny, i, j


  phi=0  
  do while (ratio>1e-5) 
    !Update the grid
    do i=1, nx
      do j=1, ny
        phi(i,j)=-(rho(i,j)-(phi(i+1,j) + phi(i-1,j))/dx**2-(phi(i,j+1)+phi(i,j-1))/dy**2 )/(2/dx**2+2/dy**2)
      end do
    end do
    !check error
    do i=1, nx
      do j=1, ny
        e_tot=e_tot+ABS((phi(i-1,j)-2*phi(i,j)+phi(i+1,j))/dx**2+(phi(i,j-1)-2*phi(i,j)+phi(i,j+1))/dy**2-rho(i,j))
        
        
        d_sum=d_sum+ABS((phi(i-1,j)-2*phi(i,j)+phi(i+1,j))/dx**2+(phi(i,j-1)-2*phi(i,j)+phi(i,j+1))/dy**2-rho(i,j))
        
      end do
    end do
    ratio=e_tot/sqrt(d_sum)

  end do
  
end subroutine
  
end module

program main

use domain_tools
use gauss_seidel

implicit none

integer, parameter :: nx=1, ny=1
real(REAL64), dimension(:), allocatable :: x, y
real(REAL64), dimension(2) :: xrange=(/-1,1/) 
real(REAL64), dimension(2) :: yrange=(/-1,1/)
real(real64), dimension(:,:), allocatable :: rho, phi, E 
real(real64) :: dx, dy
character(len=20) :: problem='single'

allocate(rho(0:nx+1, 0:ny+1))
allocate(phi(0:nx+1, 0:ny+1))
allocate(E(0:nx+1, 0:ny+1))

call create_axis(x, nx, xrange, 1)
call create_axis(y, ny, yrange, 1)

dx=real(xrange(2)-xrange(1))/(real(nx)-1.0)
dy=real(yrange(2)-yrange(1))/(real(ny)-1.0)

rho=0
phi=0

call charge(rho, nx ,ny, problem)
call potential(phi, nx, ny, rho, dx, dy)

print*, rho
print*, dx,dy
print*, phi

end program


