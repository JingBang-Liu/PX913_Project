module gauss_seidel

use kinds

implicit none

contains

subroutine charge(rho, nx, ny, problem, x, y)

  real(real64), dimension(:,:) :: rho
  real(real64), dimension(:), allocatable :: x, y 
  integer :: i,j, nx, ny
  character(len=20) :: problem
   
  print*, problem=='single' 
    if(problem=='null') then
      rho=0

    else if(problem=='single') then
      do i=1, nx
        do j=1, ny
          rho(i,j)=EXP(-(x(i)/0.1_Real64)**2-(y(j)/0.1_real64)**2)
          print*, x(i), y(j), rho(i,j)
        end do
      end do

    else if(problem=='double') then
      do i=1, nx
        do j=1, ny
          rho(i,j)=EXP(-((x(i)+0.25)/0.1)**2-((y(j)+0.25)/0.1)**2)+EXP(-((x(i)-0.75)/0.2)**2-((y(j)-0.75)/0.2)**2)
        end do
      end do
  end if
end subroutine

subroutine potential(phi, nx, ny, rho, dx, dy)

  real(real64), dimension(:, :), allocatable :: phi
  real(real64), dimension(:,:) :: rho
  real(real64) :: e_tot, d_sum, ratio=1, dx, dy
  integer :: nx, ny, i, j, iter=0

  phi=0  
  iter=0
  do i=0, nx+1
    do j=0, ny+1
      !print*, phi(i,j)
    end do
  end do
  !print*, size(phi), lbound(phi), ubound(phi)
  do while (ratio>1e-5) 
    !Update the grid
    do i=1, nx
      do j=1, ny
        phi(i,j)=-(rho(i,j)-(phi(i+1,j) + phi(i-1,j))/dx**2-(phi(i,j+1)+phi(i,j-1))/dy**2 )/(2_real64/dx**2+2_real64/dy**2)
        !print*, rho(i,j), phi(i,j)
      end do
    end do
    !check error
    e_tot=0
    d_sum=0
    do i=1, nx
      do j=1, ny
        e_tot=e_tot+ABS((phi(i-1,j)-2*phi(i,j)+phi(i+1,j))/dx**2+(phi(i,j-1)-2*phi(i,j)+phi(i,j+1))/dy**2-rho(i,j))
        
        d_sum=d_sum+((phi(i-1,j)-2*phi(i,j)+phi(i+1,j))/dx**2+(phi(i,j-1)-2*phi(i,j)+phi(i,j+1))/dy**2)**2
        
      end do
    end do
    ratio=e_tot/sqrt(d_sum)
    print*, e_tot, d_sum, ratio
    iter=iter+1
  end do
  print*, iter, e_tot, d_sum, ratio
end subroutine

subroutine field(Ex, Ey , phi, dx, dy)
  real(real64), dimension(:,:), allocatable :: phi, E
  real(real64) :: dx, dy
  
  do i=1, nx
    do j=1, ny
      Ex(i,j)=(phi(i+1,j)-phi(i-1,j))/(2*dx)
      Ey(i,j)=(phi(i,j+1)-phi(i,j-1))/(2*dy)
    end do
  end do
end subroutine
  
end module

program main

use domain_tools
use gauss_seidel

implicit none

integer, parameter :: nx=5, ny=5
real(REAL64), dimension(:), allocatable :: x, y
real(REAL64), dimension(2) :: xrange=(/-1,1/) 
real(REAL64), dimension(2) :: yrange=(/-1,1/)
real(real64), dimension(:,:), allocatable :: rho, phi, E 
real(real64) :: dx, dy
integer :: i, j
character(len=20) :: problem='single'

allocate(rho(nx, ny))
allocate(phi(0:nx+1, 0:ny+1))
allocate(Ex(1:nx, 1:ny))
allocate(Ey(1:nx, 1:ny))

call create_axis(x, nx, xrange, 1)
call create_axis(y, ny, yrange, 1)

dx=real(xrange(2)-xrange(1))/(real(nx)-1.0)
dy=real(yrange(2)-yrange(1))/(real(ny)-1.0)

rho=0
phi=0

call charge(rho, nx ,ny, problem, x, y)
print*, "-------------space------------"
call potential(phi, nx, ny, rho, dx, dy)


do i=0, nx+1
  do j=0, ny+1
!    print*,x(i), y(j), phi(i,j)
  end do
end do

end program


