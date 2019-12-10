module gauss_seidel

use kinds

implicit none

contains

!Fills the charge density rho, depending on the problem type
!'null' gives an empty grid
!single gives a single peak at the origin
!double give two peaks at (-0.25,-0.25) and (0.75, 0.75) 
subroutine charge(rho, nx, ny, problem, x, y)

  real(real64), dimension(:,:) :: rho
  real(real64), dimension(:), allocatable :: x, y 
  integer :: i,j, nx, ny
  character(len=20) :: problem
    
  if(problem=='null') then
    rho=0

  else if(problem=='single') then
    do i=1, nx
      do j=1, ny
        rho(i,j)=EXP(-(x(i)/0.1_Real64)**2-(y(j)/0.1_real64)**2)
        !print*, x(i), y(j), rho(i,j)
      end do
    end do
 
  else if(problem=='double') then
    do i=1, nx
      do j=1, ny
        rho(i,j)=EXP(-((x(i)+0.25)/0.1)**2-((y(j)+0.25)/0.1)**2)+EXP(-((x(i)-0.75)/0.2)**2-((y(j)-0.75)/0.2)**2)
      end do
    end do
    
  else
    rho=0
  end if
end subroutine

!Fills the potential phi use gauss-seidel iteration
!Iteration is done until the ratio of the total absolute error and the rms value
!of the potential is less than 1e-5. 
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
    !print*, e_tot, d_sum, ratio
    iter=iter+1
  end do
  !print*, iter, e_tot, d_sum, ratio
end subroutine

!Fills the electric field componants, Ex and Ey by numerically differentiating phi
subroutine field(Ex, Ey , phi, dx, dy, nx, ny)
  real(real64), dimension(:,:), allocatable :: phi, Ex, Ey
  real(real64) :: dx, dy
  integer :: nx, ny, i, j
  Ex=0
  Ey=0
  do i=1, nx
    do j=1, ny
      Ex(i,j)=(phi(i+1,j)-phi(i-1,j))/(2*dx)
      Ey(i,j)=(phi(i,j+1)-phi(i,j-1))/(2*dy)
    end do
  end do
end subroutine
  
end module


