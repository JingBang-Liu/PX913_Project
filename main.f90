! File name: main.f90
! Author: David Liu
! Notes: Example of verlet integration use.


PROGRAM main
  USE kinds
  USE domain_tools
  USE command_line
  USE verlet_integrator
  USE gauss_seidel
  IMPLICIT NONE

  REAL(KIND=REAL64), DIMENSION(2) :: init_pos
  REAL(KIND=REAL64), DIMENSION(2) :: init_vel
  REAL(KIND=REAL64) :: dx, dy, dt
  INTEGER :: time   
  INTEGER :: i, j
  TYPE(kinematics) :: kin_data

  REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: x, y
  REAL(KIND=REAL64), DIMENSION(2) :: xrange=(/-1,1/) 
  REAL(KIND=REAL64), DIMENSION(2) :: yrange=(/-1,1/)
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: rho, phi, Ex, Ey 
  CHARACTER(len=20) :: problem 
  LOGICAL nx_succ, ny_succ, problem_succ
 
  ! first fill in the type that stores all the data.
  ! test = the type you're filling in
  ! time = controls the length of the array.
  
  ! To use the verlet integrator, create a derived type >test<
  ! Pass it into the subroutine with:
  ! Ex, Ey = electric fields
  ! init_pos = a vector containing the initial positions
  ! init_vel = a vector containing the initial velocities
  ! dx, dy = grid discretization
  ! dt = time step
  ! time = number of time steps

  CALL parse_arg() 
  nx_succ=get_arg('nx', nx)
  ny_succ=get_arg('ny', ny)
  peoblem_succ=get_arg('problem', problem)

  if (nx_succ .AND. ny_succ .AND. problem_succ)
    print*, "Arguments Passed"
  else
    print*, "Arguments Missing"
    return
  end if
  
  CALL create_axis(x, nx, xrange, 1)
  CALL create_axis(y, ny, yrange, 1)

  dx=real(xrange(2)-xrange(1))/(real(nx)-1.0)
  dy=real(yrange(2)-yrange(1))/(real(ny)-1.0)

  CALL charge(rho, nx ,ny, problem, x, y)
  CALL potential(phi, nx, ny, rho, dx, dy)
  
  dt = 0.01
  time = 1000
     
  CALL fillType(kin_data, time)
  
  ! time has to be allocated before fillType (see initializations)
  CALL verlet(kin_data, Ex, Ey, init_pos, init_vel, dx, dy, dt, time)

  ! access the values like this:
  ! kin_data%pos_history(1, i), kin_data%vel_history(1, i), kin_data%acc_history(1, i)

  ! NetCDF section
  
     
  ! deallocates the arrays associated with the calledtype.
  CALL emptyType(kin_data)
  
END PROGRAM main
