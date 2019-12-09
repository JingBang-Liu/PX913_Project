! File name: verlet_example.f90
! Author: Aravinthen Rajkumar
! Notes: Example of verlet integration use.


PROGRAM main
  USE kinds
  USE verlet_integrator
  
  IMPLICIT NONE
  
  REAL(KIND=REAL64), DIMENSION(10, 10) :: Ex, Ey
  REAL(KIND=REAL64), DIMENSION(2) :: init_pos
  REAL(KIND=REAL64), DIMENSION(2) :: init_vel
  REAL(KIND=REAL64) :: dx, dy, dt
  INTEGER :: time = 20
  INTEGER :: i
  TYPE(kinematics) :: test


 
  ! first fill in the type that stores all the data.
  ! test = the type you're filling in
  ! time = controls the length of the array.
  CALL fillType(test, time)
  
  ! To use the verlet integrator, create a derived type >test<
  ! Pass it into the subroutine with:
  ! Ex, Ey = electric fields
  ! init_pos = a vector containing the initial positions
  ! init_vel = a vector containing the initial velocities
  ! dx, dy = grid discretization
  ! dt = time step
  ! time = number of time steps
  
  Ex = 0.0
  Ey = 0.0
  init_pos(1) = 0.0
  init_pos(2) = 0.0
  init_vel(1) = 10
  init_vel(2) = 10
  dx = 0.1
  dy = 0.1
  dt = 0.01
  ! time has to be allocated before fillType (see initializations)
  CALL verlet(test, Ex, Ey, init_pos, init_vel, dx, dy, dt, time)


  ! access the values like this:
  DO i = 1, time
     PRINT*, test%pos_history(1, i), test%vel_history(1, i), test%acc_history(1, i)
     PRINT*, test%pos_history(2, i), test%vel_history(2, i), test%acc_history(2, i)
     PRINT*, ""
  END DO
     
  ! deallocates the arrays associated with the calledtype.
  CALL emptyType(test)
  
END PROGRAM main
