! File name: main.f90
! Author: David Liu
! Notes: Example of verlet integration use.


PROGRAM main
  USE kinds
  USE verlet_integrator
  
  IMPLICIT NONE

  REAL(KIND=REAL64), DIMENSION(10, 10) :: Ex, Ey
  REAL(KIND=REAL64), DIMENSION(2) :: init_pos
  REAL(KIND=REAL64), DIMENSION(2) :: init_vel
  REAL(KIND=REAL64) :: dx, dy, dt
  INTEGER :: time   
  INTEGER :: i
  TYPE(kinematics) :: kin_data


 
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

!!! PLACEHOLDER - REPLACE WHEN GAUSS-SEIDEL IS COMPLETED !!!
  Ex = 0.0
  Ey = 0.0
  init_pos(1) = 0.0
  init_pos(2) = 0.0
  init_vel(1) = 10
  init_vel(2) = 10
  dx = 0.1
  dy = 0.1
  dt = 0.01
  time = 1000
   
!!! PLACEHOLDER - REPLACE WHEN GAUSS-SEIDEL IS COMPLETED !!!
  
  CALL fillType(kin_data, time)
  
  ! time has to be allocated before fillType (see initializations)
  CALL verlet(kin_data, Ex, Ey, init_pos, init_vel, dx, dy, dt, time)

  ! access the values like this:
  ! kin_data%pos_history(1, i), kin_data%vel_history(1, i), kin_data%acc_history(1, i)

  ! NetCDF section
  
     
  ! deallocates the arrays associated with the calledtype.
  CALL emptyType(kin_data)
  
END PROGRAM main
