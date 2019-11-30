  ! File name: verlet_integrator.f90
  ! Author: Aravinthen Rajkumar
  ! Notes: This file contains a Verlet Integration routine for particle motion.

MODULE verlet_integrator
  USE kinds
  IMPLICIT NONE
  
  TYPE kinematics
     ! global type
     ! This stores all of the kinematic data of a run.
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: pos_history ! position history
     REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: vel_history ! velocity history
     REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: acc_history ! acceleration history          
  END type kinematics
  
CONTAINS

  FUNCTION pos2cell(position, delta)
    ! delta refers to the size of the spatial discretization step.
    ! varies by dimension: dx and dy. assignment document states that we'll be given these.
    ! see comments in the "verlet" function regarding grid 
    REAL(KIND=REAL64), INTENT(IN) :: position ! position to be converted to cell
    REAL(KIND=REAL64), INTENT(IN) :: delta ! position to be converted to cell
    INTEGER :: pos2cell    
    pos2cell = FLOOR((position - 1.0_REAL64)/delta) + 1
  END FUNCTION pos2cell
      
  FUNCTION verlet(Ex, Ey, init_pos, init_vel, dx, dy, dt, time)
    ! Input:
    !    Ex, Ey - the electic fields
    !    init_pos, init_vel - the initial positions and velocities (in a vector format)
    !    dx, dy, dt - the grid partition lengths along with the magnitude of the time step
    !    time - number of time steps
    
!!!!
!!!! We need to figure out the best way to organise the grid sizes between the
!!!! verlet integrator and the gauss-seidel method.
!!!! Right now, the grid sizes are being hard-coded into the function.
!!!!    
    TYPE(kinematics) :: verlet    
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: Ex, Ey ! field matrix
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(2) :: init_pos ! initial position
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(2) :: init_vel ! initial position
    REAL(KIND=REAL64) :: dt, dx, dy ! time-step value and spatial step values.
    
    REAL(KIND=REAL64) :: pos_y, pos_x ! positions
    REAL(KIND=REAL64) :: vel_y, vel_x ! velocties
    REAL(KIND=REAL64) :: acc_y, acc_x ! accelerations
    REAL(KIND=REAL64) :: old_acc_y, old_acc_x ! history variables for accelerations
    
    INTEGER :: i, time ! number of time-steps and loop variables

    ! the initial data.
    pos_x = init_pos(1)
    pos_y = init_pos(2)
    vel_x = init_vel(1)
    vel_y = init_vel(2)
    ! this bit uses the pos2cel function to access the electric fields.
    acc_x = Ex(pos2cell(init_pos(1), dx), pos2cell(init_pos(2), dy))
    acc_y = Ey(pos2cell(init_pos(1), dx), pos2cell(init_pos(2), dy))

    DO i = 1, time
       ! store the accelerations
       old_acc_x = acc_x
       old_acc_y = acc_y

       ! calculate new positions
       pos_x = pos_x + vel_x*dt + (0.5)*(acc_x)*(acc_x)*dt*dt
       pos_y = pos_y + vel_y*dt + (0.5)*(acc_y)*(acc_y)*dt*dt
       ! accelerations have to be calculated before velocites.
       acc_x = Ex(pos2cell(pos_x, dx), pos2cell(pos_y, dy))
       acc_y = Ey(pos2cell(pos_y, dx), pos2cell(pos_y, dy))
       ! make use of old_acc here
       vel_x = vel_x + dt*(0.5*(acc_x + old_acc_x))
       vel_y = vel_y + dt*(0.5*(acc_y + old_acc_y))
       
       print*, pos_x, pos_y, vel_x, vel_y, acc_x, acc_y

    END DO
    
  END FUNCTION verlet
  
END MODULE verlet_integrator


PROGRAM main
  USE kinds
  USE verlet_integrator
  IMPLICIT NONE
  
  REAL(KIND=REAL64), DIMENSION(10, 10) :: Ex, Ey
  REAL(KIND=REAL64), DIMENSION(2) :: init_pos
  REAL(KIND=REAL64), DIMENSION(2) :: init_vel
  REAL(KIND=REAL64) :: dx, dy, dt
  INTEGER :: time = 10
  TYPE(kinematics) :: test  

  Ex = 0.0
  Ey = 0.0
  init_pos(1) = 0.0
  init_pos(2) = 0.0
  init_vel(1) = 20
  init_vel(2) = 20
  dx = 0.1
  dy = 0.1
  dt = 0.01
  
  test = verlet(Ex, Ey, init_pos, init_vel, dx, dy, dt, time)


END PROGRAM main
