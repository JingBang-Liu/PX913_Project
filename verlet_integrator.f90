  ! File name: verlet_integrator.f90
  ! Author: Aravinthen Rajkumar
  ! Notes: This file contains a Verlet Integration routine for particle motion.

MODULE verlet_integrator
  USE kinds
  IMPLICIT NONE
  
  TYPE kinematics
     ! global type
     ! This stores all of the kinematic data of a run.
     REAL(KIND=REAL64), ALLOCATABLE :: pos_history(:,:) ! position history
     REAL(KIND=REAL64), ALLOCATABLE :: vel_history(:,:) ! velocity history
     REAL(KIND=REAL64), ALLOCATABLE :: acc_history(:,:) ! acceleration history
  END type kinematics
  
CONTAINS

  SUBROUTINE fillType(data_list, T)
    ! fills the kinematics data arrays before use.
    TYPE(kinematics) :: data_list
    INTEGER :: T

    ALLOCATE(data_list%pos_history(2,T))
    ALLOCATE(data_list%vel_history(2,T))
    ALLOCATE(data_list%acc_history(2,T))        
  END SUBROUTINE fillType

  SUBROUTINE emptyType(data_list)
    ! empty the kinematics data arrays.
    TYPE(kinematics) :: data_list
    DEALLOCATE(data_list%pos_history)
    DEALLOCATE(data_list%vel_history)
    DEALLOCATE(data_list%acc_history)
  END SUBROUTINE emptyType

  FUNCTION pos2cell(position, delta)
    ! delta refers to the size of the spatial discretization step.
    ! varies by dimension: dx and dy. assignment document states that we'll be given these.
    ! see comments in the "verlet" function regarding grid

    ! Note! this gives the negative values of the grid!
    ! Quick hack: add the dimensions of Ex and Ey.
    REAL(KIND=REAL64), INTENT(IN) :: position ! position to be converted to cell
    REAL(KIND=REAL64), INTENT(IN) :: delta ! position to be converted to cell
    INTEGER :: pos2cell    
    pos2cell = FLOOR((position - 1.0_REAL64)/delta)
  END FUNCTION pos2cell
      
  SUBROUTINE Verlet(data_entry, Ex, Ey, init_pos, init_vel, dx, dy, dt, time)
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
    TYPE(kinematics) :: data_entry
    
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: Ex, Ey ! field matrix
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(2) :: init_pos ! initial position
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(2) :: init_vel ! initial position
    REAL(KIND=REAL64) :: dt, dx, dy ! time-step value and spatial step values.
    
    REAL(KIND=REAL64) :: pos_y, pos_x ! positions
    REAL(KIND=REAL64) :: vel_y, vel_x ! velocties
    REAL(KIND=REAL64) :: acc_y, acc_x ! accelerations
    REAL(KIND=REAL64) :: old_acc_y, old_acc_x ! history variables for accelerations
    
    INTEGER :: i, time, cellx, celly ! number of time-steps and loop variables
    INTEGER, DIMENSION(2) :: dims
    dims = SHAPE(Ex)
   
    ! the initial data.
    pos_x = init_pos(1)
    pos_y = init_pos(2)
    
    vel_x = init_vel(1)
    vel_y = init_vel(2)
    ! this bit uses the pos2cel function to access the electric fields.
    acc_x = Ex(pos2cell(pos_x, dx), pos2cell(pos_y, dy))
    acc_y = Ey(pos2cell(pos_x, dx), pos2cell(pos_y, dy)) 

    DO i = 1, time
       ! data input into the type, data_entry: this is the technical output.
       data_entry%pos_history(1,i) = pos_x
       data_entry%pos_history(2,i) = pos_y
       data_entry%vel_history(1,i) = vel_x
       data_entry%vel_history(2,i) = vel_y
       data_entry%acc_history(1,i) = acc_x
       data_entry%acc_history(2,i) = acc_y      
       
       ! store the accelerations
       old_acc_x = acc_x
       old_acc_y = acc_y

       ! calculate new positions
       pos_x = pos_x + vel_x*dt + (0.5)*(acc_x)*(acc_x)*dt*dt
       pos_y = pos_y + vel_y*dt + (0.5)*(acc_y)*(acc_y)*dt*dt

       ! You have to multiply the dimensions 
       cellx = pos2cell(pos_x, dx) + dims(1)
       celly = pos2cell(pos_y, dy) + dims(2)
       
       ! conditional doesn't allow us to pick values of the array that are larger than the array size.
       IF ((cellx.GT.(dims(1)-1)).OR.(celly.GT.(dims(2)-1))) THEN
          acc_x = 0.0
          acc_y = 0.0
          vel_x = 0.0
          vel_y = 0.0
       ELSE
          ! accelerations have to be calculated before velocites.
          acc_x = Ex(cellx, celly)
          acc_y = Ey(cellx, celly)
          ! make use of old_acc here
          vel_x = vel_x + dt*(0.5*(acc_x + old_acc_x))
          vel_y = vel_y + dt*(0.5*(acc_y + old_acc_y))
       END IF             
    END DO
    
  END SUBROUTINE Verlet
  
END MODULE verlet_integrator
