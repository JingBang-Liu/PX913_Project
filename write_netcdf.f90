!> @brief Module that writes charge density, potential, field strength
!>        particle trajectory, particle velocity, particle acceleration
!>        into a netCDF file.
!> @Author: JingBang Liu 
!>       
!> @note the is greatly inspired from the
!>       module write_netcdf written by CS Brady and H Ratcliffe
MODULE netcdf_write

  USE kinds
  USE netcdf
  USE verlet_integrator

  IMPLICIT NONE

  !!!! Define a type for run data so you don't have to type a lot
  TYPE :: run_data
    INTEGER :: run_data_nx, run_data_ny
    CHARACTER(LEN=25) :: run_data_problem
  END TYPE

  CONTAINS

  SUBROUTINE write_electrostatics(rho,phi,Ex,Ey,kin_data,r_d,filename,ierr)

     
    !!!! Define input variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: rho, phi, Ex, Ey
    TYPE(kinematics), INTENT(IN) :: kin_data
    TYPE(run_data) :: r_d
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: ierr
    !!! Parameters define dimensions
    !INTEGER, PARAMETER :: ndims1 = 1
    INTEGER, PARAMETER :: ndims2 = 2
    !!! Define dimensions for output variables
    CHARACTER(LEN=5), DIMENSION(ndims2) :: dims_rho=(/"rho_x","rho_y"/)
    CHARACTER(LEN=5), DIMENSION(ndims2) :: dims_phi=(/"phi_x","phi_y"/)
    CHARACTER(LEN=4), DIMENSION(ndims2) :: dims_Ex=(/"Ex_x","Ex_y"/)
    CHARACTER(LEN=4), DIMENSION(ndims2) :: dims_Ey=(/"Ey_x","Ey_y"/)
    CHARACTER(LEN=8), DIMENSION(ndims2) :: dims_pos=(/"pos_grid","pos_t   "/)
    CHARACTER(LEN=8), DIMENSION(ndims2) :: dims_vel=(/"vel_grid","vel_t   "/)
    CHARACTER(LEN=8), DIMENSION(ndims2) :: dims_acc=(/"acc_grid","acc_t   "/)
    !CHARACTER(LEN=1), DIMENSION(ndims1) :: dims_nx=(/"x"/)
    !CHARACTER(LEN=1), DIMENSION(ndims1) :: dims_ny=(/"y"/)
    !!! Define the size ids for output variables
    INTEGER, DIMENSION(ndims2) :: sizes_rho, sizes_phi, sizes_Ex, sizes_Ey
    INTEGER, DIMENSION(ndims2) :: sizes_pos, sizes_vel, sizes_acc
    !!! Define the dimension ids for output variables
    INTEGER, DIMENSION(ndims2) :: dim_ids_rho, dim_ids_phi, dim_ids_Ex, dim_ids_Ey
    INTEGER, DIMENSION(ndims2) :: dim_ids_pos, dim_ids_vel, dim_ids_acc
    !!! Define variable ids for output variables
    INTEGER :: var_id_rho, var_id_phi, var_id_Ex, var_id_Ey, var_id_pos, var_id_vel, var_id_acc
    !!! Define file id
    INTEGER :: file_id 
    INTEGER :: i

    !!! Get sizes of variables
    sizes_rho = SHAPE(rho)
    sizes_phi = SHAPE(phi)
    sizes_Ex = SHAPE(Ex)
    sizes_Ey = SHAPE(Ey)
    sizes_pos = SHAPE(kin_data%pos_history)
    sizes_vel = SHAPE(kin_data%vel_history)
    sizes_acc = SHAPE(kin_data%acc_history)
    
    !!!! Create the file, overwriting if it exists
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)

    !!!! Check if the file is created successfully
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Add run datas to global attribute
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"nx",r_d%run_data_nx)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"ny",r_d%run_data_ny)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"problem",r_d%run_data_problem)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for rho, check for any error
    DO i=1,ndims2
      ierr = nf90_def_dim(file_id,dims_rho(i),sizes_rho(i),dim_ids_rho(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id, "rho_charge_density", NF90_DOUBLE, dim_ids_rho, var_id_rho)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    
    !!!! Define dim id and var id for phi, check for any error
    DO i=1,ndims2
      ierr = nf90_def_dim(file_id,dims_phi(i),sizes_phi(i),dim_ids_phi(i))
      IF (ierr /=nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id, "phi_electric_potential", NF90_DOUBLE, dim_ids_phi, var_id_phi)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for Ex, check for any error
    DO i=1,ndims2
      ierr = nf90_def_dim(file_id,dims_Ex(i),sizes_Ex(i),dim_ids_Ex(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id,"Ex_field_intensity",NF90_DOUBLE,dim_ids_Ex,var_id_Ex)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for Ey, check for any error
    DO i=1,ndims2
      ierr = nf90_def_dim(file_id,dims_Ey(i),sizes_Ey(i),dim_ids_Ey(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id,"Ey_field_intensity",NF90_DOUBLE,dim_ids_Ey,var_id_Ey)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for particle postision, check for any error
    DO i=1,ndims2
      ierr = nf90_def_dim(file_id,dims_pos(i),sizes_pos(i),dim_ids_pos(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id,"particle_position",NF90_DOUBLE,dim_ids_pos,var_id_pos)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    
    !!!! Define dim id and var id for particle velocity, check for any error
    DO i=1,ndims2
      ierr = nf90_def_dim(file_id,dims_vel(i),sizes_vel(i),dim_ids_vel(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id,"particle_velocity",NF90_DOUBLE,dim_ids_vel,var_id_vel)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    
    !!!! Define dim id and var id for particle acceleration, check for any error
    DO i=1,ndims2
      ierr = nf90_def_dim(file_id,dims_acc(i),sizes_acc(i),dim_ids_acc(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id,"particle_acceleration",NF90_DOUBLE,dim_ids_acc,var_id_acc)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Finish defining metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    
    !!!! Write rho into file
    ierr = nf90_put_var(file_id, var_id_rho, rho) 
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "rho"
      RETURN
    END IF

    !!!! Write phi into file
    ierr = nf90_put_var(file_id, var_id_phi, phi)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "phi"
      RETURN
    END IF

    !!!! Write Ex
    ierr = nf90_put_var(file_id,var_id_Ex, Ex)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "Ex"
      RETURN
    END IF

    !!!! Write Ey
    ierr = nf90_put_var(file_id,var_id_Ey, Ey)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "Ey"
      RETURN
    END IF

    !!!! Write particle position
    ierr = nf90_put_var(file_id,var_id_pos, kin_data%pos_history)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "pos"
      RETURN
    END IF
    
    !!!! Write particle velocity
    ierr = nf90_put_var(file_id,var_id_vel, kin_data%vel_history)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "vel"
      RETURN
    END IF

    !!!! Write particle acceleration
    ierr = nf90_put_var(file_id,var_id_acc, kin_data%acc_history)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "acc"
      RETURN
    END IF


    !!!! Close the file
    ierr = nf90_close(file_id)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END SUBROUTINE write_electrostatics

END MODULE netcdf_write

