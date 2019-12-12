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

  IMPLICIT NONE

  !!!! Define a type for run data so you don't have to type a lot
  TYPE :: run_data
    INTEGER :: run_data_nx, run_data_ny
    CHARACTER(LEN=*) :: run_data_init
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
    INTEGER, PARAMETER :: ndims1 = 1
    INTEGER, PARAMETER :: ndims2 = 2
    !!! Define dimensions for output variables
    CHARACTER(LEN=*), DIMENSION(ndims2) :: dims_time=(/"axis","t"/)
    CHARACTER(LEN=1), DIMENSION(ndims2) :: dims_grid=(/"x","y"/)
    CHARACTER(LEN=1), DIMENSION(ndims1) :: dims_nx=(/"x"/)
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
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"nx",r_d%run_data_N)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"ny",r_d%run_data_N)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"problem",r_d%run_data_init)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for larr, check for any error
    DO i=1,ndims
      ierr = nf90_def_dim(file_id,dims_larr(i),sizes_larr(i),dim_ids_larr(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ierr = nf90_def_var(file_id, "cell_status", NF90_INT, dim_ids_larr, var_id_larr)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    
    !!!! Define dim id and var id for history, check for any error
    ierr = nf90_def_dim(file_id,dims_history(1),sizes_history,dim_ids_history)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "time_history", NF90_DOUBLE, dim_ids_history, var_id_history)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for x_axis, check for any error
    ierr = nf90_def_dim(file_id,dims_x_axis(1),sizes_larr(1),dim_ids_x_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id,"x_axis",NF90_INT,dim_ids_x_axis,var_id_x_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for y_axis, check for any error
    ierr = nf90_def_dim(file_id,dims_y_axis(1),sizes_larr(2),dim_ids_y_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id,"y_axis",NF90_INT,dim_ids_y_axis,var_id_y_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim id and var id for t_axis, check for any error
    ierr = nf90_def_dim(file_id,dims_t_axis(1),sizes_history,dim_ids_t_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id,"t_axis",NF90_INT,dim_ids_t_axis,var_id_t_axis)
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
    
    !!!! Write larr into file
    ierr = nf90_put_var(file_id, var_id_larr, larr) 
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Write history into file
    ierr = nf90_put_var(file_id, var_id_history, history)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Write x_axis
    ierr = nf90_put_var(file_id,var_id_x_axis, x_axis)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_var(file_id,var_id_y_axis, y_axis)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_var(file_id,var_id_t_axis, t_axis)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Close the file
    ierr = nf90_close(file_id)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    DEALLOCATE(x_axis)
    DEALLOCATE(y_axis)
    DEALLOCATE(t_axis)

  END SUBROUTINE write_electrostatics

END MODULE netcdf_write

