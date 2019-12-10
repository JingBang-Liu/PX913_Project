!> @brief Module that writes charge density, potential, field strength
!>        particle trajectory, particle velocity, particle acceleration
!>        into a netCDF file.
!> @note 
!>       
!> @note the is greatly inspired from the
!>       module write_netcdf written by CS Brady and H Ratcliffe
MODULE netcdf_write

  USE kinds
  USE netcdf

  IMPLICIT NONE

  !!!! Define a type for run data so you don't have to type a lot
  TYPE :: run_data
    INTEGER :: run_data_N, run_data_T
    CHARACTER(LEN=25) :: run_data_init
  END TYPE

  CONTAINS

  SUBROUTINE write_electrostatics(kin_data,filename,ierr)
     
    !!!! Define output variables
    TYPE(kinematics) :: kin_data
    INTEGER, DIMENSION(:), ALLOCATABLE :: x_axis, y_axis, t_axis
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, PARAMETER :: ndims = 2
    !!!! Define dimensions of output variables
    CHARACTER(LEN=1), DIMENSION(ndims) :: dims_=(/"x", "y" /)
    CHARACTER(LEN=1), DIMENSION(1) :: dims_history=(/"t"/)
    CHARACTER(LEN=1), DIMENSION(1) :: dims_x_axis=(/"X"/)
    CHARACTER(LEN=1), DIMENSION(1) :: dims_y_axis=(/"Y"/)
    CHARACTER(LEN=1), DIMENSION(1) :: dims_t_axis=(/"T"/)
    !!!! Define sizes and dimension ids for output variables
    INTEGER, DIMENSION(ndims) :: sizes_larr, dim_ids_larr
    INTEGER :: sizes_history, dim_ids_history
    INTEGER :: dim_ids_x_axis, dim_ids_y_axis, dim_ids_t_axis 
    !!!! Define file id and variables ids for output variables
    INTEGER :: file_id, var_id_larr, var_id_history, var_id_initial_grid, i
    INTEGER :: var_id_x_axis, var_id_y_axis, var_id_t_axis
    INTEGER :: ierr
    


    sizes_larr = SHAPE(larr)
    sizes_history = SIZE(history)

    !!!! Generate axis
    ALLOCATE(x_axis(sizes_larr(1)))
    ALLOCATE(y_axis(sizes_larr(2)))
    ALLOCATE(t_axis(sizes_history))
    x_axis = 1
    y_axis = 1
    t_axis = 0
    DO i=2,sizes_larr(1)
      x_axis(i) = x_axis(i-1) + 1
    END DO
    DO i=2,sizes_larr(2)
      y_axis(i) = y_axis(i-1) + 1
    END DO
    DO i=2,sizes_history
      t_axis(i) = t_axis(i-1) + 1
    END DO
      
    
    !!!! Create the file, overwriting if it exists
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)

    !!!! Check if the file is created successfully
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Add run datas to global attribute
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"N",r_d%run_data_N)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"T",r_d%run_data_N)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"init",r_d%run_data_init)
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

