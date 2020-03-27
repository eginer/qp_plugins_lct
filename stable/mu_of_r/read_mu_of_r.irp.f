
BEGIN_PROVIDER [double precision, mu_of_r_read , (n_points_final_grid,N_states)]
 implicit none
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists
  PROVIDE ezfio_filename

  if (mpi_master) then
    ! Coefs
    call ezfio_has_mu_of_r_mu_of_r_read(exists)
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_coef with MPI'
    endif
  IRP_ENDIF

  if (exists) then
    if (mpi_master) then
      call ezfio_get_mu_of_r_mu_of_r_read(mu_of_r_read)
      write(*,*) 'Read  mu_of_r_read'
    endif
    IRP_IF MPI
      call MPI_BCAST( mu_of_r_read, n_points_final_grid, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_coef with MPI'
      endif
    IRP_ENDIF
  else
   print*,'You need mu_of_r_read but there is no file in the EZFIO data base '
   stop
  endif

END_PROVIDER 
