! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/eginer/programs/qp2/src/mu_of_r/EZFIO.cfg


BEGIN_PROVIDER [ double precision, mu_of_r_array , (n_points_final_grid) ]
  implicit none
  BEGIN_DOC
! array of the values of mu(r)
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    if (size(mu_of_r_array) == 0) return

    call ezfio_has_mu_of_r_mu_of_r_array(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: mu_of_r_array ] <<<<< ..'
      call ezfio_get_mu_of_r_mu_of_r_array(mu_of_r_array)
    else
      print *, 'mu_of_r/mu_of_r_array not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mu_of_r_array, (n_points_final_grid), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mu_of_r_array with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ character*(32), mu_of_r_potential  ]
  implicit none
  BEGIN_DOC
! type of potential for the mu(r) interaction: can be [hf_coallescence | psi_coallescence | Read]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_mu_of_r_mu_of_r_potential(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: mu_of_r_potential ] <<<<< ..'
      call ezfio_get_mu_of_r_mu_of_r_potential(mu_of_r_potential)
    else
      print *, 'mu_of_r/mu_of_r_potential not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mu_of_r_potential, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mu_of_r_potential with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
