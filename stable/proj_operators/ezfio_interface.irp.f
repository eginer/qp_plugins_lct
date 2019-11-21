! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/pradines/program/qpt2TESt/qp2/src/proj_operators/EZFIO.cfg


BEGIN_PROVIDER [ double precision, thresh_pairs_mu_of_r  ]
  implicit none
  BEGIN_DOC
! Threshold for pairs of AO for the computation of mu_of_r on the AO basis
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_proj_operators_thresh_pairs_mu_of_r(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: thresh_pairs_mu_of_r ] <<<<< ..'
      call ezfio_get_proj_operators_thresh_pairs_mu_of_r(thresh_pairs_mu_of_r)
    else
      print *, 'proj_operators/thresh_pairs_mu_of_r not found in EZFIO file'
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
    call MPI_BCAST( thresh_pairs_mu_of_r, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read thresh_pairs_mu_of_r with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ double precision, thresh_int_mu_of_r  ]
  implicit none
  BEGIN_DOC
! Threshold for the two-e integral in the computation of mu_of_r on the AO basis
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_proj_operators_thresh_int_mu_of_r(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: thresh_int_mu_of_r ] <<<<< ..'
      call ezfio_get_proj_operators_thresh_int_mu_of_r(thresh_int_mu_of_r)
    else
      print *, 'proj_operators/thresh_int_mu_of_r not found in EZFIO file'
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
    call MPI_BCAST( thresh_int_mu_of_r, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read thresh_int_mu_of_r with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ character*(32), mu_of_r_potential  ]
  implicit none
  BEGIN_DOC
! type of potential for the mu(r) interaction: can be [hf_coallescence | psi_coallescence | psi_cas_ful | psi_cas_truncated | Read]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_proj_operators_mu_of_r_potential(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: mu_of_r_potential ] <<<<< ..'
      call ezfio_get_proj_operators_mu_of_r_potential(mu_of_r_potential)
    else
      print *, 'proj_operators/mu_of_r_potential not found in EZFIO file'
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
