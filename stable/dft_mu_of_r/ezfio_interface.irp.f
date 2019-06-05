! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/eginer/qp2/src/dft_mu_of_r/EZFIO.cfg


BEGIN_PROVIDER [ character*(32), mu_of_r_functional  ]
  implicit none
  BEGIN_DOC
! name of the correlation functional for the mu of r basis set correction
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dft_mu_of_r_mu_of_r_functional(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: mu_of_r_functional ] <<<<< ..'
      call ezfio_get_dft_mu_of_r_mu_of_r_functional(mu_of_r_functional)
    else
      print *, 'dft_mu_of_r/mu_of_r_functional not found in EZFIO file'
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
    call MPI_BCAST( mu_of_r_functional, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mu_of_r_functional with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ double precision, ci_threshold  ]
  implicit none
  BEGIN_DOC
! Threshold on the CI coefficients as in QMCChem
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_dft_mu_of_r_ci_threshold(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: ci_threshold ] <<<<< ..'
      call ezfio_get_dft_mu_of_r_ci_threshold(ci_threshold)
    else
      print *, 'dft_mu_of_r/ci_threshold not found in EZFIO file'
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
    call MPI_BCAST( ci_threshold, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ci_threshold with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
