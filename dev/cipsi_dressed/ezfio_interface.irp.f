! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/emmanuel/qp2/src/cipsi_dressed/EZFIO.cfg


BEGIN_PROVIDER [ logical, save_wf_after_selection  ]
  implicit none
  BEGIN_DOC
! If true, saves the wave function after the selection, before the diagonalization
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cipsi_dressed_save_wf_after_selection(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: save_wf_after_selection ] <<<<< ..'
      call ezfio_get_cipsi_dressed_save_wf_after_selection(save_wf_after_selection)
    else
      print *, 'cipsi_dressed/save_wf_after_selection not found in EZFIO file'
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
    call MPI_BCAST( save_wf_after_selection, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read save_wf_after_selection with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, seniority_max  ]
  implicit none
  BEGIN_DOC
! Maximum number of allowed open shells. Using -1 selects all determinants
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cipsi_dressed_seniority_max(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: seniority_max ] <<<<< ..'
      call ezfio_get_cipsi_dressed_seniority_max(seniority_max)
    else
      print *, 'cipsi_dressed/seniority_max not found in EZFIO file'
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
    call MPI_BCAST( seniority_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read seniority_max with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, excitation_ref  ]
  implicit none
  BEGIN_DOC
! 1: Hartree-Fock determinant, 2:All determinants of the dominant configuration
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cipsi_dressed_excitation_ref(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: excitation_ref ] <<<<< ..'
      call ezfio_get_cipsi_dressed_excitation_ref(excitation_ref)
    else
      print *, 'cipsi_dressed/excitation_ref not found in EZFIO file'
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
    call MPI_BCAST( excitation_ref, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read excitation_ref with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, excitation_max  ]
  implicit none
  BEGIN_DOC
! Maximum number of excitation with respect to the Hartree-Fock determinant. Using -1 selects all determinants
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cipsi_dressed_excitation_max(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: excitation_max ] <<<<< ..'
      call ezfio_get_cipsi_dressed_excitation_max(excitation_max)
    else
      print *, 'cipsi_dressed/excitation_max not found in EZFIO file'
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
    call MPI_BCAST( excitation_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read excitation_max with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, excitation_alpha_max  ]
  implicit none
  BEGIN_DOC
! Maximum number of excitation for alpha determinants with respect to the Hartree-Fock determinant. Using -1 selects all determinants
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cipsi_dressed_excitation_alpha_max(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: excitation_alpha_max ] <<<<< ..'
      call ezfio_get_cipsi_dressed_excitation_alpha_max(excitation_alpha_max)
    else
      print *, 'cipsi_dressed/excitation_alpha_max not found in EZFIO file'
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
    call MPI_BCAST( excitation_alpha_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read excitation_alpha_max with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, excitation_beta_max  ]
  implicit none
  BEGIN_DOC
! Maximum number of excitation for beta determinants with respect to the Hartree-Fock determinant. Using -1 selects all determinants
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cipsi_dressed_excitation_beta_max(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: excitation_beta_max ] <<<<< ..'
      call ezfio_get_cipsi_dressed_excitation_beta_max(excitation_beta_max)
    else
      print *, 'cipsi_dressed/excitation_beta_max not found in EZFIO file'
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
    call MPI_BCAST( excitation_beta_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read excitation_beta_max with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
