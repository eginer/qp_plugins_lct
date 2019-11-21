! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/pradines/program/qpt2TESt/qp2/src/rsdft_cipsi/EZFIO.cfg


BEGIN_PROVIDER [ character*(32), weight_boltz_type  ]
  implicit none
  BEGIN_DOC
! name of the exchange functional
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_rsdft_cipsi_weight_boltz_type(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: weight_boltz_type ] <<<<< ..'
      call ezfio_get_rsdft_cipsi_weight_boltz_type(weight_boltz_type)
    else
      print *, 'rsdft_cipsi/weight_boltz_type not found in EZFIO file'
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
    call MPI_BCAST( weight_boltz_type, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read weight_boltz_type with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
