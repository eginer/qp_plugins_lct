! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home_lct/dtraore/qp2/src/finite_difference_method/EZFIO.cfg


BEGIN_PROVIDER [ double precision, field_strenght  ]
  implicit none
  BEGIN_DOC
! Electric field strenght, in atomic unit, default from HalKloHel-JCP-99
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_finite_difference_method_field_strenght(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: field_strenght ] <<<<< ..'
      call ezfio_get_finite_difference_method_field_strenght(field_strenght)
    else
      print *, 'finite_difference_method/field_strenght not found in EZFIO file'
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
    call MPI_BCAST( field_strenght, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read field_strenght with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER
