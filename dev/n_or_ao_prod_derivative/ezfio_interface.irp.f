! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/pradines/program/qpt2TESt/qp2/src/n_or_ao_prod_derivative/EZFIO.cfg


BEGIN_PROVIDER [ integer, order_derivative_ontop  ]
  implicit none
  BEGIN_DOC
! Order of the expension of the spherical average on top pair density
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_n_or_ao_prod_derivative_order_derivative_ontop(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: order_derivative_ontop ] <<<<< ..'
      call ezfio_get_n_or_ao_prod_derivative_order_derivative_ontop(order_derivative_ontop)
    else
      print *, 'n_or_ao_prod_derivative/order_derivative_ontop not found in EZFIO file'
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
    call MPI_BCAST( order_derivative_ontop, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read order_derivative_ontop with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
