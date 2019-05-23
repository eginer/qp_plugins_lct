! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/pradines/program/qp2TEST/qp2/src/mu_of_r_ints/EZFIO.cfg


BEGIN_PROVIDER [ character*(32), io_mo_integrals_sr_mu_of_r  ]
  implicit none
  BEGIN_DOC
! Read/Write MO integrals with the short range interaction with mu depending on r from/to disk [ Write | Read | None ]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_mu_of_r_ints_io_mo_integrals_sr_mu_of_r(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: io_mo_integrals_sr_mu_of_r ] <<<<< ..'
      call ezfio_get_mu_of_r_ints_io_mo_integrals_sr_mu_of_r(io_mo_integrals_sr_mu_of_r)
    else
      print *, 'mu_of_r_ints/io_mo_integrals_sr_mu_of_r not found in EZFIO file'
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
    call MPI_BCAST( io_mo_integrals_sr_mu_of_r, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read io_mo_integrals_sr_mu_of_r with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

  BEGIN_PROVIDER [ logical, read_mo_integrals_sr_mu_of_r ]
 &BEGIN_PROVIDER [ logical, write_mo_integrals_sr_mu_of_r ]

  BEGIN_DOC
  ! One level of abstraction for mo_integrals_sr_mu_of_r
  END_DOC

   if (io_mo_integrals_sr_mu_of_r.EQ.'Read') then
     read_mo_integrals_sr_mu_of_r =  .True.
     write_mo_integrals_sr_mu_of_r = .False.
   else if  (io_mo_integrals_sr_mu_of_r.EQ.'Write') then
     read_mo_integrals_sr_mu_of_r = .False.
     write_mo_integrals_sr_mu_of_r =  .True.
   else if (io_mo_integrals_sr_mu_of_r.EQ.'None') then
     read_mo_integrals_sr_mu_of_r = .False.
     write_mo_integrals_sr_mu_of_r = .False.
   else
     print *, 'io_mo_integrals_sr_mu_of_r has a bad type'
     stop 1
   endif

 END_PROVIDER

BEGIN_PROVIDER [ character*(32), io_mo_integrals_mu_of_r  ]
  implicit none
  BEGIN_DOC
! Read/Write MO integrals with the long range interaction with mu depending on r from/to disk [ Write | Read | None ]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_mu_of_r_ints_io_mo_integrals_mu_of_r(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: io_mo_integrals_mu_of_r ] <<<<< ..'
      call ezfio_get_mu_of_r_ints_io_mo_integrals_mu_of_r(io_mo_integrals_mu_of_r)
    else
      print *, 'mu_of_r_ints/io_mo_integrals_mu_of_r not found in EZFIO file'
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
    call MPI_BCAST( io_mo_integrals_mu_of_r, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read io_mo_integrals_mu_of_r with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

  BEGIN_PROVIDER [ logical, read_mo_integrals_mu_of_r ]
 &BEGIN_PROVIDER [ logical, write_mo_integrals_mu_of_r ]

  BEGIN_DOC
  ! One level of abstraction for mo_integrals_mu_of_r
  END_DOC

   if (io_mo_integrals_mu_of_r.EQ.'Read') then
     read_mo_integrals_mu_of_r =  .True.
     write_mo_integrals_mu_of_r = .False.
   else if  (io_mo_integrals_mu_of_r.EQ.'Write') then
     read_mo_integrals_mu_of_r = .False.
     write_mo_integrals_mu_of_r =  .True.
   else if (io_mo_integrals_mu_of_r.EQ.'None') then
     read_mo_integrals_mu_of_r = .False.
     write_mo_integrals_mu_of_r = .False.
   else
     print *, 'io_mo_integrals_mu_of_r has a bad type'
     stop 1
   endif

 END_PROVIDER

BEGIN_PROVIDER [ character*(32), io_ao_integrals_sr_mu_of_r  ]
  implicit none
  BEGIN_DOC
! Read/Write AO integrals with the short range interaction with mu depending on r from/to disk [ Write | Read | None ]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_mu_of_r_ints_io_ao_integrals_sr_mu_of_r(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: io_ao_integrals_sr_mu_of_r ] <<<<< ..'
      call ezfio_get_mu_of_r_ints_io_ao_integrals_sr_mu_of_r(io_ao_integrals_sr_mu_of_r)
    else
      print *, 'mu_of_r_ints/io_ao_integrals_sr_mu_of_r not found in EZFIO file'
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
    call MPI_BCAST( io_ao_integrals_sr_mu_of_r, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read io_ao_integrals_sr_mu_of_r with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

  BEGIN_PROVIDER [ logical, read_ao_integrals_sr_mu_of_r ]
 &BEGIN_PROVIDER [ logical, write_ao_integrals_sr_mu_of_r ]

  BEGIN_DOC
  ! One level of abstraction for ao_integrals_sr_mu_of_r
  END_DOC

   if (io_ao_integrals_sr_mu_of_r.EQ.'Read') then
     read_ao_integrals_sr_mu_of_r =  .True.
     write_ao_integrals_sr_mu_of_r = .False.
   else if  (io_ao_integrals_sr_mu_of_r.EQ.'Write') then
     read_ao_integrals_sr_mu_of_r = .False.
     write_ao_integrals_sr_mu_of_r =  .True.
   else if (io_ao_integrals_sr_mu_of_r.EQ.'None') then
     read_ao_integrals_sr_mu_of_r = .False.
     write_ao_integrals_sr_mu_of_r = .False.
   else
     print *, 'io_ao_integrals_sr_mu_of_r has a bad type'
     stop 1
   endif

 END_PROVIDER

BEGIN_PROVIDER [ character*(32), io_ao_integrals_mu_of_r  ]
  implicit none
  BEGIN_DOC
! Read/Write AO integrals with the long range interaction with mu depending on r from/to disk [ Write | Read | None ]
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_mu_of_r_ints_io_ao_integrals_mu_of_r(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: io_ao_integrals_mu_of_r ] <<<<< ..'
      call ezfio_get_mu_of_r_ints_io_ao_integrals_mu_of_r(io_ao_integrals_mu_of_r)
    else
      print *, 'mu_of_r_ints/io_ao_integrals_mu_of_r not found in EZFIO file'
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
    call MPI_BCAST( io_ao_integrals_mu_of_r, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read io_ao_integrals_mu_of_r with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

  BEGIN_PROVIDER [ logical, read_ao_integrals_mu_of_r ]
 &BEGIN_PROVIDER [ logical, write_ao_integrals_mu_of_r ]

  BEGIN_DOC
  ! One level of abstraction for ao_integrals_mu_of_r
  END_DOC

   if (io_ao_integrals_mu_of_r.EQ.'Read') then
     read_ao_integrals_mu_of_r =  .True.
     write_ao_integrals_mu_of_r = .False.
   else if  (io_ao_integrals_mu_of_r.EQ.'Write') then
     read_ao_integrals_mu_of_r = .False.
     write_ao_integrals_mu_of_r =  .True.
   else if (io_ao_integrals_mu_of_r.EQ.'None') then
     read_ao_integrals_mu_of_r = .False.
     write_ao_integrals_mu_of_r = .False.
   else
     print *, 'io_ao_integrals_mu_of_r has a bad type'
     stop 1
   endif

 END_PROVIDER
