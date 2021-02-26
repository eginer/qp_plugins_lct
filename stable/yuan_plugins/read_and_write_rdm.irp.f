BEGIN_PROVIDER [ integer*8, dim_two_rdm  ]
  implicit none
  BEGIN_DOC
! number of elements in two-rdm
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_yuan_plugins_dim_two_rdm(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: dim_two_rdm ] <<<<< ..'
      call ezfio_get_yuan_plugins_dim_two_rdm(dim_two_rdm)
    else
      print *, 'yuan_plugins/dim_two_rdm not found in EZFIO file'
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
    call MPI_BCAST( dim_two_rdm, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read dim_two_rdm with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

subroutine read_two_rdm_and_write_to_ezfio(n,n_mo_tmp)
 implicit none
 integer, intent(in) :: n_mo_tmp
 integer*8, intent(in) :: n

 integer :: i,j,k,l
 integer*8 :: m
 double precision :: value_rdm
 double precision, allocatable :: two_rdm(:,:,:,:,:)
 character*(128) :: name_file 
 integer :: nstates
 nstates = 1

 allocate(two_rdm(n_mo_tmp,n_mo_tmp,n_mo_tmp,n_mo_tmp,nstates))
 character*(1) :: coma
 open(1, file = 'two_rdm') 
 two_rdm = 0.d0
 do m = 1_8, n
   ! a^{l}a^{k} a_i a_j
!  read(1,'(4(I3,A1),F16.13)')l,coma, k, coma, j, coma, i,coma,value_rdm
  read(1,*)l, k,  j,  i,value_rdm
  two_rdm(l,k,j,i,1) = 0.5d0 * value_rdm
 enddo
 close(1)

 name_file = 'act_2_rdm_ab_mo'
 ! Writting the two rdm on the alpha/beta
 call write_array_two_rdm(n_mo_tmp,nstates,two_rdm,name_file)
 call ezfio_set_two_body_rdm_io_two_body_rdm_ab("Read")

 ! Writting the two rdm on the spin trace
 name_file = 'act_2_rdm_spin_trace_mo'
 call write_array_two_rdm(n_mo_tmp,nstates,two_rdm,name_file)
 call ezfio_set_two_body_rdm_io_two_body_rdm_spin_trace("Read")

 deallocate(two_rdm)
end


subroutine read_one_rdm_sp_tr_and_write_to_ezfio(n_mo_tmp)
 implicit none
 integer, intent(in) :: n_mo_tmp

 double precision, allocatable :: one_rdm(:,:),one_rdm_full(:,:)
 allocate(one_rdm(n_mo_tmp,n_mo_tmp),one_rdm_full(mo_num, mo_num))
 integer :: k,l,m,n
 integer :: ncore
 ncore = mo_num - n_mo_tmp
 character*(1) :: coma
 open(1, file = 'one_rdm') 
 do l = 1, n_mo_tmp
  read(1,*)one_rdm(l,1:n_mo_tmp)
 enddo
 close(1)
 one_rdm_full = 0.d0
 do k = 1, ncore
  one_rdm_full(k,k) = 2.d0
 enddo
 m = 0
 do k = ncore+1, ncore + n_mo_tmp
  m += 1
  write(*,'(100(F16.10,X))')one_rdm(m,:)
  n = 0
  do l = ncore+1 , ncore + n_mo_tmp
   n += 1
   one_rdm_full(l,k) = one_rdm(n,m)
  enddo
 enddo
 one_rdm_full = one_rdm_full * 0.5d0

 ! Writting the one rdm on the alpha/beta
 call ezfio_set_aux_quantities_data_one_e_dm_alpha_mo(one_rdm_full)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_mo(one_rdm_full)
 call ezfio_set_density_for_dft_density_for_dft("input_density")

 deallocate(one_rdm)
end
