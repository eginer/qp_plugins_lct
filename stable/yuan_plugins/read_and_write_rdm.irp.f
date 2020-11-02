
subroutine read_two_rdm_and_write_to_ezfio(n,n_mo_tmp)
 implicit none
 integer, intent(in) :: n,n_mo_tmp

 double precision, allocatable :: two_rdm(:,:,:,:)
 allocate(two_rdm(n_mo_tmp,n_mo_tmp,n_mo_tmp,n_mo_tmp))
 integer :: i,j,k,l,m
 double precision :: value_rdm
 character*(1) :: coma
 open(1, file = 'two_rdm') 
 do m = 1, n
  read(1,'(4(I3,A1),F16.13)')l,coma, k, coma, j, coma, i,coma,value_rdm
  two_rdm(l,k,j,i) = value_rdm
 enddo
 close(1)

 ! Writting the two rdm on the alpha/beta
 call ezfio_set_two_body_rdm_two_rdm_ab_disk(two_rdm)
 call ezfio_set_two_body_rdm_io_two_body_rdm_ab("Read")

 ! Writting the two rdm on the spin trace
 call ezfio_set_two_body_rdm_two_rdm_spin_trace_disk(two_rdm)
 call ezfio_set_two_body_rdm_io_two_body_rdm_spin_trace("Read")

 deallocate(two_rdm)
end
