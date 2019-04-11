program test_mu_of_r_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  mu_of_r_potential = "hf_coallescence" 
  touch mu_of_r_potential
  call routine_print
  call routine_print_sym
end

subroutine routine_print
 implicit none
 use map_module
 integer  :: i,j,k,l
 double precision :: accu,integral
 double precision, allocatable :: integrals(:)
 double precision :: get_ao_bielec_integral_erf_mu_of_r
 integer :: ii(4),jj(4),kk(4),ll(4)
 integer(key_kind) :: i1,i2
 allocate(integrals(ao_num))
 accu = 0.d0
 do i = 1,  1
  do j = 1,  1
   do k = 1,  ao_num
                                           ! 1 2 1
!   call get_ao_bielec_integrals_erf_mu_of_r(i,j,k,ao_num,integrals) 
    do l = 1, ao_num
     !                                             1 2 1 2
!    integral = get_ao_bielec_integral_erf_mu_of_r(i,j,k,l,ao_integrals_erf_mu_of_r_map)
     
     call index_two_e_no_sym(i,j,k,l,ao_num,i1) 
     call index_reverse_two_e_no_sym(ii,jj,kk,ll,ao_num,i1)
     write(33,*)''
     write(33,*)i1
     write(33,*)i,j,k,l
     write(33,*)ii(1),jj(1),kk(1),ll(1)
     write(33,*)ii(2),jj(2),kk(2),ll(2)
     write(33,*)ii(3),jj(3),kk(3),ll(3)
     write(33,*)ii(4),jj(4),kk(4),ll(4)
!    stop
!    accu += dabs(integral - integrals(l))
    enddo
   enddo
  enddo
 enddo
 accu = accu / dble(ao_num **4)
 write(33,*)'accu = ',accu
end


subroutine routine_print_sym
 implicit none
 use map_module
 integer  :: i,j,k,l
 double precision :: accu,integral
 double precision, allocatable :: integrals(:)
 double precision :: get_ao_bielec_integral_erf_mu_of_r
 integer :: ii(8),jj(8),kk(8),ll(8)
 integer(key_kind) :: i1,i2
 allocate(integrals(ao_num))
 accu = 0.d0
 do i = 1, 1
  do j = 1, 1
   do k = 1, ao_num
                                           ! 1 2 1
!   call get_ao_bielec_integrals_erf_mu_of_r(i,j,k,ao_num,integrals) 
    do l = 1, ao_num
     !                                             1 2 1 2
!    integral = get_ao_bielec_integral_erf_mu_of_r(i,j,k,l,ao_integrals_erf_mu_of_r_map)
     
     call two_e_integrals_index(i,j,k,l,i1) 
     call two_e_integrals_index_reverse(ii,jj,kk,ll,i1)
     write(34,*)''
     write(34,*)i,j,k,l
     write(34,*)ii(1),jj(1),kk(1),ll(1)
     write(34,*)ii(2),jj(2),kk(2),ll(2)
     write(34,*)ii(3),jj(3),kk(3),ll(3)
     write(34,*)ii(4),jj(4),kk(4),ll(4)
!    accu += dabs(integral - integrals(l))
    enddo
   enddo
  enddo
 enddo
 accu = accu / dble(ao_num **4)
 write(34,*)'accu = ',accu
end

