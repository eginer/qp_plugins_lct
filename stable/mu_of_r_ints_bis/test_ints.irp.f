program test_ints
 implicit none
 call test_ao
 call test_mo
 call test_mo_sym
 call test_ao_sym
 call save_erf_mu_of_r_ao_into_erf_ints_ao
 call save_erf_mu_of_r_mo_into_erf_ints_mo
! call save_erf_two_e_integrals_ao
! call save_erf_two_e_integrals_mo
end

subroutine test_mo
 implicit none
 integer :: i,j,k,l
 double precision :: get_mo_two_e_int_erf_mu_of_r,map_int,prov_int,accu,contrib
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     map_int =  get_mo_two_e_int_erf_mu_of_r(i,j,k,l,mo_int_erf_mu_of_r_map)
     prov_int = mo_mu_of_r_two_ints_phys(i,j,k,l)
     contrib = dabs(map_int  - prov_int)
!     print*,map_int,prov_int,contrib
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu MO = ',accu/dble(mo_num**4)

end

subroutine test_ao
 implicit none
 integer :: i,j,k,l
 double precision :: get_ao_two_e_int_erf_mu_of_r,map_int,prov_int,accu,contrib
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     map_int =  get_ao_two_e_int_erf_mu_of_r(i,j,k,l,ao_int_erf_mu_of_r_map)
     prov_int = ao_mu_of_r_two_ints_phys(i,j,k,l)
     contrib = dabs(map_int  - prov_int)
!     print*,map_int,prov_int,contrib
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu AO = ',accu/dble(ao_num**4)

end

subroutine test_mo_sym
 implicit none
 integer :: i,j,k,l
 double precision :: get_mo_two_e_int_erf_mu_of_r,int_direct,int_sym,accu,contrib
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     int_direct =  get_mo_two_e_int_erf_mu_of_r(i,j,k,l,mo_int_erf_mu_of_r_map)
     int_sym    =  get_mo_two_e_int_erf_mu_of_r(j,i,l,k,mo_int_erf_mu_of_r_map)
     contrib = dabs(int_direct - int_sym)
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num**4)
end

subroutine test_ao_sym
 implicit none
 integer :: i,j,k,l
 double precision :: get_ao_two_e_int_erf_mu_of_r,int_direct,int_sym,accu,contrib
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     int_direct =  get_ao_two_e_int_erf_mu_of_r(i,j,k,l,ao_int_erf_mu_of_r_map)
     int_sym    =  get_ao_two_e_int_erf_mu_of_r(j,i,l,k,ao_int_erf_mu_of_r_map)
     contrib = dabs(int_direct - int_sym)
     accu += contrib
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(ao_num**4)
end
