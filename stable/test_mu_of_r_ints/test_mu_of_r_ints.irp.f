program test_mu_of_r_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  mu_of_r_potential = "hf_coallescence" 
  touch mu_of_r_potential
  call routine_print
end

subroutine routine_print
 implicit none
 integer  :: i,j,k,l
 double precision :: accu,integral
 double precision, allocatable :: integrals(:)
 double precision :: get_ao_bielec_integral_erf_mu_of_r
 allocate(integrals(ao_num))
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
                                           ! 1 2 1
    call get_ao_bielec_integrals_erf_mu_of_r(i,j,k,ao_num,integrals) 
    do l = 1, ao_num
     !                                             1 2 1 2
     integral = get_ao_bielec_integral_erf_mu_of_r(i,j,k,l,ao_integrals_erf_mu_of_r_map)
     accu += dabs(integral - integrals(l))
    enddo
   enddo
  enddo
 enddo
 accu = accu / dble(ao_num **4)
 print*,'accu = ',accu
end
