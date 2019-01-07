 BEGIN_PROVIDER [double precision, big_array_coulomb_integrals_erf_mu_of_r, (mo_tot_num,mo_tot_num, mo_tot_num)]
&BEGIN_PROVIDER [double precision, big_array_exchange_integrals_erf_mu_of_r,(mo_tot_num,mo_tot_num, mo_tot_num)]
 implicit none
 integer :: i,j,k,l
 double precision :: get_mo_bielec_integral_erf_mu_of_r
 double precision :: integral

 do k = 1, mo_tot_num
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
     l = j
     integral = get_mo_bielec_integral_erf_mu_of_r(i,j,k,l,mo_integrals_erf_mu_of_r_map)
     big_array_coulomb_integrals_erf_mu_of_r(j,i,k) = integral
     l = j
     integral = get_mo_bielec_integral_erf_mu_of_r(i,j,l,k,mo_integrals_erf_mu_of_r_map)
     big_array_exchange_integrals_erf_mu_of_r(j,i,k) = integral
   enddo
  enddo
 enddo


END_PROVIDER 


!BEGIN_PROVIDER [double precision, big_array_coulomb_integrals_sr, (mo_tot_num,mo_tot_num, mo_tot_num)]
!&BEGIN_PROVIDER [double precision, big_array_exchange_integrals_sr,(mo_tot_num,mo_tot_num, mo_tot_num)]
!implicit none
!integer :: i,j,k,l
!double precision :: get_mo_bielec_integral_sr
!double precision :: integral

!do k = 1, mo_tot_num
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num
!    l = j
!    integral = get_mo_bielec_integral_sr(i,j,k,l,mo_integrals_sr_map)
!    big_array_coulomb_integrals_sr(j,i,k) = integral
!    l = j
!    integral = get_mo_bielec_integral_sr(i,j,l,k,mo_integrals_sr_map)
!    big_array_exchange_integrals_sr(j,i,k) = integral
!  enddo
! enddo
!enddo


!ND_PROVIDER 

