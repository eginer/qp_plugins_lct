BEGIN_PROVIDER [double precision, rho_alpha_hf_ao, (ao_num,ao_num)]
 implicit none 
 integer :: i,j,k,ii
 rho_alpha_hf_ao = 0.d0
 do i = 1, n_valence_orb_for_hf(1)
  ii = list_valence_orb_for_hf(i,1)
  do j = 1, ao_num
   do k = 1, ao_num
    rho_alpha_hf_ao(j,k) += mo_coef(j,ii) * mo_coef(k,ii)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, rho_alpha_hf_ao_in_r, (ao_num,n_points_final_grid) ]
 implicit none
 integer :: ipoint, j,k
 rho_alpha_hf_ao_in_r = 0.d0
 do ipoint = 1, n_points_final_grid
  do j = 1, ao_num
   do k = 1, ao_num
    rho_alpha_hf_ao_in_r(j,ipoint) += rho_alpha_hf_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, rho_beta_hf_ao, (ao_num,ao_num)]
 implicit none 
 integer :: i,j,k,ii
 rho_beta_hf_ao = 0.d0
 do i = 1, n_valence_orb_for_hf(2)
  ii = list_valence_orb_for_hf(i,2)
  do j = 1, ao_num
   do k = 1, ao_num
    rho_beta_hf_ao(j,k) += mo_coef(j,ii) * mo_coef(k,ii)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, rho_beta_hf_ao_in_r, (ao_num,n_points_final_grid) ]
 implicit none
 integer :: ipoint, j,k
 rho_beta_hf_ao_in_r = 0.d0
 do ipoint = 1, n_points_final_grid
  do j = 1, ao_num
   do k = 1, ao_num
    rho_beta_hf_ao_in_r(j,ipoint) += rho_beta_hf_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, full_dens_ao, (ao_num,ao_num)]
 implicit none
 integer :: i,j,k 
 full_dens_ao = 0.d0
 do i = 1, mo_num
  do j = 1, ao_num
   do k = 1, ao_num
    full_dens_ao(k,j) += mo_coef(k,i) * mo_coef(j,i) 
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, full_dens_ao_in_r, (ao_num,n_points_final_grid) ]
 implicit none
 integer :: ipoint, j,k
 full_dens_ao_in_r = 0.d0
 do ipoint = 1, n_points_final_grid
  do j = 1, ao_num
   do k = 1, ao_num
    full_dens_ao_in_r(j,ipoint) += full_dens_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, f_hf_ab_ao, (n_points_final_grid)]
 implicit none
 integer :: alph,bet,gam,delt,i,ipoint
 double precision :: thresh, integral 
 integer :: sze,sze_max,non_zero_int
 integer, allocatable :: out_val_index(:,:)
 double precision, allocatable :: out_val(:),alpha_dens(:)
 sze = ao_num
 sze_max = ao_num * ao_num
 thresh = thresh_int_mu_of_r
 allocate(out_val(sze_max),out_val_index(2,sze_max),alpha_dens(n_points_final_grid))

 f_hf_ab_ao = 0.d0
 !  First alpha pair of AO
 do alph = 1, ao_num !
  do gam = 1, ao_num !
   call get_ao_two_e_integrals_non_zero_jl(alph,gam,thresh,sze_max,sze,out_val,out_val_index,non_zero_int) 
   do ipoint = 1, n_points_final_grid
    alpha_dens(ipoint) = rho_alpha_hf_ao_in_r(gam,ipoint) * full_dens_ao_in_r(alph,ipoint)
   enddo
   do i = 1, non_zero_int
    bet  = out_val_index(1,i)
    delt = out_val_index(2,i) 
    integral = out_val(i)
    do ipoint = 1, n_points_final_grid
     f_hf_ab_ao(ipoint) += rho_beta_hf_ao_in_r(bet,ipoint) * full_dens_ao_in_r(delt,ipoint) * integral * alpha_dens(ipoint) 
    enddo
   enddo
  enddo
 enddo


END_PROVIDER 
