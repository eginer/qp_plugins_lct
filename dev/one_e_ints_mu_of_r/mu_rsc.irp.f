
 BEGIN_PROVIDER [double precision, average_mu_rs_c     ]
&BEGIN_PROVIDER [double precision, mu_of_r_rs_c, (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_rs_c, (3,n_points_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 include 'constants.include.F'
 double precision :: weight, rho_a_hf, rho_b_hf, rho_hf, mu_rs_c
 average_mu_rs_c    = 0.d0
 double precision :: elec_a,elec_b
 double precision :: mu_plus, mu_minus, r(3),dx,mu_lda
 dx = 1.d-5
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  rho_a_hf = one_e_dm_and_grad_alpha_in_r(4,ipoint,1)
  rho_b_hf = one_e_dm_and_grad_beta_in_r(4,ipoint,1)
  rho_hf = rho_a_hf + rho_b_hf
  mu_of_r_rs_c(ipoint,1) = mu_rs_c(rho_hf) 
  average_mu_rs_c   +=  mu_of_r_rs_c(ipoint,1) * rho_hf * weight
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

  do m = 1, 3
   r(:) = final_grid_points(:,ipoint) 
   r(m) += dx 
   call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
   rho_hf = rho_a_hf + rho_b_hf
   mu_plus = mu_rs_c(rho_hf)
   r(:) = final_grid_points(:,ipoint) 
   r(m) -= dx 
   call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
   rho_hf = rho_a_hf + rho_b_hf
   mu_minus = mu_rs_c(rho_hf)

   grad_mu_of_r_rs_c(m,ipoint,1) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo

 enddo 
 average_mu_rs_c   = average_mu_rs_c  / dble(elec_a+ elec_b)

END_PROVIDER 

BEGIN_PROVIDER [ double precision, average_mu_rs_c_lda]
 implicit none
 average_mu_rs_c_lda = 0.5d0 * (average_mu_rs_c + average_mu_lda)
END_PROVIDER 

 BEGIN_PROVIDER [double precision, average_mu_grad_n   ]
 implicit none
 integer :: ipoint,i,m
 double precision :: sqpi
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 sqpi = dsqrt(dacos(-1.d0))
 cst_rs   = (4.d0 * dacos(-1.d0)/3.d0)**(-1.d0/3.d0)
 alpha_rs = 2.d0 * dsqrt((9.d0 * dacos(-1.d0)/4.d0)**(-1.d0/3.d0)) / sqpi
 average_mu_grad_n  = 0.d0 
 double precision :: elec_a,elec_b
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  rho_a_hf = 0.d0
  grad_n   = 0.d0
  rho_a_hf = one_e_dm_and_grad_alpha_in_r(4,ipoint,1)
  rho_b_hf = one_e_dm_and_grad_beta_in_r(4,ipoint,1)
  grad_n = one_e_grad_2_dm_alpha_at_r(ipoint,1) + one_e_grad_2_dm_beta_at_r(ipoint,1)
  grad_n += 2.d0 * scal_prod_grad_one_e_dm_ab(ipoint,1)
  rho_hf = rho_a_hf + rho_b_hf
  grad_n = dsqrt(grad_n)
  if(dabs(rho_hf).gt.1.d-20)then
   grad_n = grad_n/(4.d0 * rho_hf)
  else
   grad_n = 0.d0
  endif
  rs = cst_rs * rho_hf**(-1.d0/3.d0)
  g0 = g0_UEG_mu_inf(rho_a_hf,rho_b_hf)

  average_mu_grad_n += grad_n * rho_hf * weight
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight
 enddo 
 print*,'elec_a,elec_b',elec_a,elec_b
 average_mu_grad_n = average_mu_grad_n/ dble(elec_a+ elec_b)


END_PROVIDER 

double precision function mu_rs_c(rho)
 implicit none
 double precision, intent(in) :: rho
 include 'constants.include.F'
 double precision :: cst_rs,alpha_rs,rs
 cst_rs   = (4.d0 * dacos(-1.d0)/3.d0)**(-1.d0/3.d0)
 alpha_rs = 2.d0 * dsqrt((9.d0 * dacos(-1.d0)/4.d0)**(-1.d0/3.d0)) / sqpi
  
 rs = cst_rs * rho**(-1.d0/3.d0)
 mu_rs_c =  alpha_rs/dsqrt(rs) 

end

 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_rs_c, (n_points_extra_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 include 'constants.include.F'
 double precision :: weight, rho_a_hf, rho_b_hf, rho_hf, mu_rs_c
 average_mu_rs_c    = 0.d0
 double precision :: elec_a,elec_b,r(3)
 double precision :: mu_plus, mu_minus, dx,mu_lda
 dx = 1.d-5
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_extra_final_grid
  weight = final_weight_at_r_vector_extra(ipoint)
  r(:) = final_grid_points_extra(:,ipoint)
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  rho_hf = rho_a_hf + rho_b_hf
  mu_of_r_extra_grid_rs_c(ipoint,1) = mu_rs_c(rho_hf) 
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight
 enddo 

END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_hf, (n_points_extra_final_grid) ]
 implicit none 
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals 
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,sqpi,r(3)
 PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
 print*,'providing mu_of_r_extra_grid_hf ...'
 call wall_time(wall0)
 sqpi = dsqrt(dacos(-1.d0))
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf,r) & 
 !$OMP ShARED (n_points_extra_final_grid,mu_of_r_extra_grid_hf,sqpi,final_grid_points_extra) 
 do ipoint = 1, n_points_extra_final_grid
  r(:) = final_grid_points_extra(:,ipoint)
  call f_HF_valence_ab(r,r,f_hf,on_top)
  if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
    w_hf   = 1.d+10
  else 
    w_hf  = f_hf /  on_top
  endif
  mu_of_r_extra_grid_hf(ipoint) =  w_hf * sqpi * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_extra_grid_hf = ',wall1-wall0
 END_PROVIDER 





