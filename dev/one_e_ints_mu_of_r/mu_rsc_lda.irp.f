 BEGIN_PROVIDER [double precision, average_mu_lda      ]
&BEGIN_PROVIDER [double precision, mu_of_r_lda , (n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_lda , (3,n_points_final_grid,N_states) ]
 implicit none
 integer :: ipoint,i,m
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 double precision :: elec_a,elec_b
 double precision :: mu_plus, mu_minus, r(3),dx,mu_lda
 dx = 1.d-4
 average_mu_lda     = 0.d0
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  rho_a_hf = one_e_dm_and_grad_alpha_in_r(4,ipoint,1)
  rho_b_hf = one_e_dm_and_grad_beta_in_r(4,ipoint,1)
  mu_of_r_lda(ipoint,1) = mu_lda(rho_a_hf,rho_b_hf)

  average_mu_lda    +=  mu_of_r_lda(ipoint,1) * weight * rho_hf
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight

  do m = 1, 3
   r(:) = final_grid_points(:,ipoint) 
   r(m) += dx 
   call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
   mu_plus = mu_lda(rho_a_hf,rho_b_hf)
   r(:) = final_grid_points(:,ipoint) 
   r(m) -= dx 
   call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
   mu_minus = mu_lda(rho_a_hf,rho_b_hf)
   grad_mu_of_r_lda(m,ipoint,1) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo

 enddo 
 average_mu_lda    = average_mu_lda   / dble(elec_a+ elec_b)

END_PROVIDER 


double precision function mu_lda(rho_a,rho_b)
 implicit none 
 double precision, intent(in) :: rho_a,rho_b
 include 'constants.include.F'
 double precision :: g0,g0_UEG_mu_inf
 g0 = g0_UEG_mu_inf(rho_a,rho_b)
 mu_lda = - 1.d0 / (dlog(2.d0 * g0) * sqpi) 

end

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
 dx = 1.d-4
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
