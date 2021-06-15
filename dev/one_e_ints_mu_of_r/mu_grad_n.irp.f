
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
