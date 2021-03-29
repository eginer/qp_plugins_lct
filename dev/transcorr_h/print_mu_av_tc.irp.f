program print_mu_av_tc
 implicit none
 read_wf = .True.
 touch read_wf
   print*,'average_mu_lda       = ',average_mu_lda
   print*,'average_mu_rs        = ',average_mu_rs 
   print*,'average_mu_rs_c      = ',average_mu_rs_c
   print*,'average_mu_rs_c_lda  = ',average_mu_rs_c_lda
   print*,'average_mu_grad_n    = ',average_mu_grad_n
   print*,'average_mu_lda * 1/2 = ',average_mu_lda * 0.5d0
end


 BEGIN_PROVIDER [double precision, average_mu_lda      ]
&BEGIN_PROVIDER [double precision, average_mu_rs       ]
&BEGIN_PROVIDER [double precision, average_mu_rs_c     ]
&BEGIN_PROVIDER [double precision, average_mu_rs_c_lda ]
&BEGIN_PROVIDER [double precision, average_mu_grad_n   ]
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
 average_mu_lda     = 0.d0
 average_mu_rs      = 0.d0
 average_mu_rs_c    = 0.d0
 average_mu_grad_n  = 0.d0 
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
  grad_n = grad_n/(4.d0 * rho_hf)
  rs = cst_rs * rho_hf**(-1.d0/3.d0)
  g0 = g0_UEG_mu_inf(rho_a_hf,rho_b_hf)

  average_mu_rs     += 1.d0/rs * rho_hf * weight
  average_mu_rs_c   += alpha_rs/dsqrt(rs) * rho_hf * weight
  average_mu_lda    +=  - 1.d0 / (dlog(2.d0 * g0) * sqpi) * weight * rho_hf
  average_mu_grad_n += grad_n * rho_hf * weight
 enddo 
 average_mu_lda    = average_mu_lda   / dble(elec_alpha_num + elec_beta_num)
 average_mu_rs     = average_mu_rs    / dble(elec_alpha_num + elec_beta_num)
 average_mu_rs_c   = average_mu_rs_c  / dble(elec_alpha_num + elec_beta_num)
 average_mu_grad_n = average_mu_grad_n/ dble(elec_alpha_num + elec_beta_num)

 average_mu_rs_c_lda = 0.5d0 * (average_mu_lda + average_mu_rs_c)

END_PROVIDER 
