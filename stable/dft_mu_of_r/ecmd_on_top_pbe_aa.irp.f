BEGIN_PROVIDER [double precision, Ecmd_mu_of_r_pbe_n2_hf_aa, (N_states)]
 implicit none
 BEGIN_DOC
 ! Pure alpha-alpha Ecmd with the second derivative of the HF on-top pair density and the alpha/beta mu(r)
 END_DOC
 integer :: i_state,i_point
 double precision :: weight
 double precision :: constant,beta_tmp,extrapol_tmp
 constant = (-3.d0 + 2.d0 * dsqrt(2.d0)) * dsqrt(dacos(-1.d0)) / (10.d0*dsqrt(2.d0) ) 
 extrapol_tmp = 2.d0 / (3.d0 * dsqrt(dacos(-1.d0))) ! extrapolation factor for the two body density 
          
 Ecmd_mu_of_r_pbe_n2_hf_aa = 0.d0
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i_point)
   if(dabs(e_c_pbe_on_grid(1,i_point,i_state)).lt.1.d-13.or.dabs(nabla2_n2_hf_alpha_alpha(i_point)).lt.1.d-13)cycle
   beta_tmp = e_c_pbe_on_grid(1,i_point,i_state)/(constant  *  extrapol_tmp * mu_of_r_vector(i_point) * nabla2_n2_hf_alpha_alpha(i_point)) 
   Ecmd_mu_of_r_pbe_n2_hf_aa += e_c_pbe_on_grid(1,i_point,i_state) / (1.d0 + beta_tmp * mu_of_r_vector(i_point) **5) * weight 
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, Ecmd_mu_of_r_pbe_n2_hf_bb, (N_states)]
 implicit none
 BEGIN_DOC
 ! Pure alpha-alpha Ecmd with the second derivative of the HF on-top pair density and the alpha/beta mu(r)
 END_DOC
 integer :: i_state,i_point
 double precision :: weight
 double precision :: constant,beta_tmp,extrapol_tmp
 constant = (-3.d0 + 2.d0 * dsqrt(2.d0)) * dsqrt(dacos(-1.d0)) / (10.d0*dsqrt(2.d0) ) 
 extrapol_tmp = 2.d0 / (3.d0 * dsqrt(dacos(-1.d0))) ! extrapolation factor for the two body density 
          
 Ecmd_mu_of_r_pbe_n2_hf_bb = 0.d0
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i_point)
   if(dabs(e_c_pbe_on_grid(2,i_point,i_state)).lt.1.d-13.or.dabs(nabla2_n2_hf_beta_beta(i_point)).lt.1.d-13)cycle
   beta_tmp = e_c_pbe_on_grid(2,i_point,i_state)/(constant  *  extrapol_tmp * mu_of_r_vector(i_point) * nabla2_n2_hf_beta_beta(i_point)) 
   Ecmd_mu_of_r_pbe_n2_hf_bb += e_c_pbe_on_grid(2,i_point,i_state) / (1.d0 + beta_tmp * mu_of_r_vector(i_point) **5) * weight 
  enddo
 enddo

END_PROVIDER 

