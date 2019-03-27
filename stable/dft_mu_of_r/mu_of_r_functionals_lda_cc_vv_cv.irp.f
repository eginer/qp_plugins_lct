
 BEGIN_PROVIDER [double precision, e_cmd_mu_of_r_lda_vv, (N_states)]
 BEGIN_DOC
! e_cmd_mu_of_r_lda_vv = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing e_cmd_mu_of_r_lda_vv ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_cmd_mu_of_r_lda_vv = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vv_vector(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_no_core_and_grad_alpha_in_r(4,i,istate)
   rho_b(istate) = one_e_dm_no_core_and_grad_beta_in_r(4,i,istate)
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   if(isnan(ec(istate)))then
    print*,'ec is nan'
    stop
   endif
   e_cmd_mu_of_r_lda_vv(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for e_cmd_mu_of_r_lda_vv :',wall1-wall0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, e_cmd_mu_of_r_lda_cc, (N_states)]
 BEGIN_DOC
! e_cmd_mu_of_r_lda_vv = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing e_cmd_mu_of_r_lda_cc ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_cmd_mu_of_r_lda_cc = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_cc_vector(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) =  one_e_dm_and_grad_alpha_in_r(4,i,istate) - one_e_dm_no_core_and_grad_alpha_in_r(4,i,istate)
   rho_b(istate) =  one_e_dm_and_grad_beta_in_r(4,i,istate)  - one_e_dm_no_core_and_grad_beta_in_r(4,i,istate)
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   if(isnan(ec(istate)))then
    print*,'ec is nan'
    stop
   endif
   e_cmd_mu_of_r_lda_cc(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for e_cmd_mu_of_r_lda_cc :',wall1-wall0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, e_cmd_mu_of_r_lda_cv, (N_states)]
 BEGIN_DOC
! e_cmd_mu_of_r_lda_vv = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing e_cmd_mu_of_r_lda_cv ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_cmd_mu_of_r_lda_cv = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_cv_vector(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) =  (one_e_dm_and_grad_alpha_in_r(4,i,istate) - one_e_dm_no_core_and_grad_alpha_in_r(4,i,istate)) & 
                   * one_e_dm_no_core_and_grad_alpha_in_r(4,i,istate) 
   rho_a(istate) = dsqrt(rho_a(istate))
   rho_b(istate) =  (one_e_dm_and_grad_beta_in_r(4,i,istate) - one_e_dm_no_core_and_grad_beta_in_r(4,i,istate)) & 
                   * one_e_dm_no_core_and_grad_beta_in_r(4,i,istate) 
   rho_b(istate) = dsqrt(rho_b(istate))
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   if(isnan(ec(istate)))then
    print*,'ec is nan'
    stop
   endif
   e_cmd_mu_of_r_lda_cv(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for e_cmd_mu_of_r_lda_cv :',wall1-wall0
 END_PROVIDER 
