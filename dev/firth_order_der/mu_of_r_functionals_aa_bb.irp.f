
 BEGIN_PROVIDER [double precision, e_c_md_mur_aa_LDA_a, (N_states)]
 BEGIN_DOC
! Energy_c_md_mu_of_r_LDA = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_aa_LDA_a ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_c_md_mur_aa_LDA_a = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector_alpha(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_alpha_at_r(i,istate)
   rho_b(istate) = 0.d0 
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   e_c_md_mur_aa_LDA_a(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for Energy_c_md_mu_of_r_aa_LDA_a :',wall1-wall0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, e_c_md_mur_bb_LDA_b, (N_states)]
 BEGIN_DOC
! Energy_c_md_mu_of_r_LDA = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_bb_LDA_b ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_c_md_mur_bb_LDA_b = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector_beta(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = 0.d0 
   rho_b(istate) = one_e_dm_beta_at_r(i,istate) 
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   e_c_md_mur_bb_LDA_b(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for Energy_c_md_mu_of_r_bb_LDA_b :',wall1-wall0
 END_PROVIDER 

!!!!!!!!!!!!! mu_r alpha beta!!!!!!!!!!!!!!!!!!!!!!!!!!


 BEGIN_PROVIDER [double precision, e_c_md_mur_ab_LDA_a, (N_states)]
 BEGIN_DOC
! Energy_c_md_mu_of_r_LDA = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_aa_LDA_a ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_c_md_mur_ab_LDA_a = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_alpha_at_r(i,istate)
   rho_b(istate) = 0.d0 
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   e_c_md_mur_ab_LDA_a(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for Energy_c_md_mu_of_r_ab_LDA_a :',wall1-wall0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, e_c_md_mur_ab_LDA_b, (N_states)]
 BEGIN_DOC
! Energy_c_md_mu_of_r_LDA = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_ab_LDA_b ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_c_md_mur_ab_LDA_b = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = 0.d0 
   rho_b(istate) = one_e_dm_beta_at_r(i,istate) 
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   e_c_md_mur_ab_LDA_b(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for Energy_c_md_mu_of_r_ab_LDA_b :',wall1-wall0
 END_PROVIDER 


