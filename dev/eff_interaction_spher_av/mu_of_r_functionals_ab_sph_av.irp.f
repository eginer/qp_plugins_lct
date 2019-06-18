
 BEGIN_PROVIDER [double precision, e_c_md_mur_ab_sph_av_LDA, (N_states)]
 BEGIN_DOC
! Energy_c_md_mu_of_r_LDA = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_LDA ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_c_md_mur_ab_sph_av_LDA = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector_ab_sph_av(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_alpha_at_r(i,istate)
   rho_b(istate) = one_e_dm_beta_at_r(i,istate) 
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   e_c_md_mur_ab_sph_av_LDA(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for e_c_md_mur_ab_sph_av_LDA :',wall1-wall0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, e_c_md_mur_ab_sph_av_LDA_a, (N_states)]
 BEGIN_DOC
! Energy_c_md_mu_of_r_LDA = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_LDA ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_c_md_mur_ab_sph_av_LDA_a = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector_ab_sph_av(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_alpha_at_r(i,istate)
   rho_b(istate) = 0.d0
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   e_c_md_mur_ab_sph_av_LDA_a(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for e_c_md_mur_ab_sph_av_LDA_a :',wall1-wall0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, e_c_md_mur_ab_sph_av_LDA_b, (N_states)]
 BEGIN_DOC
! Energy_c_md_mu_of_r_LDA = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_LDA ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 e_c_md_mur_ab_sph_av_LDA_b = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector_ab_sph_av(i)
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_beta_at_r(i,istate)
   rho_b(istate) = 0.d0
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   e_c_md_mur_ab_sph_av_LDA_b(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for e_c_md_mur_ab_sph_av_LDA_b :',wall1-wall0
 END_PROVIDER 

