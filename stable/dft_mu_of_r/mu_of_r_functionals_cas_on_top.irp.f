

 BEGIN_PROVIDER [double precision, Energy_c_md_n_and_on_top_PBE_mu_of_r, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_n_and_on_top_PBE_mu_of_r = PBE-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction and spin polarization computed only with on-top and total density
 END_DOC
 implicit none
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:)

 allocate(eps_c_md_on_top_PBE(N_states))
 Energy_c_md_n_and_on_top_PBE_mu_of_r = 0.d0
  
 print*,'Providing Energy_c_md_n_and_on_top_PBE_mu_of_r ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i)
  mu = mu_of_r_vector(i)
  call give_epsilon_pbe_ontop_effective_spin_dens_provider(mu,i,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_n_and_on_top_PBE_mu_of_r(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_n_and_on_top_PBE_mu_of_r :',wall1-wall0

 END_PROVIDER


 BEGIN_PROVIDER [double precision, Energy_c_md_n_and_PBE_mu_of_r, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_n_and_on_top_PBE_mu_of_r = PBE-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction and spin polarization computed only with on-top and total density
 END_DOC
 implicit none
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_PBE(:)

 allocate(eps_c_md_PBE(N_states))
 Energy_c_md_n_and_PBE_mu_of_r = 0.d0
  
 print*,'Providing Energy_c_md_n_and_PBE_mu_of_r ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i)
  mu = mu_of_r_vector(i)
  call give_epsilon_pbe_effective_spin_dens_provider(mu,i,eps_c_md_PBE)
  do istate = 1, N_states
   Energy_c_md_n_and_PBE_mu_of_r(istate) += eps_c_md_PBE(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_n_and_PBE_mu_of_r :',wall1-wall0

 END_PROVIDER


 BEGIN_PROVIDER [double precision, Energy_c_md_n_and_on_top_LYP_mu_of_r, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_n_and_on_top_LYP_mu_of_r = LYP-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction and spin polarization computed only with on-top and total density
 END_DOC
 implicit none
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_LYP(:)

 allocate(eps_c_md_on_top_LYP(N_states))
 Energy_c_md_n_and_on_top_LYP_mu_of_r = 0.d0
  
 print*,'Providing Energy_c_md_n_and_on_top_LYP_mu_of_r ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i)
  mu = mu_of_r_vector(i)
  call give_epsilon_lyp_ontop_effective_spin_dens_provider(mu,i,eps_c_md_on_top_LYP)
  do istate = 1, N_states
   Energy_c_md_n_and_on_top_LYP_mu_of_r(istate) += eps_c_md_on_top_LYP(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_n_and_on_top_LYP_mu_of_r :',wall1-wall0

 END_PROVIDER



 BEGIN_PROVIDER [double precision, Energy_c_md_n_and_LYP_mu_of_r, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_n_and_LYP_mu_of_r = LYP- multi determinant functional with on top from the UEG for large mu using a mu(r) interaction and spin polarization computed only with on-top and total density
 END_DOC
 implicit none
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_LYP(:)

 allocate(eps_c_md_LYP(N_states))
 Energy_c_md_n_and_LYP_mu_of_r = 0.d0
  
 print*,'Providing Energy_c_md_n_and_LYP_mu_of_r ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i)
  mu = mu_of_r_vector(i)
  call give_epsilon_lyp_effective_spin_dens_provider(mu,i,eps_c_md_LYP)
  do istate = 1, N_states
   Energy_c_md_n_and_LYP_mu_of_r(istate) += eps_c_md_LYP(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_n_and_LYP_mu_of_r :',wall1-wall0

 END_PROVIDER


