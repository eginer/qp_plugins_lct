

 BEGIN_PROVIDER [double precision, Energy_c_md_n_and_on_top_SCAN_mu_of_r, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_n_and_on_top_SCAN_mu_of_r = SCAN-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction and spin polarization computed only with on-top and total density
 END_DOC
 implicit none
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_SCAN(:)

 allocate(eps_c_md_on_top_SCAN(N_states))
 Energy_c_md_n_and_on_top_SCAN_mu_of_r = 0.d0
  
 print*,'Providing Energy_c_md_n_and_on_top_SCAN_mu_of_r ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i)
  mu = mu_of_r_vector(i)
  call give_epsilon_scan_ontop_effective_spin_dens_provider(mu,i,eps_c_md_on_top_SCAN)
  do istate = 1, N_states
   Energy_c_md_n_and_on_top_SCAN_mu_of_r(istate) += eps_c_md_on_top_SCAN(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_n_and_on_top_SCAN_mu_of_r :',wall1-wall0

 END_PROVIDER




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



 BEGIN_PROVIDER [double precision, Energy_c_md_n_and_SCAN_mu_of_r, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_n_and_on_top_SCAN_mu_of_r = SCAN-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction and spin polarization computed only with on-top and total density
 END_DOC
 implicit none
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_SCAN(:)

 allocate(eps_c_md_SCAN(N_states))
 Energy_c_md_n_and_SCAN_mu_of_r = 0.d0
  
 print*,'Providing Energy_c_md_n_and_SCAN_mu_of_r ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i)
  mu = mu_of_r_vector(i)
  call give_epsilon_scan_effective_spin_dens_provider(mu,i,eps_c_md_SCAN)
  do istate = 1, N_states
   Energy_c_md_n_and_SCAN_mu_of_r(istate) += eps_c_md_SCAN(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_n_and_SCAN_mu_of_r :',wall1-wall0

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


 BEGIN_PROVIDER [double precision, Energy_c_md_n_and_LDA_mu_of_r, (N_states)]
 BEGIN_DOC
! Energy_c_md_LDA_mu_of_r = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_LDA_mu_of_r ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 Energy_c_md_n_and_LDA_mu_of_r = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_vector(i)
  weight= final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = effective_alpha_dm(i,istate)
   rho_b(istate) = effective_beta_dm(i,istate)
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   if(isnan(ec(istate)))then
    print*,'ec is nan'
    stop
   endif
   Energy_c_md_n_and_LDA_mu_of_r(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for Energy_c_md_n_and_LDA_mu_of_r :',wall1-wall0
 END_PROVIDER 
