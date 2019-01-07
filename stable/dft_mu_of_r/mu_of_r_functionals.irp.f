
 BEGIN_PROVIDER [double precision, Energy_c_md_mu_of_r_LDA, (N_states)]
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

 Energy_c_md_mu_of_r_LDA = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = mu_of_r_vector(i)
  weight=final_weight_functions_at_final_grid_points(i)
  do istate = 1, N_states
   rho_a(istate) = one_body_dm_alpha_at_r(i,istate)
   rho_b(istate) = one_body_dm_beta_at_r(i,istate)
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   Energy_c_md_mu_of_r_LDA(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for Energy_c_md_mu_of_r_LDA :',wall1-wall0


 END_PROVIDER 

 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_of_r_UEG, (N_states)]
&BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_of_r, (N_states)]
&BEGIN_PROVIDER [double precision, Energy_c_md_on_top_mu_of_r, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_on_top_PBE_mu_of_r_UEG_vector = PBE-on_top multi determinant functional with exact on top extracted from the UEG using a mu(r) interaction 
  ! Energy_c_md_on_top_PBE_mu_of_r_vector     = PBE-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction 
  ! Energy_c_md_on_top_u_of_r_vector     = on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction 
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu,pi
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 pi = 4.d0 * datan(1.d0)

 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 Energy_c_md_on_top_PBE_mu_of_r_UEG = 0.d0
 Energy_c_md_on_top_PBE_mu_of_r = 0.d0
 Energy_c_md_on_top_mu_of_r = 0.d0
  
 print*,'Providing Energy_c_md_mu_of_r_PBE_on_top ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight=final_weight_functions_at_final_grid_points(i)
  two_dm(:) = on_top_of_r_vector(i,:)
  mu = mu_of_r_vector(i)

  call give_epsilon_c_md_on_top_PBE_mu_corrected_UEG_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_of_r_UEG(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
  call give_epsilon_c_md_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_of_r(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
  do istate = 1, N_states
   if(mu.gt.1.d-10)then
    Energy_c_md_on_top_mu_of_r(istate) += ((-2.d0+sqrt(2.d0))*sqrt(2.d0*pi)/(3.d0*(mu**3))) * two_dm(istate) * weight 
   endif
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_on_top_PBE_mu_of_r:',wall1-wall0

 END_PROVIDER

 BEGIN_PROVIDER [double precision, mu_average, (N_states)]
 implicit none 
 BEGIN_DOC
 ! mu_average = \int dr mu(r) * n(r) /N_elec
 END_DOC
 integer :: i,istate
 double precision, allocatable :: rho_a(:), rho_b(:)
 double precision :: mu,weight
 mu_average = 0.d0
 allocate(rho_a(N_states), rho_b(N_states))
 do i = 1, n_points_final_grid
  mu = mu_of_r_vector(i)
  weight=final_weight_functions_at_final_grid_points(i)
  do istate = 1, N_states
   rho_a(istate) = one_body_dm_alpha_at_r(i,istate)
   rho_b(istate) = one_body_dm_beta_at_r(i,istate)
   mu_average(istate) +=  weight * mu_of_r_vector(i) * (rho_a(istate) + rho_b(istate))
  enddo
 enddo
 mu_average = mu_average / dble(elec_alpha_num + elec_beta_num)
 END_PROVIDER 

BEGIN_PROVIDER [double precision, Energy_c_LDA_mu_of_r, (N_states)]
 implicit none
 BEGIN_DOC 
! Energy_c_LDA_mu_of_r = similar of equation 40 of J. Chem. Phys. 149, 194301 (2018) but using a pure correlation part for the LDA (so not the ECMD)
 END_DOC
 integer :: i,istate
 double precision :: mu,weight,e_lda
 double precision, allocatable :: rho_a(:), rho_b(:)
 allocate(rho_a(N_states), rho_b(N_states))
 
 energy_c_LDA_mu_of_r = 0.d0
 
 do i = 1, n_points_final_grid
  mu = mu_of_r_vector(i)
  weight=final_weight_functions_at_final_grid_points(i)
  do istate = 1, N_states
   rho_a(istate) = one_body_dm_alpha_at_r(i,istate)
   rho_b(istate) = one_body_dm_beta_at_r(i,istate)
   call ec_only_lda_sr(mu,rho_a(istate),rho_b(istate),e_lda)
   energy_c_LDA_mu_of_r(istate) += weight * e_lda
  enddo
 enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, HF_alpha_beta_mu_of_r_bielec_energy]
&BEGIN_PROVIDER [double precision, HF_alpha_beta_bielec_energy]
 implicit none
 BEGIN_DOC
! HF_alpha_beta_mu_of_r_bielec_energy = <HF|W_ee^{mu(r)}_{alpha/beta}|HF>
! HF_alpha_beta_bielec_energy         = <HF|    W_ee_{alpha/beta}    |HF>
 END_DOC
 integer :: i,j,k
 double precision :: r(3), weight,tmp,integral_of_mu_of_r_on_HF
 double precision, allocatable :: integrals_mo(:,:),mos_array(:)
 allocate(integrals_mo(mo_tot_num,mo_tot_num),mos_array(mo_tot_num))
 HF_alpha_beta_mu_of_r_bielec_energy = 0.d0
  
 integer                        :: occ(N_int*bit_kind_size,2)
 call bitstring_to_list(ref_bitmask(1,1), occ(1,1), i, N_int)
 call bitstring_to_list(ref_bitmask(1,2), occ(1,2), i, N_int)
 
 do j= 1, elec_beta_num
  do i= 1, elec_alpha_num
    HF_alpha_beta_mu_of_r_bielec_energy += mo_bielec_integral_erf_mu_of_r_jj(occ(i,1),occ(j,2))
  enddo
 enddo

 HF_alpha_beta_bielec_energy = 0.d0
 do j= 1, elec_beta_num
  do i= 1, elec_alpha_num
    HF_alpha_beta_bielec_energy += mo_bielec_integral_jj(occ(i,1),occ(j,2))
  enddo
 enddo


END_PROVIDER 

