 BEGIN_PROVIDER [double precision, integrals_for_hf_potential_ecmd, (mo_num,mo_num,elec_alpha_num,elec_alpha_num)]
 implicit none
 integer :: i,j,m,n
 double precision :: get_mo_two_e_integral_erf

 do m = 1, elec_alpha_num ! electron 1 
  do n = 1, elec_alpha_num ! electron 2 
   do i = 1, mo_num   ! electron 1 
    do j = 1, mo_num  ! electron 2 
     integrals_for_hf_potential_ecmd(j,i,n,m) = get_mo_two_e_integral_erf(m,n,i,j,mo_integrals_erf_map) 
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER 

subroutine f_HF_ab_ecmd(r1,r2,integral_psi,two_bod)
 implicit none
 BEGIN_DOC
! f_HF_ab(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha beta spins and an HF wave function with a long range operator
! < HF | wee_{\alpha\beta}^{lr} | HF > = 0.5 * \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral_psi,two_bod
 integer :: i,j,m,n
 double precision :: mos_array_r1(mo_num)
 double precision :: mos_array_r2(mo_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 integral_psi = 0.d0
 two_bod = 0.d0
 do m = 1, elec_alpha_num
  do n = 1, elec_beta_num
   two_bod += mos_array_r1(n) * mos_array_r1(n) * mos_array_r2(m) * mos_array_r2(m) 
   do i = 1, mo_num
    do j = 1, mo_num
     integral_psi +=  integrals_for_hf_potential_ecmd(j,i,n,m) * mos_array_r1(i) * mos_array_r2(j) * mos_array_r2(n) * mos_array_r1(m) 
    enddo
   enddo
  enddo
 enddo
end


 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_vector_ecmd, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_ecmd, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_bod
 print*,'providing the mu_of_r_hf_coal_vector_ecmd ...'
 call wall_time(cpu0)
 r = 0.d0
 call f_HF_ab_ecmd(r,r,local_potential,two_bod)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_bod) & 
 !$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_vector_ecmd) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call f_HF_ab_ecmd(r,r,local_potential,two_bod)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0)then
    local_potential = 1.d-10
  else 
    local_potential = local_potential /  two_bod
  endif
  mu_of_r_hf_coal_vector_ecmd(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coal_ecmd(k,i,j) = mu_of_r_hf_coal_vector_ecmd(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coal_vector_ecmd = ',cpu1-cpu0
 END_PROVIDER



 BEGIN_PROVIDER [double precision, Energy_c_md_mu_of_r_LDA_rsdft, (N_states)]
 BEGIN_DOC
 !Energy_c_md_mu_of_r_LDA_rsdft = multi-determinantal correlation functional with mu(r) , see equation 40 in J. Chem. Phys. 149, 194301 (2018); https://doi.org/10.1063/1.5052714
 END_DOC
 implicit none 
 integer :: i,istate 
 double precision, allocatable :: rho_a(:), rho_b(:), ec(:)
 logical :: dospin
 double precision :: wall0,wall1,weight,mu
 dospin = .true. ! JT dospin have to be set to true for open shell
 print*,'Providing Energy_c_md_mu_of_r_LDA_rsdft ...'
 allocate(rho_a(N_states), rho_b(N_states), ec(N_states))

 Energy_c_md_mu_of_r_LDA_rsdft = 0.d0
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  mu = min(mu_erf,mu_of_r_hf_coal_vector_ecmd(i))
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_alpha_at_r(i,istate)
   rho_b(istate) = one_e_dm_beta_at_r(i,istate)
   call ESRC_MD_LDAERF (mu,rho_a(istate),rho_b(istate),dospin,ec(istate))
   Energy_c_md_mu_of_r_LDA_rsdft(istate) += weight * ec(istate)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time for Energy_c_md_mu_of_r_LDA_rsdft :',wall1-wall0

 END_PROVIDER 
 
 
 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_of_r_UEG_rsdft, (N_states)]
&BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_of_r_rsdft, (N_states)]
&BEGIN_PROVIDER [double precision, Energy_c_md_on_top_mu_of_r_rsdft, (N_states)]
 BEGIN_DOC
  ! Energy_c_md_on_top_PBE_mu_of_r_UEG_rsdft_vector = PBE-on_top multi determinant functional with exact on top extracted from the UEG using a mu(r) interaction 
  ! Energy_c_md_on_top_PBE_mu_of_r_rsdft_vector     = PBE-on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction 
  ! Energy_c_md_on_top_u_of_r_vector     = on_top multi determinant functional with exact on top extrapolated for large mu using a mu(r) interaction 
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu,pi
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 pi = 4.d0 * datan(1.d0)

 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 Energy_c_md_on_top_PBE_mu_of_r_UEG_rsdft = 0.d0
 Energy_c_md_on_top_PBE_mu_of_r_rsdft = 0.d0
 Energy_c_md_on_top_mu_of_r_rsdft = 0.d0
  
 print*,'Providing Energy_c_md_mu_of_r_PBE_on_top ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight=final_weight_at_r_vector(i)
  two_dm(:) = on_top_of_r_vector(i,:)
  mu = min(mu_erf,mu_of_r_hf_coal_vector_ecmd(i))

  call give_epsilon_c_md_on_top_PBE_mu_corrected_UEG_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_of_r_UEG_rsdft(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
  call give_epsilon_c_md_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_of_r_rsdft(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
  do istate = 1, N_states
   if(mu.gt.1.d-10)then
    Energy_c_md_on_top_mu_of_r_rsdft(istate) += ((-2.d0+sqrt(2.d0))*sqrt(2.d0*pi)/(3.d0*(mu**3))) * two_dm(istate) * weight 
   endif
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the Energy_c_md_on_top_PBE_mu_of_r_rsdft:',wall1-wall0

 END_PROVIDER
