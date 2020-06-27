 subroutine print_ecmd_var_energy
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 provide psi_energy ecmd_pbe_on_top_at_mu ecmd_pbe_ueg_prov 
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  write(*, '(A22,X,F32.10)') 'mu_erf_dft               = ',mu_erf_dft          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY             = ',psi_energy + ecmd_pbe_on_top_at_mu  + nuclear_repulsion 
  write(*, '(A22,X,F16.10)') 'Ecmd PBE-on-top          = ',ecmd_pbe_on_top_at_mu
  print*,  ''
  print*,  'Component of the energy ....'
  print*,  ''
 
  write(*, '(A22,X,F16.10)') '<Psi| H | Psi>           = ',psi_energy + nuclear_repulsion
  print*,'Correlation part '
  print*,'On-top functional '
  write(*, '(A22,X,F16.10)') 'Ecmd PBE-OT              = ',ecmd_pbe_on_top_at_mu
  print*,'On-top-UEG functional '
  write(*, '(A22,X,F16.10)') 'Ecmd PBE-UEG             = ',ecmd_pbe_ueg_prov
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') ' psi_energy_erf          = ',psi_energy_erf 
  print*,''
  write(*, '(A22,X,F16.10)') 'psi_energy_two_e         = ',psi_energy_two_e
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core        = ',psi_energy_h_core
  print*,  '****************************************'
  print*,  'Test for the coherence between density and wave function used'
  print*,  'psi_energy_h_core - psi_dft_energy_h_core  = ',psi_energy_h_core - psi_dft_energy_h_core
 end

 BEGIN_PROVIDER [double precision, int_exmdsrpbe_n2_exact,(N_states)]
&BEGIN_PROVIDER [double precision, int_exmdsrpbe_n2_extrapolated, (N_states)]
&BEGIN_PROVIDER [double precision, int_exmdsrpbe_n2_UEG, (N_states)]

 implicit none
 BEGIN_DOC
 ! Integrate the exmd_sr_pbe using the exact on-top and the extrapolated one
 ! As I don't know yet where to put this providers, I just let it there. Ask Emmanuel where it should go.
 END_DOC
  integer          :: m, istate, ipoint
  double precision :: weight
  double precision :: mu, mu_correction_of_on_top, on_top, on_top_extrap, on_top_UEG
  double precision :: rho, rho_a, rho_b, grad_rho_a(3),grad_rho_b(3),grad_rho_2, grad_rho_a_2, grad_rho_b_2, grad_rho_a_b
  double precision :: dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2,dexdrho2
  double precision :: ex_srmuPBE
  double precision :: g0, dg0drho
  double precision :: r(3), r_norm, f_psi
  do istate = 1, N_states
  int_exmdsrpbe_n2_exact        = 0.d0
  int_exmdsrpbe_n2_extrapolated = 0.d0
  int_exmdsrpbe_n2_UEG          = 0.d0
  do ipoint = 1, n_points_final_grid
   
   mu = mu_erf_dft 
   !call give_mu_of_r_cas(r,istate,mu,f_psi,on_top)
   weight = final_weight_at_r_vector(ipoint)
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   r_norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)

   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   rho = rho_a + rho_b
   grad_rho_a(1:3) = one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) = one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

   !on_top = on_top_cas_mu_r(ipoint,istate) ! we use mu = cas_ful
   on_top = total_cas_on_top_density(ipoint,istate)
!  We take the extrapolated on-top pair density (Eq. 29)
!  Multiplied by 2 because of difference of normalizations between the on_top of QP2 and that of JCP, 150, 084103 1-10 (2019)
   on_top_extrap = 2.d0 * mu_correction_of_on_top(mu,on_top)
   call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
   on_top_UEG = (rho**2)*g0
  !ON-TOP-EXACT 
   call exmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top,ex_srmuPBE,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2,dexdrho2)

   int_exmdsrpbe_n2_exact += ex_srmuPBE*weight

  !ON-TOP EXTRAP
   call exmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top_extrap,ex_srmuPBE,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2,dexdrho2)

   int_exmdsrpbe_n2_extrapolated += ex_srmuPBE*weight

  !ON-TOP UEG
   call exmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top_UEG,ex_srmuPBE,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2,dexdrho2)

   int_exmdsrpbe_n2_UEG += ex_srmuPBE*weight

   !testr:
   !write(33,*) r(1), on_top, on_top_extrap, on_top_UEG, on_top_UEG
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, int_ecmdsrpbe_n2_exact,(N_states)]
&BEGIN_PROVIDER [double precision, int_ecmdsrpbe_n2_extrapolated, (N_states)]
&BEGIN_PROVIDER [double precision, int_ecmdsrpbe_n2_UEG, (N_states)]

 implicit none
 BEGIN_DOC
 ! Integrate the exmd_sr_pbe using the exact on-top and the extrapolated one
 ! As I don't know yet where to put this providers, I just let it there. Ask Emmanuel where it should go.
 END_DOC
  integer          :: m, istate, ipoint
  double precision :: weight
  double precision :: mu, mu_correction_of_on_top, on_top, on_top_extrap, on_top_UEG
  double precision :: rho, rho_a, rho_b, grad_rho_a(3),grad_rho_b(3),grad_rho_2, grad_rho_a_2, grad_rho_b_2, grad_rho_a_b
  double precision :: ec_srmuPBE
  double precision :: g0, dg0drho
  double precision :: r(3), r_norm, f_psi
  double precision :: decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2

  do istate = 1, N_states
  int_ecmdsrpbe_n2_exact(istate)        = 0.d0
  int_ecmdsrpbe_n2_extrapolated(istate) = 0.d0
  int_ecmdsrpbe_n2_UEG(istate)          = 0.d0
  do ipoint = 1, n_points_final_grid
   
   mu = mu_erf_dft 
   !call give_mu_of_r_cas(r,istate,mu,f_psi,on_top)
   weight = final_weight_at_r_vector(ipoint)
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   r_norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)

   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   rho = rho_a + rho_b
   grad_rho_a(1:3) = one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) = one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

   on_top = total_cas_on_top_density(ipoint,istate)
   on_top_extrap = 2.d0 * mu_correction_of_on_top(mu,on_top)
   call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
   on_top_UEG = (rho**2)*g0

  call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top_extrap,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
  int_ecmdsrpbe_n2_extrapolated(istate) += ec_srmuPBE*weight

  call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
  int_ecmdsrpbe_n2_exact(istate) += ec_srmuPBE*weight

  call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top_UEG,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
  int_ecmdsrpbe_n2_UEG(istate) += ec_srmuPBE*weight

  enddo
 enddo
 END_PROVIDER

 subroutine print_many_energy_component
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 double precision :: energy(N_states)
 provide psi_energy ecmd_pbe_on_top_at_mu ecmd_pbe_ueg_prov ecmd_pbe_ueg_prov energy_x_md_sr_pbe
 print*,'/////////////////////////'
  print*,  '****************************************'
  write(*, '(A35,X,F32.10)') 'mu_erf_dft               = ',mu_erf_dft          
  print*,'///////////////////'
  !       <Psi | T + V_ne        +  W_ee^lr      | Psi >  + E_x^sr  +  E_c^sr         + E_H^sr 
  energy = psi_dft_energy_h_core + psi_energy_erf + energy_x_sr_pbe + energy_c_sr_pbe +  short_range_Hartree
  write(*, '(A35,X,F16.10)') 'Usual sr_pbe            = ',energy + nuclear_repulsion
  write(*, '(A35,X,F16.10)') 'kinetic energy DFT      = ',psi_dft_energy_kinetic
  write(*, '(A35,X,F16.10)') 'Nuclear repulsion DFT   = ',psi_dft_energy_nuclear_elec 
  energy = psi_energy + ecmd_pbe_on_top_at_mu  
  write(*, '(A35,X,F16.10)') 'PBE-OT energy            = ',energy + nuclear_repulsion 
  write(*, '(A35,X,F16.10)') 'Test                     = ',ecmd_pbe_on_top_at_mu + nuclear_repulsion 
  !       <Psi | T + V_ne        +  W_ee^lr      | Psi >  + E_xmd^sr   +  E_cmd^sr         + E_H^sr 
  energy = psi_dft_energy_h_core + psi_energy_erf + energy_x_md_sr_pbe + energy_c_md_sr_pbe +  short_range_Hartree
  write(*, '(A35,X,F16.10)') 'PBE-UEG MD               = ',energy + nuclear_repulsion
  print*,  ''
  print*,  'Component of the energy ....'
  print*,  ''
 
  write(*, '(A35,X,F16.10)') '<Psi| H | Psi>           = ',psi_energy + nuclear_repulsion
  write(*, '(A35,X,F16.10)') 'Short range Hartree      = ',short_range_Hartree
  print*,'Correlation part '
  print*,'On-top functional '
  write(*, '(A35,X,F16.10)') 'ecmd_pbe_on_top_at_mu     = ',ecmd_pbe_on_top_at_mu
  write(*, '(A35,X,F16.10)') 'exmdsr_pbe_on_top_exact  = ',int_exmdsrpbe_n2_exact
  write(*, '(A35,X,F16.10)') 'ecmdsr_pbe_on_top_exact  = ',int_ecmdsrpbe_n2_exact
  write(*, '(A35,X,F16.10)') 'exmdsr_pbe_on_top_extra  = ',int_exmdsrpbe_n2_extrapolated
  write(*, '(A35,X,F16.10)') 'ecmdsr_pbe_on_top_extra   = ',int_ecmdsrpbe_n2_extrapolated
  write(*, '(A35,X,F16.10)') 'exmdsr_pbe_on_top_UEG    = ',int_exmdsrpbe_n2_UEG
  write(*, '(A35,X,F16.10)') 'ecmdsr_pbe_on_top_UEG    = ',int_ecmdsrpbe_n2_UEG
  print*,'On-top-UEG functional '
  write(*, '(A35,X,F16.10)') 'Ecmd PBE-UEG              = ',ecmd_pbe_ueg_prov
  write(*, '(A35,X,F16.10)') 'energy_c_md_sr_pbe        = ',energy_c_md_sr_pbe
  print*,  '****************************************'
  print*,'Exchange    part '
  write(*, '(A35,X,F16.10)') ' Exact EXmd energy       = ',psi_energy_wee_sr - short_range_Hartree
!  write(*, '(A35,X,F16.10)') ' Exact ECmd energy       = ',psi_energy_wee_sr - short_range_Hartree
  write(*, '(A35,X,F16.10)') 'energy_x_md_sr_pbe       = ',energy_x_md_sr_pbe
  print*,  '****************************************'
  write(*, '(A35,X,F16.10)') ' psi_energy_erf          = ',psi_energy_erf 
  print*,''
  write(*, '(A35,X,F16.10)') 'psi_energy_two_e         = ',psi_energy_two_e
  write(*, '(A35,X,F16.10)') 'psi_energy_h_core        = ',psi_energy_h_core
  print*,  '****************************************'
  print*,  'Test for the coherence between density and wave function used'
  print*,  'psi_energy_h_core - psi_dft_energy_h_core  = ',psi_energy_h_core - psi_dft_energy_h_core
 end

