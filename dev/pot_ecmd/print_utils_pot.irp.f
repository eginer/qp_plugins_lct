
 subroutine print_ecmd_var_energy_barth
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
 !print*, 'mu_erf_dft                 = ',mu_erf_dft          
 !print*, 'TOTAL ENERGY PBE           = ',psi_energy + Energy_c_md_on_top_PBE_mu_vector  + nuclear_repulsion 
 !print*, 'TOTAL ENERGY mu_of_r PBE   = ',psi_energy + Energy_c_md_on_top_PBE_mu_of_r_rsdft  + nuclear_repulsion 
  print*, '********* FOR LDA Calculation*****'
  print*, 'TOTAL ENERGY LDA                  = ',psi_energy + Energy_c_md_LDA + nuclear_repulsion
  print*, 'TOTAL ENERGY mu_of_r LDA  rsdft   = ',psi_energy + Energy_c_md_mu_of_r_LDA_rsdft  + nuclear_repulsion 
  print*, 'TOTAL ENERGY mu_of_r LDA          = ',psi_energy + Energy_c_md_LDA_mu_of_r  + nuclear_repulsion 
  print*, 'TOTAL ENERGY mu_of_r PBE UEG      = ',psi_energy + Energy_c_md_on_top_PBE_mu_of_r_UEG + nuclear_repulsion 
  print*, ''
  print*, ''
  print*, 'Component of the energy ....'
  print*, 'psi_energy                        = ',psi_energy
  print*, 'nuclear_repulsion                 = ',nuclear_repulsion
  print*, 'Energy_c_md_on_top_PBE_mu_of_r_UEG = ',Energy_c_md_on_top_PBE_mu_of_r_UEG 
  print*, 'mu_average                         = ',mu_average                         
  print*, ''
  print*, 'YOLOOOOLLLOOO'
 
 !print*, 'Energy_c_md_on_top_PBE_mu_vector  = ',Energy_c_md_on_top_PBE_mu_vector
 !print*, 'Energy_c_md_on_top_PBE_mu_of_r    = ',Energy_c_md_on_top_PBE_mu_of_r_rsdft
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') 'psi_energy_two_e    = ',psi_energy_two_e
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core   = ',psi_energy_h_core
  print*,  '****************************************'
  print*,'Test for the coherence between density and wave function used'
  print*,'psi_energy_h_core - psi_dft_energy_h_core = ',psi_energy_h_core - psi_dft_energy_h_core
  print*, ''
  print*, ''
  print*,  '****************************************'
  print*, 'z_dipole_moment                  = ',z_dipole_moment
  print*,  '****************************************'

 end



 BEGIN_PROVIDER [double precision, z_dipole_moment, (N_states)]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,i,j 
 double precision :: weight, r(3)  
 double precision :: cpu0,cpu1
 call cpu_time(cpu0)
 z_dipole_moment = 0.d0

 do i = 1, mo_num  
  do j = 1, mo_num  
   do istate = 1, N_states
    z_dipole_moment(istate) += - (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_z(j,i) 
   enddo
  enddo
 enddo
 
 do istate = 1, N_states
  do i = 1,nucl_num 
   z_dipole_moment(istate) += nucl_charge(i) * nucl_coord(i,3) 
  enddo
 enddo

 call cpu_time(cpu1)
 print*,'Time for the z_dipole_moment integration :',cpu1-cpu0
END_PROVIDER




 subroutine print_z_dipole_moment_only
 implicit none
  print*, ''
  print*, ''
  print*,  '****************************************'
  print*, 'z_dipole_moment                  = ',z_dipole_moment
  print*,  '****************************************'
 end
