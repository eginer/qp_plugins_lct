
 subroutine print_ecmd_var_energy_barth
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*, 'TOTAL ENERGY mu_of_r LDA          = ',psi_energy + ecmd_lda_mu_of_r        + nuclear_repulsion 
  print*, 'TOTAL ENERGY mu_of_r PBE UEG      = ',psi_energy + ecmd_pbe_ueg_mu_of_r    + nuclear_repulsion 
  print*, ''
  print*, ''
  print*, 'Component of the energy ....'
  print*, 'psi_energy                        = ',psi_energy
  print*, 'nuclear_repulsion                 = ',nuclear_repulsion
  print*, '   ecmd_pbe_ueg_mu_of_r           = ',ecmd_pbe_ueg_mu_of_r   
  print*, '   ecmd_lda_mu_of_r               = ',ecmd_lda_mu_of_r
  print*, 'mu_average                        = ',mu_average_prov                    
  print*, ''
  print*, 'YOLOOOOLLLOOO'
 
 !print*, 'Energy_c_md_on_top_PBE_mu_vector  = ',Energy_c_md_on_top_PBE_mu_vector
 !print*, 'Energy_c_md_on_top_PBE_mu_of_r    = ',Energy_c_md_on_top_PBE_mu_of_r_rsdft
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') 'psi_energy_two_e    = ',psi_energy_two_e
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core   = ',psi_energy_h_core
  print*,  '****************************************'
  print*, ''
  print*,  '****************************************'
  print*, 'x_dipole_moment                  = ',x_dipole_moment
  print*, 'y_dipole_moment                  = ',y_dipole_moment
  print*, 'z_dipole_moment                  = ',z_dipole_moment
  print*,  '****************************************'
  print*, 'total_dipole_moment                  = ',total_dipole_moment
 end

 BEGIN_PROVIDER [double precision, x_dipole_moment, (N_states)]
&BEGIN_PROVIDER [double precision, x_dipole_moment_diag, (N_states)]
&BEGIN_PROVIDER [double precision, x_dipole_moment_off_diag, (N_states)]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,i,j 
 double precision :: weight, r(3)  
 double precision :: cpu0,cpu1,nuclei_part

 call cpu_time(cpu0)
 x_dipole_moment = 0.d0
 x_dipole_moment_diag = 0.d0
 x_dipole_moment_off_diag= 0.d0
 do istate = 1, N_states
  do i = 1, mo_num  
   x_dipole_moment_diag(istate) += -(one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)) *  mo_dipole_x(i,i)
   do j = i+1, mo_num  
    x_dipole_moment_off_diag(istate) += - 2.d0 * (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_x(j,i) 
   enddo
  enddo
  x_dipole_moment(istate) = x_dipole_moment_diag(istate) + x_dipole_moment_off_diag(istate)
 enddo
 
 print*,'electron part for x_dipole = ',x_dipole_moment
 
 nuclei_part = 0.d0
 do i = 1,nucl_num 
  nuclei_part += nucl_charge(i) * nucl_coord(i,1) 
 enddo
 print*,'nuclei   part for x_dipole = ',nuclei_part

 do istate = 1, N_states
  x_dipole_moment(istate) += nuclei_part
 enddo

 call cpu_time(cpu1)
 print*,'Time for the x_dipole_moment integration :',cpu1-cpu0
END_PROVIDER

 BEGIN_PROVIDER [double precision, y_dipole_moment, (N_states)]
&BEGIN_PROVIDER [double precision, y_dipole_moment_diag, (N_states)]
&BEGIN_PROVIDER [double precision, y_dipole_moment_off_diag, (N_states)]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,i,j 
 double precision :: weight, r(3)  
 double precision :: cpu0,cpu1,nuclei_part

 call cpu_time(cpu0)
 y_dipole_moment = 0.d0
 y_dipole_moment_diag = 0.d0
 y_dipole_moment_off_diag= 0.d0
 do istate = 1, N_states
  do i = 1, mo_num  
   y_dipole_moment_diag(istate) += -(one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)) *  mo_dipole_y(i,i)
   do j = i+1, mo_num  
    y_dipole_moment_off_diag(istate) += - 2.d0 * (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_y(j,i) 
   enddo
  enddo
  y_dipole_moment(istate) = y_dipole_moment_diag(istate) + y_dipole_moment_off_diag(istate)
 enddo
 
 print*,'electron part for y_dipole = ',y_dipole_moment
 
 nuclei_part = 0.d0
 do i = 1,nucl_num 
  nuclei_part += nucl_charge(i) * nucl_coord(i,2) 
 enddo
 print*,'nuclei   part for y_dipole = ',nuclei_part

 do istate = 1, N_states
  y_dipole_moment(istate) += nuclei_part
 enddo

 call cpu_time(cpu1)
 print*,'Time for the y_dipole_moment integration :',cpu1-cpu0
END_PROVIDER

 BEGIN_PROVIDER [double precision, z_dipole_moment, (N_states)]
&BEGIN_PROVIDER [double precision, z_dipole_moment_diag, (N_states)]
&BEGIN_PROVIDER [double precision, z_dipole_moment_off_diag, (N_states)]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,i,j 
 double precision :: weight, r(3)  
 double precision :: cpu0,cpu1,nuclei_part

 call cpu_time(cpu0)
 z_dipole_moment = 0.d0
 z_dipole_moment_diag = 0.d0
 z_dipole_moment_off_diag= 0.d0
 do istate = 1, N_states
  do i = 1, mo_num  
   z_dipole_moment_diag(istate) += -(one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)) *  mo_dipole_z(i,i)
   do j = i+1, mo_num  
    z_dipole_moment_off_diag(istate) += - 2.d0 * (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_z(j,i) 
   enddo
  enddo
  z_dipole_moment(istate) = z_dipole_moment_diag(istate) + z_dipole_moment_off_diag(istate)
 enddo
 
 print*,'electron part for z_dipole = ',z_dipole_moment
 
 nuclei_part = 0.d0
 do i = 1,nucl_num 
  nuclei_part += nucl_charge(i) * nucl_coord(i,3) 
 enddo
 print*,'nuclei   part for z_dipole = ',nuclei_part

 do istate = 1, N_states
  z_dipole_moment(istate) += nuclei_part
 enddo

 call cpu_time(cpu1)
 print*,'Time for the z_dipole_moment integration :',cpu1-cpu0
END_PROVIDER

 BEGIN_PROVIDER [double precision, total_dipole_moment, (N_states)]
 implicit none
 BEGIN_DOC
 ! Total dipole moment sqrt(x_dipole**2 + y_dipole**2 + z_dipole**2)
 END_DOC
 total_dipole_moment = dsqrt(x_dipole_moment**2 + y_dipole_moment**2 + z_dipole_moment**2)
END_PROVIDER


 subroutine print_x_dipole_moment_only
 implicit none
  print*, ''
  print*, ''
  print*,  '****************************************'
  print*, 'x_dipole_moment                  = ',x_dipole_moment
  print*,  '****************************************'
 end

 subroutine print_y_dipole_moment_only
 implicit none
  print*, ''
  print*, ''
  print*,  '****************************************'
  print*, 'y_dipole_moment                  = ',y_dipole_moment
  print*,  '****************************************'
 end

 subroutine print_z_dipole_moment_only
 implicit none
  print*, ''
  print*, ''
  print*,  '****************************************'
  print*, 'z_dipole_moment                  = ',z_dipole_moment
  print*,  '****************************************'
 end

 subroutine print_tot_dipole_moment_only
 implicit none
  print*, ''
  print*, ''
  print*,  '****************************************'
  print*, 'total_dipole_moment                  = ',total_dipole_moment
  print*,  '****************************************'
 end
