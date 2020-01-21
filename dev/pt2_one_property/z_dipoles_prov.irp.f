
 BEGIN_PROVIDER [double precision, z_dipole_moment, (N_states)]
&BEGIN_PROVIDER [double precision, z_dipole_moment_diag, (N_states)]
&BEGIN_PROVIDER [double precision, z_dipole_moment_off_diag, (N_states)]
&BEGIN_PROVIDER [double precision, elec_z_dipole_moment, (N_states)]
&BEGIN_PROVIDER [double precision, nucl_z_dipole_moment]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,i,j 
 double precision :: weight, r(3)  
 double precision :: cpu0,cpu1
 double precision :: elec_num_dm(N_states)

 call cpu_time(cpu0)
 z_dipole_moment = 0.d0
 z_dipole_moment_diag = 0.d0
 z_dipole_moment_off_diag= 0.d0
 elec_num_dm = 0.d0
 do istate = 1, N_states
  do i = 1, mo_num  
   elec_num_dm(istate) += one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)
   z_dipole_moment_diag(istate) += -(one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)) *  mo_dipole_z(i,i)
   do j = i+1, mo_num  
    z_dipole_moment_off_diag(istate) += - 2.d0 * (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_z(j,i) 
   enddo
  enddo
  z_dipole_moment(istate) = z_dipole_moment_diag(istate) + z_dipole_moment_off_diag(istate)
 enddo
 elec_z_dipole_moment = z_dipole_moment 
 print*,'electron part for z_dipole = ',z_dipole_moment
 
 nucl_z_dipole_moment = 0.d0
 do i = 1,nucl_num 
  nucl_z_dipole_moment += nucl_charge(i) * nucl_coord(i,3) 
 enddo
 print*,'nuclei   part for z_dipole = ',nucl_z_dipole_moment
 print*,'elec_num_dm = ',elec_num_dm

 do istate = 1, N_states
  z_dipole_moment(istate) += nucl_z_dipole_moment
 enddo

 call cpu_time(cpu1)
 print*,'Time for the z_dipole_moment integration :',cpu1-cpu0
END_PROVIDER


