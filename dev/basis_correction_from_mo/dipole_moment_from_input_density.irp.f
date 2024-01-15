 BEGIN_PROVIDER [ double precision, z_dipole_moment_from_input_density]
&BEGIN_PROVIDER [ double precision, y_dipole_moment_from_input_density]
&BEGIN_PROVIDER [ double precision, x_dipole_moment_from_input_density]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,i,j 
 double precision :: weight, r(3)  
 double precision :: cpu0,cpu1,nuclei_part_z,nuclei_part_y,nuclei_part_x

 provide data_one_e_dm_alpha_mo
 provide data_one_e_dm_beta_mo

! call cpu_time(cpu0)
 z_dipole_moment_from_input_density = 0.d0
 y_dipole_moment_from_input_density = 0.d0
 x_dipole_moment_from_input_density = 0.d0
 istate=1
 do i = 1, mo_num  
   z_dipole_moment_from_input_density += -(data_one_e_dm_alpha_mo(i,i,istate)+data_one_e_dm_beta_mo(i,i,istate)) *  mo_dipole_z(i,i)
   y_dipole_moment_from_input_density += -(data_one_e_dm_alpha_mo(i,i,istate)+data_one_e_dm_beta_mo(i,i,istate)) *  mo_dipole_y(i,i)
   x_dipole_moment_from_input_density += -(data_one_e_dm_alpha_mo(i,i,istate)+data_one_e_dm_beta_mo(i,i,istate)) *  mo_dipole_x(i,i)
   do j = i+1, mo_num  
    z_dipole_moment_from_input_density += - 2.d0 * (data_one_e_dm_alpha_mo(j,i,istate)+data_one_e_dm_beta_mo(j,i,istate)) *  mo_dipole_z(j,i) 
    y_dipole_moment_from_input_density += - 2.d0 * (data_one_e_dm_alpha_mo(j,i,istate)+data_one_e_dm_beta_mo(j,i,istate)) *  mo_dipole_y(j,i) 
    x_dipole_moment_from_input_density += - 2.d0 * (data_one_e_dm_alpha_mo(j,i,istate)+data_one_e_dm_beta_mo(j,i,istate)) *  mo_dipole_x(j,i) 
   enddo
  enddo
 
 print*,'electron part for z_dipole = ', z_dipole_moment_from_input_density
 print*,'electron part for y_dipole = ', y_dipole_moment_from_input_density
 print*,'electron part for x_dipole = ', x_dipole_moment_from_input_density
 
 nuclei_part_z = 0.d0
 nuclei_part_y = 0.d0
 nuclei_part_x = 0.d0
 do i = 1,nucl_num 
  nuclei_part_z += nucl_charge(i) * nucl_coord(i,3) 
  nuclei_part_y += nucl_charge(i) * nucl_coord(i,2) 
  nuclei_part_x += nucl_charge(i) * nucl_coord(i,1) 
 enddo
 print*,'nuclei   part for z_dipole = ',nuclei_part_z
 print*,'nuclei   part for y_dipole = ',nuclei_part_y
 print*,'nuclei   part for x_dipole = ',nuclei_part_x

 z_dipole_moment_from_input_density += nuclei_part_z
 y_dipole_moment_from_input_density += nuclei_part_y
 x_dipole_moment_from_input_density += nuclei_part_x

END_PROVIDER
