program write_ec_in_z
 implicit none
 read_wf = .True.
 touch read_wf
 mu_of_r_potential = "hf_coallescence"
 touch mu_of_r_potential 
 no_core_density = "pouet"
 touch no_core_density
 call routine_2
!call pouet
end

subroutine pouet
 implicit none
 print*,'elec_beta_num_no_core_grid_cyl  = ',elec_beta_num_no_core_grid_cyl
 print*,'elec_alpha_num_no_core_grid_cyl = ',elec_alpha_num_no_core_grid_cyl
 print*,'elec_beta_num_grid_cyl          = ', elec_beta_num_grid_cyl
 print*,'elec_alpha_num_grid_cyl         = ',elec_alpha_num_grid_cyl
end

subroutine routine
 implicit none
 integer :: i,j
 provide n_z_points_grid 
 print*,'n_z_points_grid     = ',n_z_points_grid
 print*,'e_c_lda_val_hf_sum      = ',e_c_lda_val_hf_sum
 print*,'e_c_lda_ful_hf_sum      = ',e_c_lda_ful_hf_sum
 print*,'Energy_c_md_mu_of_r_LDA = ',Energy_c_md_mu_of_r_LDA

 provide ezfio_filename
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen

 character*(128) :: output2
 integer :: i_unit_output2

 output=trim(ezfio_filename)//'.'//trim("ec_in_z")
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 output2=trim(ezfio_filename)//'.'//trim("ec_in_z_bis")
 output2=trim(output2)
 print*,'output2 = ',trim(output2)
 i_unit_output2 = getUnitAndOpen(output2,'w')


 do i = 1, n_z_points_grid
  write(i_unit_output,'(100(F16.10,X))')z_points_grid(i),e_c_lda_ful_hf_z_tab(i),e_c_lda_val_hf_z_tab(i),e_c_lda_ful_hf_z_sum(i),e_c_lda_val_hf_z_sum(i)
 enddo
 
 double precision :: z_min,dz,z_max,z,zm,zp,df,dzm,dzp
 integer :: nz,izm,izp,jm,jp,jz
 z_min = nucl_coord(1,3) - 2.5d0
 z_max = nucl_coord(1,3) + 7.5d0
 nz = 1000
 dz = (z_max -z_min)/dble(nz)
 z = z_min
 do i = 1, nz
  zm = z - dz 
  zp = z + dz 
  do j = 1, n_z_points_grid
   if(z_points_grid(j).gt.zm)then
    jm = max(j -1,1)
    exit
   endif
  enddo

  do j = 1, n_z_points_grid
   if(z_points_grid(j).gt.z.and.z_points_grid(j).ne.z_points_grid(jm))then
    jz = max(j -1,1)
    exit
   endif
  enddo

  do j = 1, n_z_points_grid
   if(z_points_grid(j).gt.zp.and.z_points_grid(j).ne.z_points_grid(jz))then
    jp = max(j -1,1)
    exit
   endif
  enddo
  dzp = z_points_grid(jp) - z_points_grid(jz)
  dzm = z_points_grid(jz) - z_points_grid(jm)
  df = 0.5d0 * ((e_c_lda_ful_hf_z_sum(jz) - e_c_lda_ful_hf_z_sum(jm)) / dzm +  (e_c_lda_ful_hf_z_sum(jp) - e_c_lda_ful_hf_z_sum(jz)) / dzp)
  
  write(i_unit_output2,'(100(F16.10,X))')z,df
  z += dz
 enddo

end


subroutine routine_2
 implicit none
 integer :: i,j
 provide n_z_points_grid 
 print*,''
 print*,'Energy_c_md_mu_of_r_LDA       = ',Energy_c_md_mu_of_r_LDA
 print*,'e_c_lda_hf_ful_grid_cyl_z_sum = ',e_c_lda_hf_ful_grid_cyl_z_sum
 print*,'elec_alpha_num_grid_cyl       = ',elec_alpha_num_grid_cyl
 print*,'elec_beta_num_grid_cyl        = ',elec_beta_num_grid_cyl
 print*,'elec_alpha_num                = ',elec_alpha_num
 print*,'elec_beta_num                 = ',elec_beta_num
 print*,''
 print*,'Valence quantities '
 print*,''
 print*,'Energy_c_md_mu_of_r_LDA_val   = ',Energy_c_md_mu_of_r_LDA_val
 print*,'e_c_lda_hf_val_grid_cyl_z_sum = ',e_c_lda_hf_val_grid_cyl_z_sum
 print*,'elec_alpha_num_no_core_grid   = ',elec_alpha_num_no_core_grid_cyl
 print*,'elec_beta_num_grid_cyl        = ',elec_beta_num_no_core_grid_cyl
 print*,'elec_alpha_num - n_core_orb   = ',elec_alpha_num - n_core_orb
 print*,'elec_beta_num  - n_core_orb   = ',elec_beta_num  - n_core_orb

 provide ezfio_filename
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen

 character*(128) :: output2
 integer :: i_unit_output2

 output=trim(ezfio_filename)//'.'//trim("ec_in_z")
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')
 output2=trim(ezfio_filename)//'.'//trim("ec_in_z_val")
 output2=trim(output2)
 print*,'output2 = ',trim(output2)
 i_unit_output2 = getUnitAndOpen(output2,'w')
 

 do i = 1, n_z_grid_cyl
  write(i_unit_output,'(100(F16.10,X))')z_points_grid_cylindr(i),e_c_lda_hf_ful_grid_cyl_z_tab(i),mu_hf_ful_grid_cyl_z_tab(i),dm_grid_cyl_z_tab(i)
  write(i_unit_output2,'(100(F16.10,X))')z_points_grid_cylindr(i),e_c_lda_hf_val_grid_cyl_z_tab(i),mu_hf_val_grid_cyl_z_tab(i),dm_no_core_grid_cyl_z_tab(i)
! write(i_unit_output ,'(100(F16.10,X))')z_points_grid_cylindr(i),e_c_lda_ful_sym(i)+e_c_lda_hf_ful_grid_cyl_z_tab(i),mu_hf_ful_grid_cyl_z_tab(i),e_c_lda_ful_sym(i),dm_ful_sym(i)+dm_grid_cyl_z_tab(i),dm_grid_cyl_z_tab(i)
! write(i_unit_output2,'(100(F16.10,X))')z_points_grid_cylindr(i),e_c_lda_val_sym(i)+e_c_lda_hf_val_grid_cyl_z_tab(i),mu_hf_val_grid_cyl_z_tab(i),e_c_lda_val_sym(i),dm_val_sym(i)+dm_no_core_grid_cyl_z_tab(i),dm_no_core_grid_cyl_z_tab(i)
 enddo

end
