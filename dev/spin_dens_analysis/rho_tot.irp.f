


BEGIN_PROVIDER [double precision, integrated_rho_tot_all_points, (N_z_pts)]
 BEGIN_DOC
! 
! integrated_rho(alpha,z) - integrated_rho(beta,z) for all the z points 
! chosen
!
 END_DOC
 implicit none
 integer :: i,j,k,l,i_z,h
 double precision :: z,function_integrated_delta_rho,c_k,c_j,n_i_h,accu
 integrated_rho_tot_all_points = 0.d0
 do i_z = 1, N_z_pts
  do i = 1, mo_num
    do h = 1, mo_num
     n_i_h = one_e_dm_mo_alpha_average(i,h)+one_e_dm_mo_beta_average(i,h)
     if(dabs(n_i_h).lt.1.d-10)cycle
     do j = 1, ao_num
      c_j = mo_coef(j,i)   ! coefficient of the ith MO on the jth AO
      do k = 1, ao_num
       c_k = mo_coef(k,h)   ! coefficient of the hth MO on the kth AO
       integrated_rho_tot_all_points(i_z) += c_k * c_j * n_i_h *  ao_integrated_x_y_all_points(j,k,i_z)
      enddo
     enddo
    enddo
  enddo
 enddo
!!$OMP END PARALLEL DO

 z = z_min
 accu = 0.d0
 do i = 1, N_z_pts
  accu += integrated_rho_tot_all_points(i)
  write(i_unit_integrated_rho_tot,*)z,integrated_rho_tot_all_points(i),accu
  z += delta_z
 enddo
 print*,'sum of integrated_delta_rho = ',accu

END_PROVIDER


BEGIN_PROVIDER [integer, i_unit_integrated_rho_tot]
 implicit none
 BEGIN_DOC
! fortran unit for the writing of the integrated delta_rho
 END_DOC
 integer :: getUnitAndOpen
 character*(128) :: output_i_unit_integrated_rho_tot
 output_i_unit_integrated_rho_tot=trim(ezfio_filename)//'.rho_tot'
 i_unit_integrated_rho_tot= getUnitAndOpen(output_i_unit_integrated_rho_tot,'w')

END_PROVIDER


