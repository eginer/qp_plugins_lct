 BEGIN_PROVIDER [integer, spin_dens_coord]
 implicit none
 BEGIN_DOC
 !coordinate on which you are going to plot the spin density
 !and integrate over the ohters
 !spin_dens_coord = 1  === X
 !spin_dens_coord = 2  === Y
 !spin_dens_coord = 3  === Z
 END_DOC
 spin_dens_coord = 3
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, delta_z]
&BEGIN_PROVIDER [double precision, z_min]
&BEGIN_PROVIDER [double precision, z_max]
 implicit none
!cucl2
! z_min = -5.d0
! z_max = 8.d0
!BH
z_min = -5.d0
z_max = 5.d0
 delta_z = 0.05d0
END_PROVIDER

BEGIN_PROVIDER [integer, N_z_pts]
 implicit none
 N_z_pts = int( (z_max - z_min)/delta_z )
 print*,'N_z_pts = ',N_z_pts
END_PROVIDER


BEGIN_PROVIDER [double precision, integrated_delta_rho_all_points, (N_z_pts)]
 BEGIN_DOC
! 
! integrated_rho(alpha,z) + integrated_rho(beta,z) for all the z points 
! chosen
!
 END_DOC
 implicit none
 integer :: i,j,k,l,i_z,h
 double precision :: z,function_integrated_delta_rho,c_k,c_j,n_i_h,accu
 integrated_delta_rho_all_points = 0.d0
 !$OMP PARALLEL DO DEFAULT(none) &
 !$OMP PRIVATE(i,h,j,k,c_j,c_k,n_i_h,i_z) &
 !$OMP SHARED(mo_num,ao_num,mo_coef, &
 !$OMP   ao_integrated_x_y_all_points,one_e_spin_density_mo,integrated_delta_rho_all_points,N_z_pts)          
 do i_z = 1, N_z_pts
  do i = 1, mo_num
    do h = 1, mo_num
     n_i_h = one_e_spin_density_mo(i,h)
     if(dabs(n_i_h).lt.1.d-10)cycle
     do j = 1, ao_num
      c_j = mo_coef(j,i)   ! coefficient of the ith MO on the jth AO
      do k = 1, ao_num
       c_k = mo_coef(k,h)   ! coefficient of the hth MO on the kth AO
       integrated_delta_rho_all_points(i_z) += c_k * c_j * n_i_h *  ao_integrated_x_y_all_points(j,k,i_z)
      enddo
     enddo
    enddo
  enddo
 enddo
 !$OMP END PARALLEL DO

 z = z_min
 accu = 0.d0
 do i = 1, N_z_pts
  accu += integrated_delta_rho_all_points(i)
  write(i_unit_integrated_delta_rho,*)z,integrated_delta_rho_all_points(i),accu
  z += delta_z
 enddo
 print*,'sum of integrated_delta_rho = ',accu

END_PROVIDER


!BEGIN_PROVIDER [double precision, integrated_rho_tot_all_points, (N_z_pts)]
! BEGIN_DOC
! 
! integrated_rho(alpha,z) + integrated_rho(beta,z) for all the z points 
! chosen

! END_DOC
! implicit none
! integer :: i,j,k,l,i_z,h
! double precision :: z,function_integrated_delta_rho,c_k,c_j,n_i_h,accu, accu_2
! integrated_rho_tot_all_points = 0.d0
! do i_z = 1, N_z_pts
!  do i = 1, mo_num
!    do h = 1, mo_num
!     n_i_h = one_e_dm_mo_alpha_average(i,h)+one_e_dm_mo_beta_average(i,h)
!     if(dabs(n_i_h).lt.1.d-10)cycle
!     do j = 1, ao_num
!      c_j = mo_coef(j,i)   ! coefficient of the ith MO on the jth AO
!      do k = 1, ao_num
!       c_k = mo_coef(k,h)   ! coefficient of the hth MO on the kth AO
!       integrated_rho_tot_all_points(i_z) += c_k * c_j * n_i_h *  ao_integrated_x_y_all_points(j,k,i_z)
!      enddo
!     enddo
!    enddo
!  enddo
! enddo
!!$OMP END PARALLEL DO
!
! z = z_min
! accu = 0.d0
! accu_2 = 0.d0
! do i = 1, N_z_pts
!  accu += integrated_rho_tot_all_points(i)
!  accu_2 +=  integrated_rho_tot_all_points(i)*z
!  write(i_unit_integrated_rho_tot,*)z,integrated_rho_tot_all_points(i),accu, accu_2
!  z += delta_z
! enddo
! print*,'sum of integrated_rho = ',accu
!
!END_PROVIDER





BEGIN_PROVIDER [ double precision, ao_integrated_x_y_all_points, (ao_num, ao_num, N_z_pts)]
 BEGIN_DOC
!  array of the overlap in x,y between the AO function and integrated between [z,z+dz] in the z axis
!  for all the z points that are given (N_z_pts)
 END_DOC
  implicit none
  integer :: i,j,n,l
  double precision :: f,accu
  integer :: dim1
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  integer :: i_z
  double precision :: z,SABpartial,accu_x,accu_y,accu_z
  dim1=100
 z = z_min
 do i_z = 1, N_z_pts
 !$OMP PARALLEL DO DEFAULT(none) &
 !$OMP PRIVATE(i,j,n,l,A_center,power_A,B_center,power_B,accu_z, &
 !$OMP  overlap_x,overlap_y,overlap_z,overlap,c,alpha,beta) &
 !$OMP SHARED(ao_num,nucl_coord,ao_nucl,ao_power,ao_prim_num,ao_expo_ordered_transp,ao_coef_normalized_ordered_transp, &
 !$OMP   ao_integrated_x_y_all_points,N_z_pts,dim1,i_z,z,delta_z,spin_dens_coord)           
  do j=1,ao_num
   A_center(1) = nucl_coord( ao_nucl(j), 1 )
   A_center(2) = nucl_coord( ao_nucl(j), 2 )
   A_center(3) = nucl_coord( ao_nucl(j), 3 )
   power_A(1)  = ao_power( j, 1 )
   power_A(2)  = ao_power( j, 2 )
   power_A(3)  = ao_power( j, 3 )
   do i= 1,ao_num
    B_center(1) = nucl_coord( ao_nucl(i), 1 )
    B_center(2) = nucl_coord( ao_nucl(i), 2 )
    B_center(3) = nucl_coord( ao_nucl(i), 3 )
    power_B(1)  = ao_power( i, 1 )
    power_B(2)  = ao_power( i, 2 )
    power_B(3)  = ao_power( i, 3 )

     accu_z = 0.d0
     do n = 1,ao_prim_num(j)
      alpha = ao_expo_ordered_transp(n,j)
      do l = 1, ao_prim_num(i)
       beta = ao_expo_ordered_transp(l,i)
       call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)

       c = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i) 
       if(spin_dens_coord ==1 )then
        accu_z += c* overlap_y * overlap_z * SABpartial(z,z+delta_z,A_center,B_center,power_A,power_B,alpha,beta,spin_dens_coord)
       else if (spin_dens_coord ==2 )then
        accu_z += c* overlap_x * overlap_z * SABpartial(z,z+delta_z,A_center,B_center,power_A,power_B,alpha,beta,spin_dens_coord)
       else if (spin_dens_coord ==3 )then
        accu_z += c* overlap_x * overlap_y * SABpartial(z,z+delta_z,A_center,B_center,power_A,power_B,alpha,beta,spin_dens_coord)
       endif
      enddo
     enddo
     ao_integrated_x_y_all_points(i,j,i_z) = accu_z
   enddo
  enddo
 !$OMP END PARALLEL DO
  z += delta_z
 enddo
END_PROVIDER

BEGIN_PROVIDER [integer, i_unit_integrated_delta_rho]
 implicit none
 BEGIN_DOC
! fortran unit for the writing of the integrated delta_rho
 END_DOC
 integer :: getUnitAndOpen
 character*(128) :: output_i_unit_integrated_delta_rho
 output_i_unit_integrated_delta_rho=trim(ezfio_filename)//'/delta_rho'
 i_unit_integrated_delta_rho= getUnitAndOpen(output_i_unit_integrated_delta_rho,'w')

END_PROVIDER



 BEGIN_PROVIDER [ double precision, exp_value_x_spread]
&BEGIN_PROVIDER [ double precision, exp_value_y_spread]
&BEGIN_PROVIDER [ double precision, exp_value_z_spread]
 implicit none
 integer :: i,j

 exp_value_x_spread = 0.d0
 exp_value_y_spread = 0.d0
 exp_value_z_spread = 0.d0

  do i = 1, ao_num 
   do j = 1, ao_num 
    exp_value_x_spread += ao_spread_x(j,i)*(one_e_dm_ao_alpha(j,i)+one_e_dm_ao_beta(j,i))
    exp_value_y_spread += ao_spread_y(j,i)*(one_e_dm_ao_alpha(j,i)+one_e_dm_ao_beta(j,i))
    exp_value_z_spread += ao_spread_z(j,i)*(one_e_dm_ao_alpha(j,i)+one_e_dm_ao_beta(j,i))
   enddo
  enddo
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, exp_value_x_dipole]
&BEGIN_PROVIDER [ double precision, exp_value_y_dipole]
&BEGIN_PROVIDER [ double precision, exp_value_z_dipole]
 implicit none
 integer :: i,j

 exp_value_x_dipole = 0.d0
 exp_value_y_dipole = 0.d0
 exp_value_z_dipole = 0.d0

  do i = 1, ao_num 
   do j = 1, ao_num 
    exp_value_x_dipole += ao_dipole_x(j,i)*(one_e_dm_ao_alpha(j,i)+one_e_dm_ao_beta(j,i))
    exp_value_y_dipole += ao_dipole_y(j,i)*(one_e_dm_ao_alpha(j,i)+one_e_dm_ao_beta(j,i))
    exp_value_z_dipole += ao_dipole_z(j,i)*(one_e_dm_ao_alpha(j,i)+one_e_dm_ao_beta(j,i))
   enddo
  enddo
 END_PROVIDER


 subroutine print_standard_deviation 
 implicit none
 print*,"standard_deviation_x  =",exp_value_x_spread - exp_value_x_dipole
 print*,"standard_deviation_y  =",exp_value_y_spread - exp_value_y_dipole
 print*,"standard_deviation_z  =",exp_value_z_spread - exp_value_z_dipole
 end
