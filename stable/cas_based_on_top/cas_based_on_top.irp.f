program cas_based_on_top_density
  implicit none
  BEGIN_DOC
! TODO : Small example to use the different quantities in this plugin
  END_DOC

  !! You force QP2 to read the wave function in the EZFIO folder
  !! It is assumed that all Slater determinants in the wave function 
  !! belongs to an active space defined by core, inactive and active list of orbitals 
  read_wf = .True.
  touch read_wf
  call routine_test_cas_based_on_top_density
end

subroutine routine_test_cas_based_on_top_density
 implicit none
 integer :: i_point,istate
 double precision :: r(3),prov_dm,on_top_in_r
 double precision :: accu(N_states),accu_2(N_states)
 accu = 0.d0
 accu_2 = 0.d0
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  do istate = 1, N_states
   ! provider to get the on-top on the becke grid 
   prov_dm = total_cas_on_top_density(i_point,istate)
   ! subroutine to get the on-top in any points in space
   call give_on_top_in_r_one_state(r,istate,on_top_in_r)
   accu(istate) += dabs(prov_dm - on_top_in_r) * final_weight_at_r_vector(i_point)
   accu_2(istate) += prov_dm * final_weight_at_r_vector(i_point)
  enddo
 enddo
 print*,'difference between provider and routine = ',accu
 print*,'integral of the on-top                  = ',accu_2
 print*,'integral_on_top                         = ',integral_on_top(:)
end

subroutine write_on_top_in_real_space
 implicit none
 double precision :: zmax,dz,r(3),on_top_in_r,total_density,zcenter,dist

 integer :: nz,i,istate
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename
 output=trim(ezfio_filename)//'.on_top'
 print*,'output = ',trim(output)                                                                                         
 i_unit_output = getUnitAndOpen(output,'w')


 zmax = 2.0d0
 print*,'nucl_coord(1,3) = ',nucl_coord(1,3)
 print*,'nucl_coord(2,3) = ',nucl_coord(2,3)
 dist = dabs(nucl_coord(1,3) - nucl_coord(2,3))
 zmax += dist 
 zcenter = (nucl_coord(1,3) + nucl_coord(2,3))*0.5d0
 print*,'zcenter = ',zcenter
 print*,'zmax    = ',zmax
 nz = 1000
 dz = zmax / dble(nz)
 r(:) = 0.d0 
 r(3) = zcenter -zmax * 0.5d0 
 print*,'r(3)    = ',r(3)
 istate = 1
 do i = 1, nz
  call give_on_top_in_r_one_state(r,istate,on_top_in_r)
  call give_cas_density_in_r(r,total_density)
  write(i_unit_output,*)r(3),on_top_in_r,total_density
  r(3) += dz
 enddo


end

