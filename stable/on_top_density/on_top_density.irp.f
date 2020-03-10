program on_top_density
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
!  call routine
  call print_on_top
end

subroutine print_on_top
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

subroutine routine
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
   prov_dm = core_inact_act_on_top_of_r(i_point,istate)
   call give_on_top_in_r_one_state(r,istate,on_top_in_r)
   accu(istate) += dabs(prov_dm - on_top_in_r) * final_weight_at_r_vector(i_point)
   accu_2(istate) += prov_dm * final_weight_at_r_vector(i_point)
  enddo
 enddo
 print*,'accu   = ',accu
 print*,'accu_2 = ',accu_2
end
