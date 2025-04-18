program print_on_top_z
 implicit none
 io_two_body_rdm_ab = "Read"
 touch io_two_body_rdm_ab 
 no_core_density = .True.
 touch no_core_density
 call write_on_top_in_real_space 
 call integrated_on_top
end


subroutine write_on_top_in_real_space
 implicit none
 BEGIN_DOC
! This routines is a simple example of how to plot the on-top pair density on a simple 1D grid 
 END_DOC
 double precision :: zmax,dz,r(3),on_top_in_r,total_density,zcenter,dist
 double precision :: core_dens, inact_dens,act_dens(2,1)

 integer :: nz,i,istate
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename
 output=trim(ezfio_filename)//'.on_top'
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')


 zmax = 5.0d0
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

  write(i_unit_output,*)" z,    on-top(z)"
 do i = 1, nz
  call give_on_top_in_r_one_state(r,istate,on_top_in_r)
!  call give_cas_density_in_r(core_dens,inact_dens,act_dens,total_density,r)
  write(i_unit_output,*)r(3),on_top_in_r
  r(3) += dz
 enddo

end

subroutine integrated_on_top
 implicit none
 integer :: ipoint, istate
 double precision :: weight,accu
 istate = 1
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  accu += weight * total_cas_on_top_density(ipoint,istate)
 enddo
 print*,'integrated on-top pair density : ',accu
end
