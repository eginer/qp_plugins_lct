program plot_density
  implicit none
  read_wf = .True.
  touch read_wf
  call routine
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
end

subroutine routine
 implicit none
 double precision :: u(3),norm_u,du(3),accu,dx,r(3),dist,dm_a,dm_b,on_top_in_r
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.density'                                                                                                    
 i_unit_output = getUnitAndOpen(output,'w')

 integer :: i,nx,istate
 integer :: i1,i2
 ! vector that goes from i2 to i1 
 i2 = 3 ! oxygen
 i1 = 2 ! hydrogen 
 ! u = i2 - i1
 ! i2 - u = i1 
 u(:) = nucl_coord_transp(:,i2) - nucl_coord_transp(:,i1)
 norm_u = 0.d0
 do i = 1, 3
  norm_u += u(i)**2.d0
 enddo
 norm_u = dsqrt(norm_u)
 print*,'norm_u = ',norm_u
 print*,'nucl_dist',nucl_dist(i1,i2)
 du(:) = u(:)/norm_u
 nx = 1000
! dx = (norm_u+1.5d0)/nx
 dx = (norm_u)/nx
 du(:) = du(:) * dx
 accu = 0.d0
 print*,'dx  = ',dx

 r(:) = nucl_coord_transp(:,i2)
 istate = 1
 write(i_unit_output,*)' #  distance from O  on top           density'
 do i = 1, nx
  call give_on_top_in_r_one_state(r,istate,on_top_in_r)
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  dist  = dsqrt((r(1) - nucl_coord_transp(1,i2))**2.d0)
  dist += dsqrt((r(2) - nucl_coord_transp(2,i2))**2.d0)
  dist += dsqrt((r(3) - nucl_coord_transp(3,i2))**2.d0)
  write(i_unit_output,'(100(F16.10,X))')dist,2.d0 * on_top_in_r,dm_a+dm_b
  r(:) -= du(:)
 enddo
 print*,'r = '
 print*,r
 print*,'r(1)'
 print*, nucl_coord_transp(:,i1)
 print*,'integral_on_top(istate) = ',2.d0 * integral_on_top(istate)
 print*,'v_ne_psi_energy(istate) = ',v_ne_psi_energy(istate)
 


end
