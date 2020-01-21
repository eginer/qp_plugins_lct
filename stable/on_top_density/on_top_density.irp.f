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
 double precision :: zmax,dz,r(3),on_top_in_r,total_density
 integer :: nz,i,istate
 zmax = 5.d0
 nz = 1000
 dz = zmax / dble(nz)
 r(:) = 0.d0 
 r(3) = -zmax * 0.5d0 
 istate = 1
 do i = 1, nz
  call give_on_top_in_r_one_state(r,istate,on_top_in_r)
  call give_cas_density_in_r(r,total_density)
  write(33,*)r(3),on_top_in_r,total_density
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
   prov_dm = core_inact_act_on_top_of_r_new(i_point,istate)
   call give_on_top_in_r_one_state(r,istate,on_top_in_r)
   accu(istate) += dabs(prov_dm - on_top_in_r) * final_weight_at_r_vector(i_point)
   accu_2(istate) += prov_dm * final_weight_at_r_vector(i_point)
  enddo
 enddo
 print*,'accu   = ',accu
 print*,'accu_2 = ',accu_2
end
