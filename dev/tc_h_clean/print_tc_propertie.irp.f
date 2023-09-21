program print_tc_properties
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf 

! call prepare_psi_coef
 call print_prop

end

!subroutine prepare_psi_coef
! implicit none
! double precision, allocatable :: reigvec_tc_tmp(:),leigvec_tc_tmp(:)
! integer :: j
! allocate(reigvec_tc_tmp(N_det),leigvec_tc_tmp(N_det))
!  do j = 1, N_det
!   reigvec_tc_tmp(j) = reigvec_tc(j,1)
!   leigvec_tc_tmp(j) = leigvec_tc(j,1)
!  enddo 
!  call set_psi_det_to_rl_eigv(reigvec_tc_tmp,leigvec_tc_tmp)
!end
!
!subroutine set_psi_det_to_rl_eigv(reigvec_tc_tmp,leigvec_tc_tmp)
! implicit none
! double precision, intent(in) :: reigvec_tc_tmp(N_det),leigvec_tc_tmp(N_det)
! integer :: j
! N_states = 2
! touch N_states 
! do j = 1, N_det
!  psi_coef(j,1) = reigvec_tc_tmp(j)
!  psi_coef(j,2) = leigvec_tc_tmp(j)
! enddo
! touch psi_coef 
!end 

subroutine print_prop
 implicit none
 integer :: i
 double precision :: accu 
 accu = 0.D0
 do i = 1,N_det
  print*,'psi_coef',i,psi_coef(i,1),psi_coef(i,2)
  accu += psi_coef(i,1) * psi_coef(i,2)
 enddo
 print*,'accu               = ',accu
  print*,'one_e_tm_mo_norm   = ',one_e_tm_mo_norm
  print*,'*********************'
  print*,'x_dipole_tc_moment = ',x_dipole_tc_moment
  print*,'y_dipole_tc_moment = ',y_dipole_tc_moment
  print*,'z_dipole_tc_moment = ',z_dipole_tc_moment
end
