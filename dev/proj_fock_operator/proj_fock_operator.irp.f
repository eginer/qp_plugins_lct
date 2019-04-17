program proj_fock_operator
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print*,'mo_num = ',mo_num
  call pouet
end

subroutine pouet
  implicit none
  integer :: k,j,nx
  nx = 10000
  double precision :: dx,r(3),rmax
  double precision :: kin_proj_op_in_r_on_phi_k, v_ne_proj_op_in_r_on_phi_k
  double precision :: kin_op_in_r_on_phi_k, v_ne_op_in_r_on_phi_k
  double precision :: proj_kin,proj_v_ne,kin,v_ne
  double precision :: mos_array(mo_num)
  double precision :: mo_two_e_integral
  rmax = 5.d0
  dx =  rmax/(dble(nx))
  r = 0.d0

  k = 1 ! which MO you try
  print*,'epsilon k   = ',Fock_matrix_diag_mo(k)
  print*,'h       k   = ',mo_kinetic_integrals(k,k) + mo_integrals_n_e(k,k) + mo_F_two_e_operator(k,k)

  r(3) = -rmax+1.d-4
  do j = 1, 2 * nx
   call give_all_mos_at_r(r,mos_array) 

   proj_kin = kin_proj_op_in_r_on_phi_k(k,r)/mos_array(k) 
   kin = kin_op_in_r_on_phi_k(k,r)/mos_array(k) 

   proj_v_ne = v_ne_proj_op_in_r_on_phi_k(k,r) / mos_array(k) 
   v_ne = v_ne_op_in_r_on_phi_k(k,r) / mos_array(k) 

   write(33,'(100(F16.10,X))')r(3),kin,proj_kin,v_ne,proj_v_ne,mos_array(k)
   r(3) += dx 
  enddo
end
