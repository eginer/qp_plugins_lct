BEGIN_PROVIDER [integer, occ_two_e_psi, (2, N_det)]
 implicit none
 use bitmasks
 integer :: i,j,n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 if(elec_num.gt.2)then
  print*,'this provider works only for two electron WF'
  stop
 endif
 do i = 1, N_det
  call bitstring_to_list_ab(psi_det(1,1,i), occ, n_occ_ab, N_int)
  occ_two_e_psi(1,i) = occ(1,1)
  occ_two_e_psi(2,i) = occ(1,2)
 enddo

END_PROVIDER 

subroutine get_two_e_psi_at_r1r2(r1,r2,psi)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: psi(N_states)
 double precision, allocatable :: mos_array_r2(:),mos_array_r1(:)
 allocate(mos_array_r1(mo_num),mos_array_r2(mo_num))
 call give_all_mos_at_r(r1,mos_array_r1)  
 call give_all_mos_at_r(r2,mos_array_r2)  
 integer :: i,istate,i_up,i_down
 psi = 0.d0
 do istate = 1, N_states
  do i = 1, N_det
   i_up   = occ_two_e_psi(1,i)
   i_down = occ_two_e_psi(2,i)
   psi(istate) += mos_array_r1(i_up) * mos_array_r2(i_down) * psi_coef(i,istate)
  enddo
 enddo

end
