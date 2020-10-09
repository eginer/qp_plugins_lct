subroutine get_two_e_psi_general_at_r1r2(r1,r2,ndet,nstates,occ_psi,coefs,psi)
 implicit none
 double precision, intent(in) :: r1(3),r2(3),coefs(ndet,nstates)
 integer, intent(in) :: occ_psi(2,ndet),ndet,nstates
 double precision, intent(out):: psi(N_states)
 double precision, allocatable :: mos_array_r2(:),mos_array_r1(:)
 allocate(mos_array_r1(mo_num),mos_array_r2(mo_num))
 call give_all_mos_at_r(r1,mos_array_r1)  
 call give_all_mos_at_r(r2,mos_array_r2)  
 integer :: i,istate,i_up,i_down
 psi = 0.d0
 do istate = 1, nstates
  do i = 1, ndet
   i_up   = occ_psi(1,i)
   i_down = occ_psi(2,i)
   psi(istate) += mos_array_r1(i_up) * mos_array_r2(i_down) * coefs(i,istate)
  enddo
 enddo
end
