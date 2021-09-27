

subroutine get_delta_no_store(psidet,psicoef,ndet,delta)
 implicit none
 use bitmasks
 double precision, intent(in)   :: psicoef(ndet)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet
 double precision, intent(out)  :: delta(ndet) 
 double precision :: hij,htilde_ij,delta_mat,hmono,heff,hderiv,hthree
 integer :: i,j
 PROVIDE scalar_mu_r_pot_physicist_mo deriv_mu_r_pot_physicist_mo
 PROVIDE three_body_3_index three_body_3_index_exch_12 three_body_3_index_exch_13 three_body_3_index_exch_23
 PROVIDE three_body_5_index three_body_5_index_exch_13 three_body_5_index_exch_32
 PROVIDE three_body_4_index three_body_4_index_exch_12 three_body_4_index_exch_12_part

 delta = 0.d0
 i=1
 j=1
 call htilde_mat(psi_det(1,1,i),psi_det(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
 call  i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
 !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(dynamic,8) &
 !$OMP PRIVATE(i,j,delta_mat,hmono,heff,hderiv,hthree,htilde_ij,hij)
  do i = 1, N_det
   do j = 1, N_det
    call htilde_mat(psi_det(1,1,i),psi_det(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    call  i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
    delta_mat = htilde_ij - hij 
    delta(i) = delta(i) + psicoef(j) * delta_mat
   enddo
  enddo
 !$OMP END PARALLEL DO

end
