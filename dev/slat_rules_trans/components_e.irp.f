
subroutine print_e_comp_transcorr
 implicit none
 print*,'tot_htrans_expect     = ',tot_htrans_expect(1) + nuclear_repulsion
 print*,'tot_htrans_expect_bis = ',tot_htrans_expect_bis(1) + nuclear_repulsion
 print*,'erf_expect+eff_expect = ',erf_expect(1)+eff_expect(1)
 print*,'deriv_expect          = ',deriv_expect(1)
 print*,'three-body            = ',hthree_body_expect(1)
 print*,'hcore_expect          = ',hcore_expect(1)

end

 BEGIN_PROVIDER [double precision, tot_htrans_expect_bis, (N_states)]
&BEGIN_PROVIDER [double precision, tot_htrans_expect, (N_states)]
&BEGIN_PROVIDER [double precision, erf_expect, (N_states)]
&BEGIN_PROVIDER [double precision, eff_expect, (N_states)]
&BEGIN_PROVIDER [double precision, deriv_expect, (N_states)]
&BEGIN_PROVIDER [double precision, hcore_expect, (N_states)]
&BEGIN_PROVIDER [double precision, hthree_body_expect, (N_states)]
 implicit none
 integer :: i,j
 erf_expect = 0.d0
 eff_expect = 0.d0
 deriv_expect = 0.d0
 hcore_expect = 0.d0
 hthree_body_expect = 0.d0
 tot_htrans_expect_bis = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   tot_htrans_expect_bis(1) += htilde_matrix_elmt(j,i) * reigvec_trans(j,1) * reigvec_trans(i,1)
   erf_expect(1)   += reigvec_trans(j,1) * reigvec_trans(i,1) * htilde_matrix_elmt_erf(j,i)
   eff_expect(1)   += reigvec_trans(j,1) * reigvec_trans(i,1) * htilde_matrix_elmt_eff(j,i)
   deriv_expect(1) += reigvec_trans(j,1) * reigvec_trans(i,1) * htilde_matrix_elmt_deriv(j,i)
   hcore_expect(1) += reigvec_trans(j,1) * reigvec_trans(i,1) * htilde_matrix_elmt_hcore(j,i)
   hthree_body_expect(1) += reigvec_trans(j,1) * reigvec_trans(i,1) * htilde_matrix_elmt_hthree(j,i)
  enddo
 enddo
 tot_htrans_expect(1) = erf_expect(1) + eff_expect(1) + deriv_expect(1) + hcore_expect(1) + hthree_body_expect(1)
! tot_htrans_expect(1) *= 1.d0/dsqrt(reigvec_trans_norm(1))

END_PROVIDER 
