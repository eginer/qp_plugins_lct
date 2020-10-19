program print_energy_comp_transcorr
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine

end

subroutine routine
 implicit none
 print*,'tot_htrans_expect     = ',tot_htrans_expect(1) + nuclear_repulsion
 print*,'erf_expect+eff_expect = ',erf_expect(1)+eff_expect(1)
 print*,'deriv_expect          = ',deriv_expect(1)
 print*,'hcore_expect          = ',hcore_expect(1)

end

 BEGIN_PROVIDER [double precision, erf_expect, (N_states)]
&BEGIN_PROVIDER [double precision, tot_htrans_expect, (N_states)]
&BEGIN_PROVIDER [double precision, eff_expect, (N_states)]
&BEGIN_PROVIDER [double precision, deriv_expect, (N_states)]
&BEGIN_PROVIDER [double precision, hcore_expect, (N_states)]
 implicit none
 integer :: i,j
 erf_expect = 0.d0
 eff_expect = 0.d0
 deriv_expect = 0.d0
 hcore_expect = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   erf_expect(1)   += psi_coef(j,1) * psi_coef(i,1) * htilde_matrix_elmt_erf(j,i)
   eff_expect(1)   += psi_coef(j,1) * psi_coef(i,1) * htilde_matrix_elmt_eff(j,i)
   deriv_expect(1) += psi_coef(j,1) * psi_coef(i,1) * htilde_matrix_elmt_deriv(j,i)
   hcore_expect(1) += psi_coef(j,1) * psi_coef(i,1) * htilde_matrix_elmt_hcore(j,i)
  enddo
 enddo
 tot_htrans_expect(1) = erf_expect(1) + eff_expect(1) + deriv_expect(1) + hcore_expect(1)

END_PROVIDER 
