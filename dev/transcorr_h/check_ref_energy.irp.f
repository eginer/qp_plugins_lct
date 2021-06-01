program chec_ref_energy
 implicit none
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 read_wf = .True.
 touch read_wf
! call routine 
 call test_diag
end

subroutine routine
 implicit none
 double precision :: hmono,heff,hderiv,hthree,htot
 call direct_diag_htilde_mu_mat(ref_bitmask,hmono,heff,hderiv,hthree,htot)
 print*,'hmono       = ',hmono
 print*,'heff        = ',heff
 print*,'hderiv      = ',hderiv
 print*,'twobody     = ',heff+hderiv
 print*,'hthree      = ',hthree
 print*,'htot        = ',htot
 print*,'core_energy = ',core_energy

end

subroutine test_diag
 implicit none
 integer :: i
 double precision :: hmono,heff,hderiv,hthree,htot
 double precision :: hbis_mono,hbis_eff,hbis_deriv,hbis_three,hbis_tot
 do i = 1, N_det
  call diag_htilde_mu_mat_3_index(psi_det(1,1,i),hmono,heff,hderiv,hthree,htot)
  call direct_diag_htilde_mu_mat(psi_det(1,1,i),hbis_mono,hbis_eff,hbis_deriv,hbis_three,hbis_tot)
  print*,'htot, hbis_tot',htot, hbis_tot
  if(dabs(hthree - hbis_three).gt.1.d-12)then
   print*,'Error !'
   print*,'Good, test, error '
   print*,hbis_three, hthree, dabs(hthree - hbis_three)
  endif

 enddo


end
