program chec_ref_energy
 implicit none
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 call routine 
end

subroutine routine
 implicit none
 double precision :: hmono,herf,heff,hderiv,hthree,htot
 call direct_diag_htilde_mu_mat(ref_bitmask,hmono,herf,heff,hderiv,hthree,htot)
 print*,'hmono       = ',hmono
 print*,'herf        = ',herf
 print*,'heff        = ',heff
 print*,'hderiv      = ',hderiv
 print*,'twobody     = ',herf+heff+hderiv
 print*,'hthree      = ',hthree
 print*,'htot        = ',htot
 print*,'core_energy = ',core_energy

end
