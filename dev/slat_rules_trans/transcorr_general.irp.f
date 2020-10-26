program transcorr_h
 implicit none
 read_wf = .True.
 touch read_wf
 call provide_all
 call print_energy
 call print_e_comp_transcorr
 call print_eigv
 call print_pert
! call plot_on_top_left_right
! call print_psi_exc_psi_trans
 call write_left_right
end

subroutine provide_all
 use bitmasks
 integer(bit_kind) :: key_i(N_int,2), key_j(N_int,2)
 integer :: i,j,degree
 double precision :: hij,s2,hmono,herf,heff,hderiv,htot
 double precision :: accu
 accu = 0.d0
 key_i(:,:) = psi_det(:,:,1)
 call diag_htilde_mat(key_i,hmono,herf,heff,hderiv,htot)
 provide eigval_trans
end

