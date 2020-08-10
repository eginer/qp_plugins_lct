program pouet
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 use bitmasks
 integer(bit_kind) :: key_i(N_int,2)
 integer(bit_kind) :: key_j(N_int,2)
 double precision :: hij,s2,hmono,herf,heff,hderiv
 key_i(:,:) = psi_det(:,:,1)
 call i_H_j_s2(key_i,key_i,N_int,hij,s2)
 call diag_htilde_mat(key_i,hmono,herf,heff,hderiv)
 print*,'hij          = ',hij
 print*,'hmono + herf = ',hmono + herf
 print*,'heff         = ',heff 
 print*,'hderiv       = ',hderiv
end
