program pouet
 implicit none
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 use bitmasks
 integer(bit_kind) :: key_i(N_int,2), key_j(N_int,2)
 integer :: i,j,degree
 double precision :: hij,s2,hmono,herf,heff,hderiv,htot
 double precision :: accu
 accu = 0.d0
 key_i(:,:) = psi_det(:,:,1)
 call diag_htilde_mat(key_i,hmono,herf,heff,hderiv,htot)
 do i = 1, N_det
  key_i(:,:) = psi_det(:,:,i)
  do j = 1, N_det
   key_j(:,:) = psi_det(:,:,j)
   call i_H_j_s2(key_i,key_j,N_int,hij,s2)
   call htilde_mat(key_j,key_i,hmono,herf,heff,hderiv,htot)
   print*,'hij,htot',hij,htot
   accu += dabs(hij - htot )
  enddo
 enddo
 print*,'accu = ',accu
end
