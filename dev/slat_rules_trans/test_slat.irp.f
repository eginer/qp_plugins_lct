program pouet
 implicit none
 read_wf = .True.
 touch read_wf
! call routine
! call test_pert
 call test_eigv
 integer :: ith
 ith = 1
! call test_left_right_eigenvalues(ith)
! call test_overlap_matrix
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
 print*,''
 print*,'H      eigenvalue'
 i = 1
 print*,'E0       = ',CI_energy(i)
 print*,''
 print*,'Htilde eigenvalue'
 i = 1
  print*,'E0_tilde = ',eigval_trans(i) + nuclear_repulsion
end

subroutine test_eigv
 implicit none
 integer :: i,j
 do i = 1, N_det
  print*,reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1)),leigvec_trans(i,1)/dsqrt(leigvec_trans_norm(1)),psi_coef(i,1)
 enddo
 double precision :: accu1,accu2,e
 accu1 = 0.d0
 accu2 = 0.d0
 do i = 1, N_det
  accu2 += reigvec_trans(i,1) * reigvec_trans(i,1)
  do j = 1, N_det
   accu1 += htilde_matrix_elmt(j,i) * reigvec_trans(i,1) * reigvec_trans(j,1)
  enddo
 enddo
 accu2 = dsqrt(accu2)
 e = accu1/accu2
 print*,'right eigenvector'
 print*,'e                = ',e
 accu1 = 0.d0
 accu2 = 0.d0
 do i = 1, N_det
  accu2 += leigvec_trans(i,1) * leigvec_trans(i,1)
  do j = 1, N_det
   accu1 += htilde_matrix_elmt(j,i) * leigvec_trans(i,1) * leigvec_trans(j,1)
  enddo
 enddo
 accu2 = dsqrt(accu2)
 e = accu1/accu2
 print*,'left eigenvector'
 print*,'e                = ',e
 print*,'eigval_trans(1)  = ',eigval_trans(1)
 do i = 1, N_det
  psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
 enddo
 touch psi_coef
 call save_wavefunction
end

subroutine test_pert
 implicit none
 integer :: i,j
 double precision :: accu,hmono,herf,heff,hderiv,htot
 double precision :: accu_mono,accu_double
 double precision :: pert_mono,pert_double
 double precision :: h00,hii,htotbis,phase
 integer          :: degree,exc(0:2,2,2)
 integer          :: h1, p1, h2, p2, s1, s2
 accu = 0.d0
 accu_mono = 0.d0
 accu_double = 0.d0
 pert_mono = 0.d0
 pert_double = 0.d0
 call htilde_mat(psi_det(1,1,1),psi_det(1,1,1),hmono,herf,heff,hderiv,h00)
 do i = 1, N_det
  ! <i|H|0>
  call htilde_mat(psi_det(1,1,i),psi_det(1,1,1),hmono,herf,heff,hderiv,htot)

  ! <i|H|i>
  call htilde_mat(psi_det(1,1,i),psi_det(1,1,i),hmono,herf,heff,hderiv,hii)

  ! <0|H|i>
  call htilde_mat(psi_det(1,1,1),psi_det(1,1,i),hmono,herf,heff,hderiv,htotbis)

  call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,1),degree,N_int)
  if(degree==1)then
   call get_single_excitation(psi_det(1,1,i),psi_det(1,1,1),exc,phase,N_int)
   call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
   print*,'h1,p1',h1,p1
   accu_mono += htotbis * reigvec_trans(i,1)/reigvec_trans(1,1)
   pert_mono += htot*htotbis/(h00 - hii)
  else if (degree==2)then
   call debug_det(psi_det(1,1,i),N_int)
   call debug_det(psi_det(1,1,1),N_int)
   call get_double_excitation(psi_det(1,1,i),psi_det(1,1,1),exc,phase,N_int)
   call decode_exc(exc,2,h1,p1,h2,p2,s1,s2)
   print*,'h1,p1,h2,p2',h1,p1,h2,p2
   accu_double += htotbis * reigvec_trans(i,1)/reigvec_trans(1,1)
   pert_double += htot*htotbis/(h00 - hii)
  endif
  accu += htotbis * reigvec_trans(i,1)/reigvec_trans(1,1)
 enddo
 print*,''
 print*,''
 print*,'accu         = ',accu
 print*,'eigval_trans = ',eigval_trans(1)
 print*,'<d0|H|d0>    = ',h00 
 print*,'E corr       = ',eigval_trans(1) - h00
 print*,'accu_mono    = ',accu_mono
 print*,'accu_double  = ',accu_double
 print*,'pert_mono    = ',pert_mono
 print*,'pert_double  = ',pert_double


end
