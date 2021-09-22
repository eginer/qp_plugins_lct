

subroutine get_delta_no_store(psidet,psicoef,ndet,delta)
 implicit none
 use bitmasks
 double precision, intent(in)   :: psicoef(ndet)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet
 double precision, intent(out)  :: delta(ndet) 
 double precision :: hij,htilde_ij,delta_mat,hmono,heff,hderiv,hthree
 integer :: i,j
 delta = 0.d0
  do j = 1, N_det
   do i = 1, N_det
    call htilde_mat(psi_det(1,1,i),psi_det(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    call  i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
    delta_mat = htilde_ij - hij 
    delta(i) += psicoef(j) * delta_mat
   enddo
  enddo

end
