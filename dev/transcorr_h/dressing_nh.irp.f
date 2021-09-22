
subroutine get_dressed_matrix(u0,h_dressed,idress)
 use bitmasks
 BEGIN_DOC
! You enter with u0, a good guess to the right eigenvector of the TC Hamiltonian
!
! You get out with a dressed symmetric matrix taking the effect of (Htilde - H) |u0>
 END_DOC
 implicit none
 integer, intent(in) :: idress
 double precision, intent(in) :: u0(N_det)
 double precision, intent(out):: h_dressed(N_det,N_det)
 double precision, allocatable :: delta_u0(:), delta_mat(:,:)
 double precision :: a
 integer :: i
 a = 1.d0
 allocate(delta_u0(N_det),delta_mat(N_det,N_det))
 delta_mat = htilde_matrix_elmt - h_matrix_all_dets ! Delta = Htilde - H
 !!!!!!!!!!!!! Computing the dressing vector 
 delta_u0 = 0.d0
 call h_non_hermite(delta_u0,u0(idress),delta_mat,a,1,N_det)  ! delta_u0 = Delta |u0> 
 delta_u0 *= 1.d0/u0(idress)
 !!!!!!!!!!!!! Computing the dressing matrix 
 h_dressed = h_matrix_all_dets
 h_dressed(idress,idress) += delta_u0(idress) 
 
 do i = 2, N_det
  h_dressed(idress,idress) -= delta_u0(i)/u0(idress) * u0(i)
  h_dressed(idress,i) += delta_u0(i)
  h_dressed(i,idress) += delta_u0(i)
 enddo
end

