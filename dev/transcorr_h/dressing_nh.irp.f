
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
! h_dressed = 0.d0
 h_dressed(idress,idress) += delta_u0(idress) 
 
 do i = 1,idress-1
  h_dressed(idress,idress) -= delta_u0(i)/u0(idress) * u0(i)
  h_dressed(idress,i) += delta_u0(i)
  h_dressed(i,idress) += delta_u0(i)
 enddo

 do i = idress+1, N_det
  h_dressed(idress,idress) -= delta_u0(i)/u0(idress) * u0(i)
  h_dressed(idress,i) += delta_u0(i)
  h_dressed(i,idress) += delta_u0(i)
 enddo
end

subroutine get_h_p_delta_psi(psicoef,psidet,ndet,N_st,h_p_delta_psi)
 implicit none
 use bitmasks
 double precision, intent(in)   :: psicoef(ndet,N_st)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet,N_st
 double precision, intent(out)  :: h_p_delta_psi(ndet,N_st) 


 double precision, allocatable :: delta_u0(:,:),hpsi(:,:)
 allocate(delta_u0(ndet,N_st),hpsi(ndet,N_st))
 integer :: idress,i,j,istate
 idress = 1
 call get_delta_and_h_no_store(psidet,psicoef,ndet,N_st,delta_u0,hpsi)
 h_p_delta_psi = hpsi
! h_p_delta_psi = 0.d0
 do istate = 1, N_st
  delta_u0(:,istate) *= 1.d0/psicoef(idress,istate)
 
  h_p_delta_psi(idress,istate) += psicoef(idress,istate) * delta_u0(idress,istate)
  do i = 1, idress-1
   h_p_delta_psi(idress,istate) -= delta_u0(i,istate) * psicoef(i,istate)
   h_p_delta_psi(i,istate) +=  delta_u0(i,istate) * psicoef(idress,istate)
   h_p_delta_psi(idress,istate) +=  delta_u0(i,istate) * psicoef(i,istate)
  enddo

  do i = idress+1, ndet
   h_p_delta_psi(idress,istate) -= delta_u0(i,istate) * psicoef(i,istate)
   h_p_delta_psi(i,istate) +=  delta_u0(i,istate) * psicoef(idress,istate)
   h_p_delta_psi(idress,istate) +=  delta_u0(i,istate) * psicoef(i,istate)
  enddo
 
 enddo

end

subroutine get_delta_and_h_no_store(psidet,psicoef,ndet,Nst,delta,hpsi)
 implicit none
 use bitmasks
 double precision, intent(in)   :: psicoef(ndet,Nst)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet,Nst
 double precision, intent(out)  :: delta(ndet,Nst), hpsi(ndet,Nst)
 double precision :: hij,htilde_ij,delta_mat,hmono,heff,hderiv,hthree
 integer :: i,j,istate
 PROVIDE scalar_mu_r_pot_physicist_mo deriv_mu_r_pot_physicist_mo
 PROVIDE three_body_3_index three_body_3_index_exch_12 three_body_3_index_exch_13 three_body_3_index_exch_23
 PROVIDE three_body_5_index three_body_5_index_exch_13 three_body_5_index_exch_32
 PROVIDE three_body_4_index three_body_4_index_exch_12 three_body_4_index_exch_12_part

 delta = 0.d0
 hpsi = 0.d0
 i=1
 j=1
 call htilde_mat(psi_det(1,1,i),psi_det(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
 call  i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
 !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(dynamic,8) &
 !$OMP PRIVATE(i,j,delta_mat,hmono,heff,hderiv,hthree,htilde_ij,hij)
  do i = 1, N_det
   do j = 1, N_det
    call htilde_mat(psi_det(1,1,i),psi_det(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    call i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
    delta_mat = htilde_ij - hij 
    do istate = 1, Nst
     delta(i,istate) = delta(i,istate) + psicoef(j,istate) * delta_mat
     hpsi(i,istate)  = hpsi(i,istate)  + psicoef(j,istate) * hij
    enddo
   enddo
  enddo
 !$OMP END PARALLEL DO

end
