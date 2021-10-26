subroutine set_dress_vec_s(u_in,dets_in,sze,N_st,dress_vec)
 implicit none
  use bitmasks
 double precision, intent(in)   :: u_in(sze)
 integer(bit_kind), intent(in)  :: dets_in(N_int,2,sze)
 integer, intent(in)            :: sze,N_st
 double precision, intent(out)  :: dress_vec(sze,N_st)
 integer :: i,ii,k,j, l
 double precision :: f, tmp
 double precision, allocatable :: delta(:)

 allocate(delta(sze))
 call get_delta_no_store_no_provide(dets_in,u_in,sze,delta)

 dress_vec(:,:) = 0.d0

 l = 1
 do j = 1, n_det
   if (j == l) cycle
   dress_vec(j,1)  = delta(j) 
   dress_vec(l,1) -= u_in(j) * delta(j) / u_in(l)
 enddo
 dress_vec(l,1) += delta(l) 
 dress_vec(l,1) *= 0.5d0
end


subroutine get_delta_no_store_no_provide(psidet,psicoef,ndet,delta)
 implicit none
 use bitmasks
 double precision, intent(in)   :: psicoef(ndet)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet
 double precision, intent(out)  :: delta(ndet) 
 double precision :: hij,htilde_ij,delta_mat,hmono,heff,hderiv,hthree
 integer :: i,j

 delta = 0.d0
 i=1
 j=1
 call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
 call i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
 !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(dynamic,8) &
 !$OMP PRIVATE(i,j,delta_mat,hmono,heff,hderiv,hthree,htilde_ij,hij)
  do i = 1, N_det
   do j = 1, N_det
    call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    call i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
    delta_mat = htilde_ij - hij 
    delta(i) = delta(i) + psicoef(j) * delta_mat
   enddo
  enddo
 !$OMP END PARALLEL DO

end

subroutine htilde_psi_no_store_no_provide(psidet,psicoef,ndet,htilde_psi)
 implicit none
 use bitmasks
 double precision, intent(in)   :: psicoef(ndet)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet
 double precision, intent(out)  :: htilde_psi(ndet) 
 double precision :: hij,htilde_ij,htilde_psi_mat,hmono,heff,hderiv,hthree
 integer :: i,j

 htilde_psi = 0.d0
 i=1
 j=1
 call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
 !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(dynamic,8) &
 !$OMP PRIVATE(i,j,htilde_psi_mat,hmono,heff,hderiv,hthree,htilde_ij)
  do i = 1, N_det
   do j = 1, N_det
    call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    htilde_psi(i) = htilde_psi(i) + psicoef(j) * htilde_ij
   enddo
  enddo
 !$OMP END PARALLEL DO

end

