subroutine get_dressing_tc_for_dav(u_in,dets_in,sze,N_st,dagger,dress_vec)
 implicit none
 BEGIN_DOC
! you enter with a wave function of determinants dets_in and coefs u_in
!
! you get out with a dressing vector to dress_vec be used in the davidson dressed
 END_DOC
  use bitmasks
 double precision, intent(in)   :: u_in(sze)
 integer(bit_kind), intent(in)  :: dets_in(N_int,2,sze)
 integer, intent(in)            :: sze,N_st
 logical, intent(in)            :: dagger
 double precision, intent(out)  :: dress_vec(sze,N_st)
 integer :: i,ii,k,j, l
 double precision :: f, tmp
 double precision, allocatable :: delta(:)

 allocate(delta(sze))
 if(dagger)then
  call get_delta_tc_dagger_psi(dets_in,u_in,sze,delta)
 else
  call get_delta_tc_psi(dets_in,u_in,sze,delta)
 endif

 dress_vec(:,:) = 0.d0

 l = 1
 do j = 1, sze
   if (j == l) cycle
   dress_vec(j,1)  = delta(j) 
   dress_vec(l,1) -= u_in(j) * delta(j) / u_in(l)
 enddo
 dress_vec(l,1) += delta(l) 
 dress_vec(l,1) *= 0.5d0
end


subroutine get_delta_tc_psi(psidet,psicoef,ndet,delta)
 implicit none
 BEGIN_DOC
! you get in with a wave function psidet,psicoef and you get out with 
!
! |delta> = (Htilde - H) |Psi>
 END_DOC
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
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(delta,ndet,psidet,psicoef,N_int) &
 !$OMP PRIVATE(i,j,delta_mat,hmono,heff,hderiv,hthree,htilde_ij,hij)
  do i = 1, ndet
   do j = 1, ndet
    call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    call i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
    delta_mat = htilde_ij - hij 
    delta(i) = delta(i) + psicoef(j) * delta_mat
   enddo
  enddo
 !$OMP END PARALLEL DO

end

subroutine get_delta_av_tc_psi(psidet,psicoef,ndet,delta)
 implicit none
 BEGIN_DOC
! you get in with a wave function psidet,psicoef and you get out with 
!
! |delta> = (Htilde - H) |Psi>
 END_DOC
 use bitmasks
 double precision, intent(in)   :: psicoef(ndet)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet
 double precision, intent(out)  :: delta(ndet) 
 double precision :: hij,htilde_ij,delta_mat,hmono,heff,hderiv,hthree,htilde_ji
 integer :: i,j

 delta = 0.d0
 i=1
 j=1
 call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
 call i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(delta,ndet,psidet,psicoef,N_int) &
 !$OMP PRIVATE(i,j,delta_mat,hmono,heff,hderiv,hthree,htilde_ij,htilde_ji,hij)
  do i = 1, ndet
   do j = 1, ndet
    call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    call htilde_mu_mat(psidet(1,1,j),psidet(1,1,i),hmono,heff,hderiv,hthree,htilde_ji)
    call i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
!    delta_mat = htilde_ij - hij 
    delta_mat = 0.5d0 * (htilde_ij + htilde_ji)
    delta(i) = delta(i) + psicoef(j) * delta_mat
   enddo
  enddo
 !$OMP END PARALLEL DO

end

subroutine get_e_components_htilde(psidet,psicoef,ndet,hmono_av,heff_av,hderiv_av,hthree_av,htot_av)
 use bitmasks
 implicit none
 double precision, intent(in)   :: psicoef(ndet)
 integer(bit_kind), intent(in)  :: psidet(N_int,2,ndet)
 integer, intent(in)            :: ndet
 double precision, intent(out)  :: hmono_av,heff_av,hderiv_av,hthree_av,htot_av
 double precision :: hij,htot,htilde_psi_mat,hmono,heff,hderiv,hthree,u_dot_v
 double precision, allocatable :: hmono_vec(:),heff_vec(:),hderiv_vec(:),hthree_vec(:)
 integer :: i,j
 allocate(hmono_vec(ndet),heff_vec(ndet),hderiv_vec(ndet),hthree_vec(ndet))

 hmono_av  = 0.d0
 heff_av   = 0.d0
 hderiv_av = 0.d0
 hthree_av = 0.d0

 hmono_vec = 0.d0
 heff_vec = 0.d0
 hderiv_vec = 0.d0
 hthree_vec = 0.d0
 i=1
 j=1
 call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htot)
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(hmono_vec,heff_vec,hderiv_vec,hthree_vec,ndet,psidet,psicoef) &
 !$OMP PRIVATE(i,j,htilde_psi_mat,hmono,heff,hderiv,hthree,htot)
  do i = 1, ndet
   do j = 1, ndet
    call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htot)
    hmono_vec(i)  += psicoef(j) * hmono  
    heff_vec(i)   += psicoef(j) * heff   
    hderiv_vec(i) += psicoef(j) * hderiv 
    hthree_vec(i) += psicoef(j) * hthree 
   enddo
  enddo
 !$OMP END PARALLEL DO
 do i = 1, ndet
  hmono_av  += psicoef(i) * hmono_vec(i)
  heff_av   += psicoef(i) * heff_vec(i)
  hderiv_av += psicoef(i) * hderiv_vec(i)
  hthree_av += psicoef(i) * hthree_vec(i)
 enddo
! hmono_av  = u_dot_v(psicoef,hmono_vec,ndet)
! heff_av   = u_dot_v(psicoef,heff_vec,ndet)
! hderiv_av = u_dot_v(psicoef,hderiv_vec,ndet)
! hthree_av = u_dot_v(psicoef,hthree_vec,ndet)

 htot_av = hmono_av + heff_av + hderiv_av + hthree_av


end

subroutine get_htilde_psi(psidet,psicoef,ndet,htilde_psi)
 BEGIN_DOC
! you get in with a wave function psidet,psicoef and you get out with 
!
! |delta> = Htilde |Psi>
 END_DOC
 use bitmasks
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
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(htilde_psi,ndet,psidet,psicoef) & 
 !$OMP PRIVATE(i,j,htilde_psi_mat,hmono,heff,hderiv,hthree,htilde_ij)
  do i = 1, ndet
   do j = 1, ndet
    call htilde_mu_mat(psidet(1,1,i),psidet(1,1,j),hmono,heff,hderiv,hthree,htilde_ij)
    htilde_psi(i) = htilde_psi(i) + psicoef(j) * htilde_ij
   enddo
  enddo
 !$OMP END PARALLEL DO

end

subroutine get_delta_tc_dagger_psi(psidet,psicoef,ndet,delta)
 implicit none
 BEGIN_DOC
! you get in with a wave function psidet,psicoef and you get out with 
!
! |delta> = (Htilde - H)^{DAGGER} |Psi>
!
 END_DOC
 use bitmasks
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
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(delta,ndet,psidet,psicoef,N_int) & 
 !$OMP PRIVATE(i,j,delta_mat,hmono,heff,hderiv,hthree,htilde_ij,hij)
  do i = 1, ndet
   do j = 1, ndet
    ! just changed i<-->j with respect to get_delta_tc_psi 
    call htilde_mu_mat(psidet(1,1,j),psidet(1,1,i),hmono,heff,hderiv,hthree,htilde_ij)
    call i_H_j(psidet(1,1,i),psidet(1,1,j),N_int,hij)
    delta_mat = htilde_ij - hij 
    delta(i) = delta(i) + psicoef(j) * delta_mat
   enddo
  enddo
 !$OMP END PARALLEL DO

end

subroutine get_htilde_dagger_psi(psidet,psicoef,ndet,htilde_psi)
 implicit none
 BEGIN_DOC
! you get in with a wave function psidet,psicoef and you get out with 
!
! |delta> = (Htilde)^DAGGER |Psi>
 END_DOC
 use bitmasks
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
 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
 !$OMP SHARED(htilde_psi,ndet,psidet,psicoef)&
 !$OMP PRIVATE(i,j,htilde_psi_mat,hmono,heff,hderiv,hthree,htilde_ij)
  do i = 1, ndet
   do j = 1, ndet
    ! just changed i<-->j with respect to htilde_psi
    call htilde_mu_mat(psidet(1,1,j),psidet(1,1,i),hmono,heff,hderiv,hthree,htilde_ij)
    htilde_psi(i) = htilde_psi(i) + psicoef(j) * htilde_ij
   enddo
  enddo
 !$OMP END PARALLEL DO

end

