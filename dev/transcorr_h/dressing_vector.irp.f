! BEGIN_PROVIDER [ double precision, dressing_column_h, (N_det,N_states) ]
!&BEGIN_PROVIDER [ double precision, dressing_column_s, (N_det,N_states) ]
! implicit none
! BEGIN_DOC
! ! \Delta_{state-specific}. \Psi
! ! Diagonal element is divided by 2 because Delta = D + D^t
! END_DOC
!
! integer :: i,ii,k,j, l
! double precision :: f, tmp
! double precision, allocatable :: delta(:)
!
! allocate(delta(N_det))
! call get_delta_no_store(psi_det,psi_coef,N_det,delta)
!
! dressing_column_h(:,:) = 0.d0
! dressing_column_s(:,:) = 0.d0
!
! l = 1
! do j = 1, n_det
!   if (j == l) cycle
!   dressing_column_h(j,1)  = delta(j) 
!   dressing_column_h(l,1) -= psi_coef(j,1) * delta(j) / psi_coef(l,1)
! enddo
! dressing_column_h(l,1) += delta(l) 
! dressing_column_h(l,1) *= 0.5d0
!
!END_PROVIDER


subroutine delta_to_dressing_vector(psicoef,delta,ndet,dressing,idress)
 implicit none 
 BEGIN_DOC
! you start with a delta(i) = <I|Delta |psicoef> 
!
! and you get out with a dressing(i) which can be used in a dressed davidson routine 
 END_DOC
 double precision, intent(in)   :: psicoef(ndet),delta(ndet)
 integer, intent(in)            :: ndet,idress
 double precision, intent(out)  :: dressing(ndet) 
! call get_delta_no_store(psidet,psicoef,ndet,delta)
 integer :: l,j
 l = idress
 dressing = 0.d0
 do j = 1, ndet
   if (j == l) cycle
   dressing(j)  = delta(j) 
   dressing(l) -= psicoef(j) * delta(j) / psicoef(l)
 enddo
 dressing(l) += delta(l)
 dressing(l) *= 0.5d0

end

