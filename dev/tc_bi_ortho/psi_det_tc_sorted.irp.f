use bitmasks

BEGIN_PROVIDER [ double precision, psi_average_norm_contrib_tc, (psi_det_size) ]                                                                
  implicit none
  BEGIN_DOC
  ! Contribution of determinants to the state-averaged density.
  END_DOC
  integer                        :: i,j,k
  double precision               :: f
    
  psi_average_norm_contrib_tc(:) = 0.d0
  do k=1,N_states
    do i=1,N_det
      psi_average_norm_contrib_tc(i) = psi_average_norm_contrib_tc(i) +    &
          dabs(psi_l_coef_bi_ortho(i,k)*psi_l_coef_bi_ortho(i,k))*state_average_weight(k)
    enddo
  enddo
  f = 1.d0/sum(psi_average_norm_contrib_tc(1:N_det))
  do i=1,N_det
    psi_average_norm_contrib_tc(i) = psi_average_norm_contrib_tc(i)*f
  enddo
END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_tc, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_tc, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ double precision, psi_average_norm_contrib_sorted_tc, (psi_det_size) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_tc_order, (psi_det_size) ]
   implicit none
   BEGIN_DOC
   ! Wave function sorted by determinants contribution to the norm (state-averaged)
   !
   ! psi_det_sorted_tc_order(i) -> k : index in psi_det
   END_DOC
   integer                        :: i,j,k
   integer, allocatable           :: iorder(:)
   allocate ( iorder(N_det) )
   do i=1,N_det
     psi_average_norm_contrib_sorted_tc(i) = -psi_average_norm_contrib_tc(i)
     iorder(i) = i
   enddo
   call dsort(psi_average_norm_contrib_sorted_tc,iorder,N_det)
   do i=1,N_det
     do j=1,N_int
       psi_det_sorted_tc(j,1,i) = psi_det(j,1,iorder(i))
       psi_det_sorted_tc(j,2,i) = psi_det(j,2,iorder(i))
     enddo
    psi_average_norm_contrib_sorted_tc(i) = -psi_average_norm_contrib_sorted_tc(i)
    psi_det_sorted_tc_order(iorder(i)) = i
   enddo
   double precision :: accu
   do k=1,N_states
    accu = 0.d0
    do i=1,N_det
     psi_coef_sorted_tc(i,k) = dsqrt(dabs(psi_l_coef_bi_ortho(iorder(i),k)*psi_l_coef_bi_ortho(iorder(i),k)))
     accu += psi_coef_sorted_tc(i,k)**2
    enddo
    accu = 1.d0/dsqrt(accu)
    do i=1,N_det
     psi_coef_sorted_tc(i,k) *= accu
    enddo
   enddo

   psi_det_sorted_tc(:,:,N_det+1:psi_det_size) = 0_bit_kind
   psi_coef_sorted_tc(N_det+1:psi_det_size,:) = 0.d0
   psi_average_norm_contrib_sorted_tc(N_det+1:psi_det_size) = 0.d0
   psi_det_sorted_tc_order(N_det+1:psi_det_size) = 0

   deallocate(iorder)

END_PROVIDER

 BEGIN_PROVIDER [double precision, reigvec_tc_bi_orth_sorted, (psi_det_size, N_states)]
&BEGIN_PROVIDER [double precision, leigvec_tc_bi_orth_sorted, (psi_det_size, N_states)]

   implicit none
   integer                       :: i, j, k
   reigvec_tc_bi_orth_sorted = 0.d0
   leigvec_tc_bi_orth_sorted = 0.d0
   do i = 1, N_det
    reigvec_tc_bi_orth_sorted(i,1) = reigvec_tc_bi_orth(psi_det_sorted_tc_order(i),1)
    leigvec_tc_bi_orth_sorted(i,1) = leigvec_tc_bi_orth(psi_det_sorted_tc_order(i),1)
   enddo

END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_tc_bit, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_tc_bit, (psi_det_size,N_states) ]
   implicit none
   BEGIN_DOC
   ! Determinants on which we apply $\langle i|H|psi \rangle$ for perturbation.
   ! They are sorted by determinants interpreted as integers. Useful
   ! to accelerate the search of a random determinant in the wave
   ! function.
   END_DOC

   call sort_dets_by_det_search_key(N_det, psi_det, psi_coef, size(psi_coef,1),       &
       psi_det_sorted_tc_bit, psi_coef_sorted_tc_bit, N_states)

END_PROVIDER

