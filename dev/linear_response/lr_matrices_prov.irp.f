 BEGIN_PROVIDER[double precision, A_IJ, (N_det-1 , N_det-1, N_states)]
 &BEGIN_PROVIDER[double precision, A_IJ_K_equal_0, (N_det-1, N_det-1, N_states)]
 &BEGIN_PROVIDER[double precision, B_IJ, (N_det-1, N_det-1, N_states)]
 &BEGIN_PROVIDER[double precision, B_IJ_K_equal_0, (N_det-1, N_det-1, N_states)]
 &BEGIN_PROVIDER[double precision, S_IJ, (N_det-1, N_det-1, N_states)]
 implicit none
 BEGIN_DOC
 !
 END_DOC 
 integer :: i_state
 integer :: det_I, det_J
 integer :: det_II, det_JJ
 double precision :: contrib_to_AIJ, contrib
 
 print*, "psi_energy(1)=",psi_energy(1)

 A_IJ = 0.d0
 A_IJ_K_equal_0 = 0.d0
 B_IJ = 0.d0
 B_IJ_K_equal_0 = 0.d0
 S_IJ = 0.d0

 det_II = 0
! print*, "hf_index=",hf_index
 do i_state=1, N_states
  do det_I = 1, N_det
   if(det_I == hf_index) cycle
   det_II+=1
   det_JJ=0
   do det_J = 1, N_det
    if (det_J == hf_index) cycle
    det_JJ+=1
    contrib_to_AIJ = 0.d0
    if (det_II.eq.det_JJ)then
     contrib_to_AIJ = psi_energy(1)  
    endif
    A_IJ(det_II, det_JJ, i_state) = H_matrix_all_dets(det_I,det_J) - contrib_to_AIJ + kernel_lr(det_I,det_J,i_state)  
    A_IJ_K_equal_0(det_II, det_JJ, i_state) = H_matrix_all_dets(det_I,det_J) - contrib_to_AIJ 
    B_IJ(det_II, det_JJ, i_state) = kernel_lr(det_I,det_J,i_state)  
    contrib = 0.d0
    if (det_II.eq.det_JJ)then
     contrib = 1.d0
    endif
    S_IJ(det_II, det_JJ, i_state) = contrib - psi_coef(det_I,i_state) * psi_coef(det_J,i_state)
   enddo
  enddo
 enddo

 END_PROVIDER

!BEGIN_PROVIDER[double precision, A_IJ, (N_det-1, N_det-1, N_states)]
!&BEGIN_PROVIDER[double precision, A_IJ_K_equal_0, (N_det-1, N_det-1, N_states)]
!&BEGIN_PROVIDER[double precision, B_IJ, (N_det-1, N_det-1, N_states)]
!&BEGIN_PROVIDER[double precision, B_IJ_K_equal_0, (N_det-1, N_det-1, N_states)]
!&BEGIN_PROVIDER[double precision, S_IJ, (N_det-1, N_det-1, N_states)]
!implicit none
!BEGIN_DOC
!!
!END_DOC 
!integer :: i_state
!integer :: det_I, det_J
!double precision :: contrib_to_AIJ, contrib
!
!A_IJ = 0.d0
!A_IJ_K_equal_0 = 0.d0
!B_IJ = 0.d0
!B_IJ_K_equal_0 = 0.d0
!S_IJ = 0.d0
!do i_state=1, N_states
! do det_I = 1, N_det-1
!  do det_J = 1, N_det-1
!   contrib_to_AIJ = 0.d0
!   if (det_I.eq.det_J)then
!    contrib_to_AIJ = psi_energy(1)  
!   endif
!   A_IJ(det_I, det_J, i_state) = H_matrix_all_dets(det_I+1,det_J+1) - contrib_to_AIJ + kernel_lr(det_I+1,det_J+1,i_state)  
!   A_IJ_K_equal_0(det_I, det_J, i_state) = H_matrix_all_dets(det_I+1,det_J+1) - contrib_to_AIJ 
!   B_IJ(det_I, det_J, i_state) = kernel_lr(det_I+1,det_J+1,i_state)  
!   contrib = 0.d0
!   if (det_I.eq.det_J)then
!    contrib = 1.d0
!   endif
!   S_IJ(det_I, det_J, i_state) = contrib - psi_coef(det_I+1,i_state) * psi_coef(det_J+1,i_state)
!  enddo
! enddo
!enddo

!END_PROVIDER

