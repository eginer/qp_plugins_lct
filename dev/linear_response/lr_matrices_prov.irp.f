 BEGIN_PROVIDER[double precision, A_IJ, (N_det-1 , N_det-1, N_states)]
 implicit none
 BEGIN_DOC
 !
 END_DOC 
 integer :: i_state
 integer :: det_I, det_J
 integer :: det_II, det_JJ
 double precision :: contrib_to_AIJ
 double precision, allocatable :: kernel_lr_for_prov(:,:,:),H_matrix_for_prov(:,:)
 
 allocate(kernel_lr_for_prov(N_det,N_det,N_states),H_matrix_for_prov(N_det,N_det))
 kernel_lr_for_prov = kernel_lr
 H_matrix_for_prov = H_matrix_all_dets

 A_IJ = 0.d0

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
    A_IJ(det_II, det_JJ, i_state) = H_matrix_for_prov(det_I,det_J) - contrib_to_AIJ + kernel_lr_for_prov(det_I,det_J,i_state)  
   enddo
  enddo
 enddo

 deallocate(kernel_lr_for_prov,H_matrix_for_prov)
 END_PROVIDER

 BEGIN_PROVIDER[double precision, A_IJ_K_equal_0, (N_det-1, N_det-1, N_states)]
 implicit none
 BEGIN_DOC
 !
 END_DOC 
 integer :: i_state
 integer :: det_I, det_J
 integer :: det_II, det_JJ
 double precision :: contrib_to_AIJ
 double precision, allocatable :: H_matrix_for_prov(:,:)

 allocate(H_matrix_for_prov(N_det,N_det))
 H_matrix_for_prov = H_matrix_all_dets

 A_IJ_K_equal_0 = 0.d0

!print*, 'H_matrix_all_dets'
!do det_I = 1, N_det
! write(*,'(100(f10.5,x))'), H_matrix_all_dets(1:N_det,det_I)
!enddo

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
    A_IJ_K_equal_0(det_II, det_JJ, i_state) = H_matrix_for_prov(det_I,det_J) - contrib_to_AIJ 
   enddo
  enddo
 enddo

 deallocate(H_matrix_for_prov)

 END_PROVIDER

 BEGIN_PROVIDER[double precision, B_IJ, (N_det-1, N_det-1, N_states)]
 implicit none
 BEGIN_DOC
 !
 END_DOC 
 integer :: i_state
 integer :: det_I, det_J
 integer :: det_II, det_JJ
 double precision, allocatable :: kernel_lr_for_prov(:,:,:)
 
 allocate(kernel_lr_for_prov(N_det,N_det,N_states))
 kernel_lr_for_prov = kernel_lr

 B_IJ = 0.d0

!print*, 'H_matrix_all_dets'
!do det_I = 1, N_det
! write(*,'(100(f10.5,x))'), H_matrix_all_dets(1:N_det,det_I)
!enddo

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
    B_IJ(det_II, det_JJ, i_state) = kernel_lr_for_prov(det_I,det_J,i_state)  
   enddo
  enddo
 enddo

 deallocate(kernel_lr_for_prov)
 END_PROVIDER

 BEGIN_PROVIDER[double precision, S_IJ, (N_det-1, N_det-1, N_states)]
 implicit none
 BEGIN_DOC
 !
 END_DOC 
 integer :: i_state
 integer :: det_I, det_J
 integer :: det_II, det_JJ
 double precision :: contrib
 
 S_IJ = 0.d0

!print*, 'H_matrix_all_dets'
!do det_I = 1, N_det
! write(*,'(100(f10.5,x))'), H_matrix_all_dets(1:N_det,det_I)
!enddo

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
    contrib = 0.d0
    if (det_II.eq.det_JJ)then
     contrib = 1.d0
    endif
    S_IJ(det_II, det_JJ, i_state) = contrib - psi_coef(det_I,i_state) * psi_coef(det_J,i_state)
   enddo
  enddo
 enddo

 END_PROVIDER

