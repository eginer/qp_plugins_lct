program print
 read_wf = .True.
 touch read_wf
 call provide_all_stuffs
end
subroutine provide_all_stuffs
 implicit none
 provide ref_fobo_hamiltonian_matrix dressing_ref_fobo_hamiltonian
 integer :: i,j,istate
 double precision, allocatable :: psi_restart_ref_fobo_normalized(:),psi_ref_fobo_zeroth_order(:)
 double precision, allocatable :: eigvalues(:),eigvectors(:,:)
 double precision, allocatable :: H_dressed(:,:)
 double precision, allocatable :: H_print(:,:)
 double precision :: accu_norm
 integer, allocatable :: array_print(:)
 allocate (H_dressed(N_det_ref_fobo,N_det_ref_fobo))
 allocate (H_print(N_det_ref_fobo,N_det_ref_fobo))
 allocate (psi_restart_ref_fobo_normalized(N_det_ref_fobo))
 allocate (psi_ref_fobo_zeroth_order(N_det_ref_fobo))
 print*,'# nuclear_repulsion = ',nuclear_repulsion 
 allocate (eigvalues(N_det_ref_fobo))
 allocate (array_print(N_det_ref_fobo))
 allocate (eigvectors(N_det_ref_fobo,N_det_ref_fobo))


 do i = 1, N_det_ref_fobo
  array_print(i) = i
 enddo
 
 do istate= 1, N_states 
  do i = 1, N_det_ref_fobo 
   do j = 1, N_det_ref_fobo
    H_print(i,j) = ref_fobo_hamiltonian_matrix(j,i)
   enddo
  enddo
  do i = 1, N_det_ref_fobo
   H_print(i,i) -= ref_fobo_hamiltonian_matrix(1,1)
  enddo
  print*,'ref_fobo Hamiltonian matrix emelent = ',ref_fobo_hamiltonian_matrix(1,1)
  print*,'ISTATE = ',istate
  accu_norm = 0.d0
  do i = 1, N_det_ref_fobo 
   accu_norm += psi_ref_fobo_coef(i,1) * psi_ref_fobo_coef(i,1)
  enddo
  print*,'accu_norm = ',accu_norm
  accu_norm = 1.d0/dsqrt(accu_norm)
  do i = 1, N_det_ref_fobo 
   psi_restart_ref_fobo_normalized(i) = psi_ref_fobo_coef(i,istate)* accu_norm
  enddo
  print*,'Composition of the wave function'
  do i = 1, N_det_ref_fobo
   if(idx_ref_fobo_cas_lmct_mlct(i,1)==0)then ! CAS 
    print*, i,'CAS'
    call debug_det(psi_ref_fobo(1,1,i),N_int)
   else if (idx_ref_fobo_cas_lmct_mlct(i,1) == +1)then ! LMCT 
    print*, i,'LMCT',idx_ref_fobo_cas_lmct_mlct(i,2)
    call debug_det(psi_ref_fobo(1,1,i),N_int)
    
   else if (idx_ref_fobo_cas_lmct_mlct(i,1) == -1)then ! MLCT 
    print*, i,'MLCT',idx_ref_fobo_cas_lmct_mlct(i,2)
    call debug_det(psi_ref_fobo(1,1,i),N_int)
 
   endif
  enddo

  print*, ''
  print*,'Summary of the composition of the wave function '
  do i = 1, N_det_ref_fobo
   if(idx_ref_fobo_cas_lmct_mlct(i,1)==0)then ! CAS 
    print*, i,'CAS '
   else if (idx_ref_fobo_cas_lmct_mlct(i,1) == +1)then ! LMCT 
    print*, i,'LMCT',idx_ref_fobo_cas_lmct_mlct(i,2)
   else if (idx_ref_fobo_cas_lmct_mlct(i,1) == -1)then ! MLCT 
    print*, i,'MLCT',idx_ref_fobo_cas_lmct_mlct(i,2)
   endif
  enddo

  print*,'-------------------'
  print*,'-------------------'
  print*,'CAS MATRIX         '
  print*,''
  write(*,'(A4,3X,100(I8   ,4X))')'    ',array_print(:)
  do i = 1, N_det_ref_fobo
   write(*,'(I4,3X,100(F8.5 ,4X))')i, H_print(i,:)
  enddo

  print*,''
  print*,'-------------------'
  print*,'-------------------'
  print*,'CAS MATRIX DRESSING'
  print*,''
  
  write(*,'(A4,3X,100(I8   ,4X))')'    ',array_print(:)
  do i = 1, N_det_ref_fobo
   write(*,'(I4,3X,100(F8.5 ,4X))')i, dressing_ref_fobo_hamiltonian(i,:,istate)
  enddo
  print*,''
  print*,'-------------------'
  print*,'-------------------'
  call lapack_diagd(eigvalues,eigvectors,ref_fobo_hamiltonian_matrix,n_det_ref_fobo,n_det_ref_fobo)
  do i = 1, N_det_ref_fobo
   psi_ref_fobo_zeroth_order(i) = eigvectors(i,istate)
  enddo
  
  do i = 1, 3
   print*,'eigvalues',i,eigvalues(i) +  nuclear_repulsion
  enddo
  do i = 2, 3
   print*,'DE = ',eigvalues(1) - eigvalues(i)
  enddo



  do i = 1, N_det_ref_fobo
   do j = 1, N_det_ref_fobo
    H_dressed(j,i) = ref_fobo_hamiltonian_matrix(j,i) + dressing_ref_fobo_hamiltonian(j,i,istate)
    H_print(i,j) += dressing_ref_fobo_hamiltonian(j,i,istate)
   enddo
  enddo
  print*,''
  print*,'-------------------'
  print*,'-------------------'
  print*,'TOTAL DRESSED H MATRIX '
  print*,''
  write(*,'(A4,3X,100(I8   ,4X))')'    ',array_print(:)
  do i = 1, N_det_ref_fobo
   write(*,'(I4,3X,100(F8.5 ,4X))') i,H_print(i,:)
  enddo
  print*,''

  do i = 1, N_det_ref_fobo
   H_print(i,i) -= dressing_ref_fobo_hamiltonian(1,1,istate)
  enddo

  print*,''
  print*,'-------------------'
  print*,'-------------------'
  print*,'TOTAL DRESSED H MATRIX (shifted by the first element)'
  print*,''
  write(*,'(A4,3X,100(I8   ,4X))')'    ',array_print(:)
  do i = 1, N_det_ref_fobo
   write(*,'(I4,3X,100(F8.5 ,4X))') i,H_print(i,:)
  enddo
  print*,''

  print*,'-------------------'
  print*,'-------------------'
  print*,'TOTAL DRESSED H MATRIX in cm-1 (shifted by the first element)'
  print*,''
  write(*,'(A4,3X,100(I8   ,4X))')'    ',array_print(:)
  do i = 1, N_det_ref_fobo
   write(*,'(I4,3X,100(F9.1 ,4X))') i,H_print(i,:) * 219474.63d0
  enddo
  print*,''
  print*,''
  print*,''
  print*,''
 
 
  print*, ''
  print*,'Ref determinant            = ',  nuclear_repulsion + ref_fobo_hamiltonian_matrix(1,1)
  print*,'Ref determinant + dressed  = ',H_print(1,1) +  nuclear_repulsion + ref_fobo_hamiltonian_matrix(1,1)
  print*,'Zeroth-order energy        = ',eigvalues(1) +  nuclear_repulsion 
  print*,'Stabilization zeroth-order = ', eigvalues(1) - ref_fobo_hamiltonian_matrix(1,1)
  print*,'Dressed matrix eigenvalue  = ',energies_ref_fobo_dressed(istate) + nuclear_repulsion
! print*,'Variational energy         = ',psi_energy(istate) + nuclear_repulsion
  print*,'Stabilization second-order = ', energies_ref_fobo_dressed(1) - ref_fobo_hamiltonian_matrix(1,1)
  print*, ''
  print*, 'Composition of the WF '
  print*, ' Type      Number   Exact         Dressed       Naked'
  do i = 1, N_det_ref_fobo
   if(idx_ref_fobo_cas_lmct_mlct(i,1)==0)then ! CAS 
    write(*,'(A5,X,A3,3X,2X,I3,3X,100(F10.7 ,4X))') ' CAS ',' 0 ',i,psi_ref_fobo_coef(i,istate)/psi_ref_fobo_coef(1,istate), psi_ref_fobo_coef_dressed(i,istate)/psi_ref_fobo_coef_dressed(1,istate),psi_ref_fobo_zeroth_order(i)/psi_ref_fobo_zeroth_order(1)
   else if (idx_ref_fobo_cas_lmct_mlct(i,1) == +1)then ! LMCT 
    write(*,'(A5,X,I3,3X,2X,I3,3X,100(F10.7 ,4X))') ' LMCT',idx_ref_fobo_cas_lmct_mlct(i,2),i,psi_ref_fobo_coef(i,istate)/psi_ref_fobo_coef(1,istate), psi_ref_fobo_coef_dressed(i,istate)/psi_ref_fobo_coef_dressed(1,istate),psi_ref_fobo_zeroth_order(i)/psi_ref_fobo_zeroth_order(1)
   else if (idx_ref_fobo_cas_lmct_mlct(i,1) == -1)then ! MLCT 
    write(*,'(A5,X,I3,3X,2X,I3,3X,100(F10.7 ,4X))') ' MLCT',idx_ref_fobo_cas_lmct_mlct(i,2),i,psi_ref_fobo_coef(i,istate)/psi_ref_fobo_coef(1,istate), psi_ref_fobo_coef_dressed(i,istate)/psi_ref_fobo_coef_dressed(1,istate),psi_ref_fobo_zeroth_order(i)/psi_ref_fobo_zeroth_order(1)
   endif
  enddo
 enddo
 
 deallocate (H_dressed)
 deallocate (H_print)
 deallocate (psi_restart_ref_fobo_normalized)
 deallocate (psi_ref_fobo_zeroth_order)

 deallocate (eigvalues,array_print)
 deallocate (eigvectors)

end
