program pouet
 implicit none
! call test_general_lapack
 call test_dav
end

subroutine test_dav
 implicit none
 integer :: sze,nstates
 double precision, allocatable :: reigv(:,:),leigv(:,:),eigval(:)
 external hcalc_r_tmp
 external hcalc_l_tmp
 sze = 2
 nstates = 1
 allocate(reigv(sze,nstates),leigv(sze,nstates),eigval(nstates))
 reigv = 0.d0
 leigv = 0.d0
 reigv(1,1) = 1.d0
 leigv(1,1) = 1.d0
 eigval(1) = h_non_hermit(1,1)
 call non_hermit_dav(sze,nstates,leigv,reigv,eigval,hcalc_l_tmp,hcalc_r_tmp)
 print*,'eigval_ht = ',eigval_ht(1)
 print*,'eigval    = ',eigval(1)
 integer :: i
 do i = 1, sze
  write(*,'(I4,X,2(F16.10,X))')i,reigv(i,1),reigvec_ht(i,1)
 enddo

end

subroutine test_general_lapack
 implicit none
 double precision, allocatable :: B(:,:),A(:,:)
 double precision, allocatable :: reigvec(:,:),leigvec(:,:),eigval(:)
 double precision, allocatable :: non_ortho_basis(:,:)
 integer :: n_real_eigv,sze,i,j,k,l
 double precision, allocatable :: S(:,:)
 sze = 10
 allocate(B(sze,sze),A(sze,sze),non_ortho_basis(sze,sze))
 allocate(reigvec(sze,sze),leigvec(sze,sze),eigval(sze))
 A = 0.d0
 B = 0.d0
 do i = 1, sze
  B(i,i) = 1.d0
  do j = 1, sze
   A(j,i) = h_non_hermit(j,i)
  enddo
 enddo
 print*,'With the generalized eigenvalue problem matrix '
 call non_hrmt_general_real_diag(sze,A,B,reigvec,leigvec,n_real_eigv,eigval)
 print*,''
 do i = 1, n_real_eigv
  print*,'i ',i,eigval(i)
 enddo
 allocate(S(n_real_eigv,n_real_eigv))
 do i = 1, n_real_eigv
  do j = 1, n_real_eigv
   S(j,i) = 0.d0
   do k = 1, sze
    S(j,i) += reigvec(k,i) * leigvec(k,j)
   enddo
  enddo
 enddo
 print*,'Left-Right Overlap matrix of the eigenstates '
 do i = 1, n_real_eigv
  write(*,'(100(F16.10,X))')S(i,:)
 enddo

 print*,'With the usual eigenvalue problem matrix '
 call non_hrmt_real_diag(sze,A,reigvec,leigvec,n_real_eigv,eigval)
 print*,''
 do i = 1, n_real_eigv
  print*,'i ',i,eigval(i)
 enddo

 print*,''
 print*,''
 print*,'Non orthonormal basis'
 print*,''
 non_ortho_basis = 0.d0
 do i = 1, sze
  do j = 1, sze
    if(i==j)then
     non_ortho_basis(j,i) = 1.d0
    else
     non_ortho_basis(j,i) = 1.d0/dabs(dble(i-j))
    endif
  enddo
 enddo
 A = 0.d0
 B = 0.d0
 do i = 1, sze
  do j = 1, sze
   do k = 1, sze
    B(j,i) += non_ortho_basis(k,i) * non_ortho_basis(k,j)
    do l = 1, sze
     A(j,i) += non_ortho_basis(k,i) * non_ortho_basis(l,k) * h_non_hermit(l,j)
    enddo
   enddo
  enddo
 enddo

 print*,'With the generalized eigenvalue problem matrix '
 print*,''
 print*,''
 print*,'Overlap matrix'
 print*,''
 do i = 1, sze
  write(*,'(100(F16.10,X))')B(i,:)
 enddo
 print*,''
 print*,''
 call non_hrmt_general_real_diag(sze,A,B,reigvec,leigvec,n_real_eigv,eigval)
 print*,''
 do i = 1, n_real_eigv
  print*,'i ',i,eigval(i)
 enddo
 do i = 1, n_real_eigv
  do j = 1, n_real_eigv
   S(j,i) = 0.d0
   do k = 1, sze
    do l = 1, sze
     S(j,i) += leigvec(l,j) * B(l,k) * reigvec(k,i) 
    enddo
   enddo
  enddo
 enddo
 print*,'Left-Right Overlap matrix of the eigenstates '
 do i = 1, n_real_eigv
  write(*,'(100(F16.10,X))')S(i,:)
 enddo
end
