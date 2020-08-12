BEGIN_PROVIDER [double precision, htilde_matrix_elmt, (N_det,N_det)]
 implicit none
 integer :: i,j
 double precision :: hmono,herf,heff,hderiv,htot
 do i = 1, N_det
  do j = 1, N_det
  ! < J | Htilde | I >
!   call htilde_mat(psi_det(1,1,j),psi_det(1,1,i),hmono,herf,heff,hderiv,htot)
   call htilde_mat(psi_det(1,1,i),psi_det(1,1,j),hmono,herf,heff,hderiv,htot)
   htilde_matrix_elmt(j,i) = htot
  enddo
 enddo
! htilde_matrix_elmt = H_matrix_all_dets
END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigval_trans, (N_det)]
&BEGIN_PROVIDER [integer, n_good_trans_eigval]
&BEGIN_PROVIDER [double precision, reigvec_trans, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, leigvec_trans, (N_det,N_det)]
&BEGIN_PROVIDER [double precision, reigvec_trans_norm, (N_det)]
&BEGIN_PROVIDER [double precision, leigvec_trans_norm, (N_det)]
 implicit none
 character*1 :: JOBVL,JOBVR
 JOBVL = "V" ! computes the left  eigenvectors 
 JOBVR = "V" ! computes the right eigenvectors 
 integer     :: n,lda,ldvl,ldvr,LWORK,INFO
 double precision, allocatable :: A(:,:),WR(:),WI(:),Vl(:,:),VR(:,:)
 double precision, allocatable :: WORK(:)
 integer :: i,j,k
 integer :: n_good
 integer, allocatable :: list_good(:), iorder(:)
 double precision, allocatable :: ei(:)
 ! Eigvalue(n) = WR(n) + i * WI(n)
 n = n_det
 allocate(A(n,n),WR(n),WI(n),VL(n,n),VR(n,n))
 lda  = n
 ldvl = n
 ldvr = n
 A = htilde_matrix_elmt
 allocate(WORK(1))
 LWORK = -1 ! to ask for the optimal size of WORK
 call dgeev('V','V',n,A,lda,WR,WI,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
 if(INFO.gt.0)then
  print*,'dgeev failed !!',INFO
  stop
 endif
 LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK 
 deallocate(WORK)
 allocate(WORK(LWORK))
 ! Actual diagonalization 
 A = htilde_matrix_elmt
 call dgeev('V','V',n,A,lda,WR,WI,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
 if(INFO.ne.0)then
  print*,'dgeev failed !!',INFO
  stop
 endif
! double precision :: accu
! double precision, allocatable :: matrix(:,:)
! allocate(matrix(n,n))
! do i = 1, n
!  do j = 1, n
!   accu = 0.d0
!   do k = 1, n
!    accu += VL(k,j) * VR(k,i)
!   enddo
!   matrix(j,i) = accu
!  enddo
! enddo
! print*,'Overlap L/R matrix'
! do i = 1, n
!  write(*,'(1000(F16.10,X))')matrix(i,:)
! enddo

 ! You track the real eigenvalues 
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.1.d-12)then
   n_good += 1
  endif
 enddo
 allocate(list_good(n_good),iorder(n_good),ei(n_good))
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.1.d-12)then
   n_good += 1
   list_good(n_good) = i
   ei(n_good) = WR(i)
  endif
 enddo
 n_good_trans_eigval = n_good 
 do i = 1, n_good
  iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(ei,iorder,n_good)
 double precision :: accu1, accu2
 do i = 1, n_good_trans_eigval
  eigval_trans(i) = ei(i)
  accu1 = 0.d0
  accu2 = 0.d0
  do j = 1, n_det
   reigvec_trans(j,i) = VR(j,list_good(iorder(i)))
   leigvec_trans(j,i) = Vl(j,list_good(iorder(i)))
   accu1 += reigvec_trans(j,i) * reigvec_trans(j,i)
   accu2 += leigvec_trans(j,i) * leigvec_trans(j,i)
  enddo
  reigvec_trans_norm(i) = accu1
  leigvec_trans_norm(i) = accu2
 enddo

END_PROVIDER 
