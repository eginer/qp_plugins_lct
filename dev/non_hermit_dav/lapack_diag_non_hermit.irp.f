subroutine lapack_diag_non_sym(n,A,WR,WI,VL,VR)
 implicit none
  BEGIN_DOC
! You enter with a general non hermitian matrix A(n,n) 
!
! You get out with the real WR and imaginary part WI of the eigenvalues 
!
! Eigvalue(n) = WR(n) + i * WI(n)
!
! And the left VL and right VR eigenvectors 
!
! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
!
! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n)
 double precision, intent(out):: WR(n),WI(n),VL(n,n),VR(n,n)
  character*1 :: JOBVL,JOBVR
  JOBVL = "V" ! computes the left  eigenvectors 
  JOBVR = "V" ! computes the right eigenvectors 
  integer     :: lda,ldvl,ldvr,LWORK,INFO
  double precision, allocatable :: WORK(:)
  integer :: n_good
  lda  = n
  ldvl = n
  ldvr = n
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
  call dgeev('V','V',n,A,lda,WR,WI,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
  if(INFO.ne.0)then
   print*,'dgeev failed !!',INFO
   stop
  endif
end 


subroutine non_hrmt_real_diag(n,A,reigvec,leigvec,n_real_eigv,eigval)
 implicit none
 BEGIN_DOC
! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
!
! of a non hermitian matrix A(n,n)
!
! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
 END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n)
 double precision, intent(out) :: reigvec(n,n),leigvec(n,n),eigval(n)
 double precision, allocatable :: Aw(:,:)
 integer, intent(out) :: n_real_eigv
 print*,'Computing the left/right eigenvectors ...'
 character*1 :: JOBVL,JOBVR
 JOBVL = "V" ! computes the left  eigenvectors 
 JOBVR = "V" ! computes the right eigenvectors 
 double precision, allocatable :: WR(:),WI(:),Vl(:,:),VR(:,:)
 integer :: i,j,k
 integer :: n_good
 integer, allocatable :: list_good(:), iorder(:)
 double precision :: thr
 thr = 1.d-5
 ! Eigvalue(n) = WR(n) + i * WI(n)
 allocate(WR(n),WI(n),VL(n,n),VR(n,n),Aw(n,n))
 Aw = A
 call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
 ! You track the real eigenvalues 
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.thr)then
   n_good += 1
  else
   print*,'Found an imaginary component to eigenvalue'
   print*,'Re(i) + Im(i)',WR(i),WI(i)
  endif
 enddo
 allocate(list_good(n_good),iorder(n_good))
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.thr)then
   n_good += 1
   list_good(n_good) = i
   eigval(n_good) = WR(i)
  endif
 enddo
 n_real_eigv = n_good 
 do i = 1, n_good
  iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(eigval,iorder,n_good)
 print*,'n_real_eigv = ',n_real_eigv
 print*,'n           = ',n
 do i = 1, n_real_eigv
  print*,i,'eigval(i) = ',eigval(i) 
  do j = 1, n
   reigvec(j,i) = VR(j,list_good(iorder(i)))
   leigvec(j,i) = Vl(j,list_good(iorder(i)))
  enddo
 enddo
end

subroutine lapack_diag_general_non_sym(n,A,B,WR,beta,WI,VL,VR)
 implicit none
  BEGIN_DOC
! You enter with a general non hermitian matrix A(n,n) and another B(n,n)
!
! You get out with the real WR and imaginary part WI of the eigenvalues 
!
! Eigvalue(n) = (WR(n) + i * WI(n))/(beta(n))
!
! And the left VL and right VR eigenvectors 
!
! VL(i,j) = <i|Psi_left(j)>  :: projection on the basis element |i> on the jth left  eigenvector 
!
! VR(i,j) = <i|Psi_right(j)> :: projection on the basis element |i> on the jth right eigenvector 
  END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n),B(n,n)
 double precision, intent(out):: WR(n),WI(n),beta(n),VL(n,n),VR(n,n)
  character*1 :: JOBVL,JOBVR
  JOBVL = "V" ! computes the left  eigenvectors 
  JOBVR = "V" ! computes the right eigenvectors 
  integer     :: lda,ldvl,ldvr,LWORK,INFO
  double precision, allocatable :: WORK(:)
  integer :: n_good
  lda  = n
  ldvl = n
  ldvr = n
  allocate(WORK(1))
  LWORK = -1 ! to ask for the optimal size of WORK
  call dggev('V','V',n,A,lda,B,lda,WR,WI,beta,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
  if(INFO.gt.0)then
   print*,'dgeev failed !!',INFO
   stop
  endif
  LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK 
  deallocate(WORK)
  allocate(WORK(LWORK))
  ! Actual diagonalization 
  call dggev('V','V',n,A,lda,B,lda,WR,WI,beta,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
  if(INFO.ne.0)then
   print*,'dgeev failed !!',INFO
   stop
  endif
end 

subroutine non_hrmt_general_real_diag(n,A,B,reigvec,leigvec,n_real_eigv,eigval)
 implicit none
 BEGIN_DOC
! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
!
! of a non hermitian matrix A(n,n) and B(n,n) 
!
! A reigvec = eigval * B * reigvec
!
! (A)^\dagger leigvec = eigval * B * leigvec
!
! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
 END_DOC
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n),B(n,n)
 double precision, intent(out) :: reigvec(n,n),leigvec(n,n),eigval(n)
 double precision, allocatable :: Aw(:,:),Bw(:,:)
 integer, intent(out) :: n_real_eigv
 print*,'Computing the left/right eigenvectors ...'
 character*1 :: JOBVL,JOBVR
 JOBVL = "V" ! computes the left  eigenvectors 
 JOBVR = "V" ! computes the right eigenvectors 
 double precision, allocatable :: WR(:),WI(:),Vl(:,:),VR(:,:),beta(:)
 integer :: i,j,k
 integer :: n_good
 integer, allocatable :: list_good(:), iorder(:)
 ! Eigvalue(n) = WR(n) + i * WI(n)
 allocate(WR(n),WI(n),VL(n,n),VR(n,n),Aw(n,n),beta(n),Bw(n,n))
 Aw = A
 Bw = B
 call lapack_diag_general_non_sym(n,A,B,WR,beta,WI,VL,VR)
 ! You track the real eigenvalues 
 n_good = 0
 do i = 1, n
  print*,'beta(i) = ',beta(i)
  if(dabs(WI(i)).lt.1.d-12)then
   n_good += 1
  else
!   print*,'Found an imaginary component to eigenvalue'
!   print*,'Re(i) + Im(i)',WR(i),WI(i)
  endif
 enddo
 allocate(list_good(n_good),iorder(n_good))
 n_good = 0
 do i = 1, n
  if(dabs(WI(i)).lt.1.d-12)then
   n_good += 1
   list_good(n_good) = i
   eigval(n_good) = WR(i)/beta(i)
  endif
 enddo
 n_real_eigv = n_good 
 do i = 1, n_good
  iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(eigval,iorder,n_good)
 print*,'n_real_eigv = ',n_real_eigv
 print*,'n           = ',n
 do i = 1, n_real_eigv
  print*,i,'eigval(i) = ',eigval(i) 
  do j = 1, n
   reigvec(j,i) = VR(j,list_good(iorder(i)))
   leigvec(j,i) = Vl(j,list_good(iorder(i)))
  enddo
 enddo
end
