
 BEGIN_PROVIDER [double precision, one_e_j_eigval_trans, (mo_num)]
&BEGIN_PROVIDER [integer, n_good_one_e_trans_eigval]
&BEGIN_PROVIDER [double precision, one_e_j_reigvec_trans, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, one_e_j_leigvec_trans, (mo_num,mo_num)]
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
 n = mo_num
 allocate(A(n,n),WR(n),WI(n),VL(n,n),VR(n,n))
 lda  = n
 ldvl = n
 ldvr = n
 A = h_tilde_one_j
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
 A = h_tilde_one_j
 call dgeev('V','V',n,A,lda,WR,WI,VL,ldvl,VR,ldvr,WORK,LWORK,INFO)
 if(INFO.ne.0)then
  print*,'dgeev failed !!',INFO
  stop
 endif

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
 n_good_one_e_trans_eigval = n_good 
 do i = 1, n_good
  iorder(i) = i
 enddo
 ! You sort the real eigenvalues 
 call dsort(ei,iorder,n_good)
 double precision :: accu1, accu2
 do i = 1, n_good_one_e_trans_eigval
  one_e_j_eigval_trans(i) = ei(i)
  accu1 = 0.d0
  accu2 = 0.d0
  do j = 1, mo_num
   one_e_j_reigvec_trans(j,i) = VR(j,list_good(iorder(i)))
   one_e_j_leigvec_trans(j,i) = Vl(j,list_good(iorder(i)))
   accu1 += one_e_j_reigvec_trans(j,i) * one_e_j_reigvec_trans(j,i)
   accu2 += one_e_j_leigvec_trans(j,i) * one_e_j_leigvec_trans(j,i)
  enddo
!  one_e_j_reigvec_trans_norm(i) = accu1
!  one_e_j_leigvec_trans_norm(i) = accu2
 enddo
END_PROVIDER 
