subroutine lapack_diag_non_sym(n,A,WR,WI,VL,VR)
 implicit none
 integer, intent(in) :: n
 double precision, intent(in) :: A(n,n)
 double precision, intent(out):: WR(n),WI(n),VL(n,n),VR(n,n)
  print*,'Computing the left/right eigenvectors ...'
  character*1 :: JOBVL,JOBVR
  JOBVL = "V" ! computes the left  eigenvectors 
  JOBVR = "V" ! computes the right eigenvectors 
  integer     :: lda,ldvl,ldvr,LWORK,INFO
  double precision, allocatable :: WORK(:)
  integer :: n_good
  ! Eigvalue(n) = WR(n) + i * WI(n)
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
