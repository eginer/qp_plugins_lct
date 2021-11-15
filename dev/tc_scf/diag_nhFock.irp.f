 BEGIN_PROVIDER [ double precision, Reig_tcFock_ao, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Leig_tcFock_ao, (mo_num,mo_num) ]

  BEGIN_DOC
  ! left & right eigenvectors of the nh Fock matrix 
  ! in the |AO| basis obtained with level shift.
  END_DOC

  implicit none
  integer                       :: i, j, iorb, jorb
  double precision, allocatable :: F(:,:), Rvec(:,:), Lvec(:,:), diag(:)
  double precision, allocatable :: tmp(:,:)

  allocate( F(mo_num,mo_num) )
  F(1:mo_num,1:mo_num) = Fock_matrix_tc_mo_tot(1:mo_num,1:mo_num)

  allocate( diag(mo_num) , Rvec(mo_num,mo_num) , Lvec(mo_num,mo_num) )
  call real_nhdiag_RLeigvec_dgeev(F, mo_num, diag, Rvec, Lvec)
  deallocate( diag )

  ! Left.T x F
  allocate( tmp(mo_num,mo_num) )
  call dgemm( 'T', 'N', mo_num, mo_num, mo_num, 1.d0 &
            , Lvec, size(Lvec,1), F, size(F,1)       &
            , 0.d0, tmp, size(tmp,1) )
  deallocate( F , Lvec )

  ! F x Right
  allocate( F(mo_num,mo_num) )
  call dgemm( 'N', 'N', mo_num, mo_num, mo_num, 1.d0 &
            , tmp, size(tmp,1), Rvec, size(Rvec,1)   &
            , 0.d0, F, size(F,1) )
  deallocate( tmp , Rvec )

  ! Insert level shift here
  do i = elec_beta_num+1, elec_alpha_num
    F(i,i) = F(i,i) + 0.5d0 * level_shift
  enddo
  do i = elec_alpha_num+1, mo_num
    F(i,i) = F(i,i) + level_shift
  enddo

  allocate( diag(mo_num) , Rvec(mo_num,mo_num) , Lvec(mo_num,mo_num) )
  call real_nhdiag_RLeigvec_dgeev(F, mo_num, diag, Rvec, Lvec)

  deallocate(F, diag)

  Reig_tcFock_ao(1:mo_num,1:mo_num) = Rvec(1:mo_num,1:mo_num)
  Leig_tcFock_ao(1:mo_num,1:mo_num) = Lvec(1:mo_num,1:mo_num)

  deallocate( Rvec , Lvec )

END_PROVIDER



 BEGIN_PROVIDER [ double precision, Reig_tcFock_mo, (ao_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Leig_tcFock_mo, (ao_num,mo_num) ]

  BEGIN_DOC
  ! right eigenvectors of the nh Fock matrix 
  ! in the |MO| basis obtained with level shift.
  END_DOC

  implicit none

  call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0       &
            , mo_coef       , size(mo_coef,1)              &
            , Reig_tcFock_ao, size(Reig_tcFock_ao,1)       &
            , 0.d0, Reig_tcFock_mo, size(Reig_tcFock_mo,1) )

  call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0       &
            , mo_coef       , size(mo_coef,1)              &
            , Leig_tcFock_ao, size(Leig_tcFock_ao,1)       &
            , 0.d0, Leig_tcFock_mo, size(Leig_tcFock_mo,1) )

END_PROVIDER


!__________________________________________________________________________________________________
!
! DGEEV computes for an N-by-N real nonsymmetric matrix A, the eigenvalues and, optionally, 
! the left and/or right eigenvectors.
!
! The right eigenvector v(j) of A satisfies:     A * v(j) = lambda(j) * v(j)
! where lambda(j) is its eigenvalue.
! The left eigenvector u(j) of A satisfies :     u(j)**H * A = lambda(j) * u(j)**H
! where u(j)**H denotes the conjugate-transpose of u(j).
!
! The computed eigenvectors are normalized to have Euclidean norm equal to 1 
! and largest component real.
!
! Sven Hammarling:
!  Both dgeevx and dgeev compute eigenvalues and, optionally, eigenvectors. Both routines 
!  use the same underlying computational routines to achieve this so you should not see any 
!  difference in efficiency. But dgeevx can also compute condition numbers for you which can 
!  be very helpful in giving you information about the sensitivity of your problem. 
!  The down side of computing the eigenvector condition numbers is that additional N by N 
!  workspace is required, which only really matters if your matrices are very large.
!

subroutine real_nhdiag_Reigvec_dgeev(H, N, E, Rvec)

  implicit none

  integer         , intent(in)  :: N
  double precision, intent(in)  :: H(N,N)
  double precision, intent(out) :: E(N), Rvec(N,N)

  integer                       :: i, INFO, LWORK 
  double precision, allocatable :: A(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: WORK(:)

  allocate( A(N,N) )
  ! On entry, the N-by-N matrix A. On exit, A has been overwritten.
  A(1:N,1:N) = H(1:N,1:N)

  allocate( WR(N) , WI(N) )
  ! WR and WI contain the real and imaginary parts, respectively, of the 
  ! computed eigenvalues. Complex conjugate pairs of eigenvalues appear 
  ! consecutively with the eigenvalue having the positive imaginary part first.

  allocate( VL(N,N) , VR(N,N) )
  ! VL is not referenced.
  ! The right eigenvectors v(j) are stored one after another in the columns 
  ! of VR, in the same order as their eigenvalues.
  ! If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
  ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, 
  ! then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1)

  LWORK = -1
  ! If LWORK = -1, then a workspace query is assumed; the routine only calculates 
  ! the optimal size of the WORK array, returns this value as the first entry of 
  ! the WORK array, and no error message related to LWORK is issued by XERBLA.
  allocate( WORK(1) )
  ! On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

  call dgeev( "N", "V", N, A, size(A,1)              &
            , WR, WI, VL, size(VL,1), VR, size(VR,1) &
            , WORK, LWORK, INFO) 

  if( INFO .ne. 0 ) then
    print *, ' Problem is subroutine real_nhdiag_Reigvec_dgeev'
    print *, ' INFO != 0 in the first call of dgeev'
    print *, ' INFO  =', INFO
    stop
  endif

  LWORK = WORK(1)
  ! LWORK is the dimension of the array WORK. 
  ! LWORK >= 4*N. For good performance, LWORK must generally be larger.
  deallocate( WORK )
  allocate( WORK(LWORK) )
  call dgeev( "N", "V", N, A, size(A,1)              &
            , WR, WI, VL, size(VL,1), VR, size(VR,1) &
            , WORK, LWORK, INFO) 

  if( INFO .ne. 0 ) then
    print *, ' Problem is subroutine real_nhdiag_Reigvec_dgeev'
    print *, ' INFO != 0 in the second call of dgeev'
    print *, ' INFO  =', INFO
    stop
  endif

  deallocate( A , WORK )
  deallocate( VL )
  
  do i = 1, N
    if( dabs(WI(i)) .gt. (1d-5) ) then
      print *, ' Problem is subroutine real_nhdiag_Reigvec_dgeev'
      print *, ' non real eigenvalue !'
      print *, i, WI(i)
      stop
    endif
    E   (    i) = WR(    i)
    Rvec(1:N,i) = VR(1:N,i)
  enddo

  deallocate( WR , WI , VR )

  return
end subroutine real_nhdiag_Reigvec_dgeev


subroutine real_nhdiag_RLeigvec_dgeev(H, N, E, Rvec, Lvec)

  implicit none

  integer         , intent(in)  :: N
  double precision, intent(in)  :: H(N,N)
  double precision, intent(out) :: E(N), Rvec(N,N), Lvec(N,N)

  integer                       :: i, INFO, LWORK 
  double precision, allocatable :: A(:,:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: WORK(:)

  allocate( A(N,N) )
  ! On entry, the N-by-N matrix A. On exit, A has been overwritten.
  A(1:N,1:N) = H(1:N,1:N)

  allocate( WR(N) , WI(N) )
  ! WR and WI contain the real and imaginary parts, respectively, of the 
  ! computed eigenvalues. Complex conjugate pairs of eigenvalues appear 
  ! consecutively with the eigenvalue having the positive imaginary part first.

  allocate( VL(N,N) , VR(N,N) )
  ! The left eigenvectors u(j) are stored one after another in the columns 
  ! of VL, in the same order as their eigenvalues.
  ! If the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column of VL.
  ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, 
  ! then u(j) = VL(:,j) + i*VL(:,j+1) and  u(j+1) = VL(:,j) - i*VL(:,j+1).
  !
  ! The right eigenvectors v(j) are stored one after another in the columns 
  ! of VR, in the same order as their eigenvalues.
  ! If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
  ! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, 
  ! then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1)

  LWORK = -1
  ! If LWORK = -1, then a workspace query is assumed; the routine only calculates 
  ! the optimal size of the WORK array, returns this value as the first entry of 
  ! the WORK array, and no error message related to LWORK is issued by XERBLA.
  allocate( WORK(1) )
  ! On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

  call dgeev( "V", "V", N, A, size(A,1)              &
            , WR, WI, VL, size(VL,1), VR, size(VR,1) &
            , WORK, LWORK, INFO) 

  if( INFO .ne. 0 ) then
    print *, ' Problem is subroutine real_nhdiag_RLeigvec_dgeev'
    print *, ' INFO != 0 in the first call of dgeev'
    print *, ' INFO  =', INFO
    stop
  endif

  LWORK = WORK(1)
  ! LWORK is the dimension of the array WORK. 
  ! LWORK >= 4*N. For good performance, LWORK must generally be larger.
  deallocate( WORK )
  allocate( WORK(LWORK) )
  call dgeev( "V", "V", N, A, size(A,1)              &
            , WR, WI, VL, size(VL,1), VR, size(VR,1) &
            , WORK, LWORK, INFO) 

  if( INFO .ne. 0 ) then
    print *, ' Problem is subroutine real_nhdiag_RLeigvec_dgeev'
    print *, ' INFO != 0 in the second call of dgeev'
    print *, ' INFO  =', INFO
    stop
  endif

  deallocate( A , WORK )
  
  do i = 1, N
    if( dabs(WI(i)) .gt. (1d-5) ) then
      print *, ' Problem is subroutine real_nhdiag_RLeigvec_dgeev'
      print *, ' non real eigenvalue !'
      print *, i, WI(i)
      stop
    endif
    E   (    i) = WR(    i)
    Rvec(1:N,i) = VR(1:N,i)
    Lvec(1:N,i) = VL(1:N,i)
  enddo

  deallocate( WR , WI , VR , VL )

  return
end subroutine real_nhdiag_RLeigvec_dgeev
!__________________________________________________________________________________________________
