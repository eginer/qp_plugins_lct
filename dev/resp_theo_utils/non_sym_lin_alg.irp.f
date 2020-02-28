subroutine non_sym_lapack_diagd(eigvalues_r,eigvalues_i,eigvectors_r,H,nmax,n)
  implicit none
  BEGIN_DOC
  ! Computes the right eigenvectors of matrix H
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n,nmax
  double precision, intent(out)  :: eigvectors_r(nmax,n)
  double precision, intent(out)  :: eigvalues_r(n),eigvalues_i(n)
  double precision, intent(in)   :: H(nmax,n)
  double precision,allocatable   :: eigenvectors_l(:,:)
  double precision,allocatable   :: work(:)
  double precision,allocatable   :: A(:,:)
  integer                        :: lwork, info, i,j,l,k, liwork

  allocate(A(nmax,n))

  A=H
  lwork = max(1000,2*n*n + 6*n+ 1)
  allocate (work(lwork),eigenvectors_l(n,n))

  lwork = -1
  call dgeev( 'N', 'V', n, A, nmax, eigvalues_r, eigvalues_i, eigenvectors_l, n, eigvectors_r, nmax,  work, lwork, info )
  if (info < 0) then
    print *, irp_here, ': DGEEV: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  lwork  = int( work( 1 ) )
  deallocate (work)

  allocate (work(lwork))
  call dgeev( 'N', 'V', n, A, nmax, eigvalues_r, eigvalues_i, eigenvectors_l, n, eigvectors_r, nmax,  work, lwork, info )
  deallocate(work)

  if (info < 0) then
    print *, irp_here, ': DGEEV: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0  ) then
    write(*,*)'DGEEV Failed'
    stop 1
  end if
  deallocate(A)
end

