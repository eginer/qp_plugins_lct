
subroutine get_inv_half_svd(matrix, n, matrix_inv_half)

  BEGIN_DOC
  !   :math:`X = S^{-1/2}` obtained by SVD
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: matrix(n,n)
  double precision, intent(out) :: matrix_inv_half(n,n)

  integer                       :: num_linear_dependencies
  integer                       :: LDA, LDC
  integer                       :: info, i, j, k
  double precision, parameter   :: threshold = 1.d-6
  double precision, allocatable :: U(:,:),Vt(:,:), D(:)

  LDA = size(matrix, 1)
  LDC = size(matrix_inv_half, 1)
  if(LDA .ne. LDC) then
    print*, ' LDA != LDC'
    stop
  endif

  print*, ' n   = ', n
  print*, ' LDA = ', LDA
  print*, ' LDC = ', LDC

  allocate(U(LDC,n), Vt(LDA,n), D(n))

  call svd(matrix, LDA, U, LDC, D, Vt, LDA, n, n)

  print*, ' U'
  do i = 1, n
   write(*, '(100(F16.10,x))') U(:,i)
  enddo

  print*, ' Vt'
  do i = 1, n
   write(*, '(100(F16.10,x))') Vt(:,i)
  enddo

  print*, ' D'
  do i = 1,n
   write(*, '(100(F16.10,x))') D(i)
  enddo

  num_linear_dependencies = 0
  do i = 1, n
    if(abs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  write(*,*) ' linear dependencies', num_linear_dependencies

  matrix_inv_half = 0.d0
  do k = 1, n
    if(D(k) /= 0.d0) then
      do j = 1, n
        do i = 1, n
          matrix_inv_half(i,j) = matrix_inv_half(i,j) + U(i,k) * D(k) * Vt(k,j)
          !matrix_inv_half(i,j) = matrix_inv_half(i,j) + U(i,k) * D(k) * Vt(j,k)
        enddo
      enddo
    endif
  enddo

  print*,'S inv half'
  do i = 1, n
   write(*, '(1000(F16.10,X))') matrix_inv_half(i,:)
  enddo

  double precision, allocatable :: Stmp(:,:), Stmp2(:,:)
  allocate( Stmp(n,n), Stmp2(n,n) )
 
  ! S^-1/2 x S
  Stmp = 0.d0
  call dgemm( 'N', 'N', n, n, n, 1.d0                                            &
            , matrix_inv_half, size(matrix_inv_half, 1), matrix, size(matrix, 1) &
            , 0.d0, Stmp, size(Stmp, 1) )

  ! ( S^-1/2 x S ) x S^-1/2
  !Stmp2 = 0.d0
  !call dgemm( 'N', 'N', n, n, n, 1.d0                                        &
  !          , Stmp, size(Stmp, 1), matrix_inv_half, size(matrix_inv_half, 1) &
  !          , 0.d0, Stmp2, size(Stmp2, 1) )

  ! S^-1/2 x ( S^-1/2 x S )
  Stmp2 = 0.d0
  call dgemm( 'N', 'N', n, n, n, 1.d0                                        &
            , matrix_inv_half, size(matrix_inv_half, 1), Stmp, size(Stmp, 1) &
            , 0.d0, Stmp2, size(Stmp2, 1) )
 
  double precision :: accu_d,accu_nd
  accu_nd = 0.d0
  accu_d  = 0.d0
  do i = 1, n
    do j = 1, n
      if(i==j) then
       accu_d += Stmp2(j,i)
      else 
       accu_nd = accu_nd + Stmp2(j,i) * Stmp2(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)
  print*, ' after S^-1/2: sum of off-diag S elements = ', accu_nd
  print*, ' after S^-1/2: sum of     diag S elements = ', accu_d
  do i = 1, n
   write(*,'(1000(F16.10,X))') Stmp2(i,:)
  enddo

  !double precision  :: thresh
  !thresh = 1.d-10
  !if( accu_nd.gt.thresh .or. dabs(accu_d-dble(n)).gt.thresh) then
  !  stop
  !endif

end subroutine get_inv_half_svd

! ---

