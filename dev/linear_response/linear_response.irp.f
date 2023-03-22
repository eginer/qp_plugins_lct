program linear_response
  implicit none
  BEGIN_DOC
  ! call the diagonlize_lr routine to solve the linear response equation
  END_DOC
  read_wf = .true.
  touch read_wf

  call diagonalize_lr 
end


subroutine diagonalize_lr
 implicit none
 BEGIN_DOC
 ! Solve the generalized eigenvalue equation M x = lambda * N * x where 
 ! - M is the dim*dim (dim = 2*(N_det-1)) matrix where the diagonal blocks (Ndet-1)*(Ndet-1) are the A matrix
 !   and the off diagonal blocks (Ndet-1)*(Ndet-1) are the B matrix.
 ! - N is the dim*dim matrix where the diagonal blocks Ndet*Ndet are the overlap (and minus the overlap) matrix
 !   and the off-diagonal blocks are zero.
 ! - x are the eigenvectors with associated eigenvalues lambda.
 !
 END_DOC
 use omp_lib
 
 double precision, allocatable :: AB_matrix(:,:), S_matrix(:,:)
 
 integer :: det_I, det_J, dim, dim_half, i, j

 double precision, allocatable :: alphar(:), alphai(:), beta(:), mat_a(:,:), mat_b(:,:)
 double precision, allocatable :: eigenvectors_left(:,:), eigenvectors(:,:), eigenvalues(:) 
 double precision, allocatable :: work(:)
 integer :: lwork, info
 
 double precision :: cpu_time_1, cpu_time_2
 double precision :: omp_time_1, omp_time_2
 integer, allocatable ::  ipiv(:)
 integer :: n,istate
 
 !###############################################################
 print*, "Allocate matrices"
 call cpu_time (cpu_time_1)
 omp_time_1 = OMP_get_wtime()

 dim_half = N_det - 1
 dim = 2*dim_half
 lwork = 10*dim

 allocate(ipiv(dim))
 allocate(work(lwork))
 allocate(eigenvectors_left(dim,dim),eigenvectors(dim,dim),&
         &alphar(dim),alphai(dim), beta(dim))

 allocate(mat_a(dim, dim), mat_b(dim, dim), AB_matrix(dim, dim), S_matrix(dim, dim))
 allocate(eigenvalues(dim))
 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Allocate matrices - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Allocate matrices - WALLCLOCK TIME =", omp_time_2 - omp_time_1


 !###############################################################
 print*, "Fill S matrix + inversion"
 call cpu_time (cpu_time_1)
 omp_time_1 = omp_get_wtime()

 S_matrix = 0.d0

 do det_I = 1, dim_half 
  do det_J = det_I, dim_half 
    S_matrix (det_I, det_J) = S_IJ(det_I, det_J, 1)  
    S_matrix (det_J, det_I) = S_matrix (det_I, det_J)
    
    S_matrix (det_I + dim_half, det_J + dim_half) = - S_matrix (det_I, det_J)
    S_matrix (det_J + dim_half, det_I + dim_half) = S_matrix (det_I + dim_half, det_J + dim_half)
  enddo
 enddo

 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Fill S matrix - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Fill S matrix - WALLCLOCK TIME =", omp_time_2 - omp_time_1

 !###############################################################
 print*, "Fill AB matrix"
 call cpu_time (cpu_time_1)
 omp_time_1 = omp_get_wtime()

 AB_matrix = 0.d0
 do det_I = 1, dim_half 
  do det_J = 1, dim_half 
    AB_matrix (det_I, det_J) = A_IJ(det_I, det_J, 1)
    AB_matrix (det_I + dim_half, det_J + dim_half) = AB_matrix (det_I, det_J)
     
    AB_matrix (det_I, det_J + dim_half) = B_IJ(det_I, det_J, 1) 
    AB_matrix (det_I + dim_half, det_J) = AB_matrix (det_I, det_J + dim_half)
    
  enddo
 enddo

 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Fill AB matrix - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Fill AB matrix - WALLCLOCK TIME =", omp_time_2 - omp_time_1

 !###############################################################
 print*, "Eigenvalue problem"
 call cpu_time (cpu_time_1)
 omp_time_1 = omp_get_wtime()
! print*,"AB matrix"
! do det_I = 1, dim_half
!  write(*,'(100(f10.5,x))'), AB_matrix(1:dim_half, det_I)
! enddo
! print*, "S matrix"
! do det_I = 1, dim_half
!  write(*,'(100(f10.5,x))'), S_matrix(1:dim_half, det_I)
! enddo


  mat_a = AB_matrix
  mat_b = S_matrix
  call dggev("N", "V", dim, mat_a, dim, mat_b, dim, alphar, alphai, beta, &
  & eigenvectors_left,dim,eigenvectors,dim, work, lwork, info)

 deallocate(mat_a, mat_b, AB_matrix, S_matrix)

 eigenvalues = alphar/beta
 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Eigenvalue problem - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Eigenvalue problem - WALLCLOCK TIME =", omp_time_2 - omp_time_1

!print*, "Eigenvalues (alphareal/beta)"
!do det_I=1, dim
! print*, eigenvalues(det_I)
!enddo

 integer, allocatable :: iorder(:)
 allocate(iorder(dim))
 double precision, allocatable :: eigenvectors_sort(:,:), eigenvectors_transp(:,:)
 allocate(eigenvectors_sort(dim, dim))
 allocate(eigenvectors_transp(dim_half, dim_half))
 do i = 1, dim
   iorder(i) = i
 enddo
 call dsort(eigenvalues,iorder,dim) ! sort eigenvalues : dim_half of negatives , dim_half of positives in mirror
!call dlasrt("I", dim, eigenvalues,info)
 print*, "Eigenvalues (alphareal/beta)"
 do det_I=1, dim
  print*, eigenvalues(det_I)
 enddo
 do i = 1, dim
  do j = 1, dim
   ! the first dim_half are of negative eigenvalues, the second dim_half are the one interesting 
   eigenvectors_sort(j,i) = eigenvectors(j,iorder(i)) ! sort the corresponding eigenvectors 
  enddo
 enddo
!print*,'eigenvectors '
!do i = 1, dim
! write(*,'(I2,X,100(F8.5,X))')i,eigenvectors_sort(i,1:dim)
!enddo
  do i = 1, dim_half
   do j = 1, dim_half
    ! eigenvectors(i,j) = <det_i|Psi_j>
    ! eigenvectors_transp(j,i) = <Psi_j| det_i> 
    eigenvectors_transp(j,i) = eigenvectors_sort(i,j+dim_half) ! we obtain the interesting one
   enddo
  enddo
!print*,'eigenvectors transp'
!do i = 1, dim_half
! write(*,'(I2,X,100(F8.5,X))')i,eigenvectors_transp(i,1:dim_half)
!enddo
 
  !! Computation of S2
  double precision, allocatable :: S2_values_lr(:), norm_vec(:)
  double precision :: s2_ij,s2_ii
  allocate(S2_values_lr(dim_half), norm_vec(dim_half))
  S2_values_lr = 0.d0
  norm_vec = 0.d0
  i=0
  do det_I = 1, N_det
   if(det_I == hf_index) cycle
   i+= 1
!  call get_s2(psi_det(1,1,det_I),psi_det(1,1,det_I),N_int,s2_ii)
   do istate = 1, dim_half
!   S2_values_lr(istate) += s2_ii * eigenvectors_transp(istate,i) * eigenvectors_transp(istate,i)
    norm_vec(istate) += eigenvectors_transp(istate,i) * eigenvectors_transp(istate,i)
   enddo
   j = 0
   do det_J = 1, N_det
    if(det_J == hf_index) cycle
     j += 1
     call get_s2(psi_det(1,1,det_I),psi_det(1,1,det_J),N_int,s2_ij)
     ! s2 = <det_i| S^2 | det_j>
     do istate = 1, dim_half
      S2_values_lr(istate) += 1.d0 * s2_ij * eigenvectors_transp(istate,i) * eigenvectors_transp(istate,j)
     enddo
   enddo
  enddo
  do i = 1, dim_half
   S2_values_lr(i) *= 1.d0/(norm_vec(i))
  enddo
 !do istate = 1, dim_half
 ! print*,'S2_values_lr',istate,S2_values_lr(istate)
 !enddo
 print*, " Half eigenvalues "
 do istate=1, dim_half
  print*, eigenvalues(istate+dim_half), S2_values_lr(istate)
 enddo


 deallocate(eigenvectors_left, eigenvectors)
 deallocate(work, alphar, alphai, beta)
end subroutine


