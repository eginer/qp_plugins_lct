program linear_response_k_equal_0
  implicit none
  BEGIN_DOC
  ! call the diagonlize_lr routine to solve the linear response equation
  END_DOC
  read_wf = .true.
  touch read_wf

  call diagonalize_lr_k_equal_0 
end


subroutine diagonalize_lr_k_equal_0
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
 integer :: n
 
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
    AB_matrix (det_I, det_J) = A_IJ_K_equal_0(det_I, det_J, 1)
    AB_matrix (det_I + dim_half, det_J + dim_half) = AB_matrix (det_I, det_J)
     
    AB_matrix (det_I, det_J + dim_half) = B_IJ_K_equal_0(det_I, det_J, 1) 
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

  print*,"AB matrix"
  do det_I = 1, dim_half
   write(*,'(100(f10.5,x))'), AB_matrix(1:dim_half, det_I)
  enddo
  print*, "S matrix"
  do det_I = 1, dim_half
   write(*,'(100(f10.5,x))'), S_matrix(1:dim_half, det_I)
  enddo

  mat_a = AB_matrix
  mat_b = S_matrix
  call dggev("N", "V", dim, mat_a, dim, mat_b, dim, alphar, alphai, beta, &
  & eigenvectors_left,dim,eigenvectors,dim, work, lwork, info)
 
 deallocate(eigenvectors_left, eigenvectors)
 deallocate(mat_a, mat_b, AB_matrix, S_matrix)

 eigenvalues = alphar/beta
 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Eigenvalue problem - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Eigenvalue problem - WALLCLOCK TIME =", omp_time_2 - omp_time_1

call dlasrt("I", dim, eigenvalues,info)
 print*, "Eigenvalues (alphareal/beta)"
 do det_I=1, dim
  print*, eigenvalues(det_I)
 enddo

 deallocate(work, alphar, alphai, beta)
end subroutine


