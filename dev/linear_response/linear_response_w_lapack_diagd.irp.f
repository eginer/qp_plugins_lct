program linear_response_w_lapack_diagd
  implicit none
  BEGIN_DOC
  ! call the diagonlize_lr routine to solve the linear response equation
  END_DOC
  read_wf = .true.
  touch read_wf

  print*, "WARNING : the script of the routine associated with this program &
  & has not been updated with major corrections needed in the original program &
  & linear_response."
  call diagonalize_lr_w_lapack_diagd 
end

subroutine diagonalize_lr_w_lapack_diagd
 implicit none
 BEGIN_DOC
 ! Solve the generalized eigenvalue equation M x = lambda * N * x where 
 ! - M is the dim*dim (dim = 2*N_det) matrix where the diagonal blocks Ndet*Ndet are the A matrix
 !   and the off diagonal blocks Ndet*det are the B matrix.
 ! - N is the dim*dim matrix where the diagonal blocks Ndet*Ndet are the overlap (and minus the overlap) matrix
 !   and the off-diagonal blocks are zero.
 ! - x are the eigenvectors with associated eigenvalues lambda.
 !
 ! Because N is invertible, we solve the following standard eigenvalue problem
 ! Mprim * x = lambda * x where Mprim = N^(-1) * M.
 ! Thus, we use the utils/lapack_diag routine.
 !
 END_DOC
 use omp_lib
 
 double precision, allocatable :: AB_matrix(:,:), S_matrix(:,:)
 double precision, allocatable :: copy_AB_matrix(:,:), copy_S_matrix(:,:), S_inv_time_AB(:,:)
 
 integer :: det_I, det_J, dim, i, j
 double precision, allocatable :: eig_omega_lr(:) 

 double precision, allocatable :: alphar(:), alphai(:), beta(:), mat_a(:,:), mat_b(:,:)
 double precision, allocatable :: eigenvectors_left(:,:), eigenvectors(:,:) 
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

 dim = 2*N_det
 lwork = 10*dim

 allocate(ipiv(dim))
 allocate(work(lwork))
 allocate(eigenvectors_left(dim,dim),eigenvectors(dim,dim),&
         &eig_omega_lr(dim),alphar(dim),alphai(dim), beta(dim))

 allocate(mat_a(dim, dim), mat_b(dim, dim), AB_matrix(dim, dim), S_matrix(dim, dim))
 allocate(copy_AB_matrix(dim,dim), copy_S_matrix(dim,dim),S_inv_time_AB(dim,dim))
 
 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Allocate matrices - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Allocate matrices - WALLCLOCK TIME =", omp_time_2 - omp_time_1


 !###############################################################
 print*, "Fill S matrix + inversion"
 call cpu_time (cpu_time_1)
 omp_time_1 = omp_get_wtime()

 AB_matrix = 0.d0
 S_matrix = 0.d0

 do det_I = 1, N_det 
  do det_J = det_I, N_det
    S_matrix (det_I, det_J) = S_IJ(det_I, det_J, 1)  
    S_matrix (det_J, det_I) = S_matrix (det_I, det_J)
    
    S_matrix (det_I + N_det, det_J + N_det) = - S_matrix (det_I, det_J)
    S_matrix (det_J + N_det, det_I + N_det) = S_matrix (det_I + N_det, det_J + N_det)
  enddo
 enddo

 ! Inversion of S_matrix. Inversion is saved in copy_S_matrix
 copy_S_matrix = S_matrix
 call dgetrf(dim, dim, copy_S_matrix, dim, ipiv, info)
!if (info == 0) then
!   print *, "S is invertible"
!else
!   print *, "S is not invertible"
!end if
   
 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Fill S matrix + inversion - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Fill S matrix + inversion - WALLCLOCK TIME =", omp_time_2 - omp_time_1

 !###############################################################
 print*, "Fill AB matrix"
 call cpu_time (cpu_time_1)
 omp_time_1 = omp_get_wtime()

 do det_I = 1, N_det 
  do det_J = 1, N_det
    AB_matrix (det_I, det_J) = A_IJ(det_I, det_J, 1)
    AB_matrix (det_I + N_det, det_J + N_det) = AB_matrix (det_I, det_J)
     
    AB_matrix (det_I, det_J + N_det) = B_IJ(det_I, det_J, 1) 
    AB_matrix (det_I + N_det, det_J) = AB_matrix (det_I, det_J + N_det)
    
  enddo
 enddo
 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Fill AB matrix - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Fill AB matrix - WALLCLOCK TIME =", omp_time_2 - omp_time_1

 !###############################################################
 print*, "Multiply S^-1 matrix AB matrix"
 call cpu_time (cpu_time_1)
 omp_time_1 = omp_get_wtime()

 copy_AB_matrix = AB_matrix
 call dsymm("L", "U", dim, dim, 1.d0, copy_S_matrix, dim, copy_AB_matrix, dim,&
 & 0.d0, S_inv_time_AB, dim)

 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Multiply S^-1 matrix AB matrix - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Multiply S^-1 matrix AB matrix - WALLCLOCK TIME =", omp_time_2 - omp_time_1

 !###############################################################
 print*, "Eigenvalue problem"
 call cpu_time (cpu_time_1)
 omp_time_1 = omp_get_wtime()
 
 call lapack_diagd(alphar,eigenvectors,S_inv_time_AB,dim,dim)

 call cpu_time (cpu_time_2)
 omp_time_2 = OMP_get_wtime()
 print*, "END Eigenvalue problem - CPUTIME =", cpu_time_2 - cpu_time_1
 print*, "END Eigenvalue problem - WALLCLOCK TIME =", omp_time_2 - omp_time_1

call dlasrt("I", dim, alphar,info)
 do det_I=1, dim
  print*, alphar(det_I)
 enddo

 deallocate(work)
end subroutine

