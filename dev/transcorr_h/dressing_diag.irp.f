
program test_dressing
 implicit none
 read_wf = .True.
 touch read_wf
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 170
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 extra_grid_type_sgn = 1 
 touch extra_grid_type_sgn 
 my_extra_grid_becke = .False.
 touch my_extra_grid_becke 
 print*,'Warning : the Becke grid parameters are automatically set to '
 print*,'my_n_pt_a_grid = ',my_n_pt_a_grid
 print*,'my_n_pt_r_grid = ',my_n_pt_r_grid
 print*,'If you want to modify them, you have to modify the following file '
 print*,'qp2/plugins/qp_plugins_lct/dev/transcorr_h/transcorr_general.irp.f'
 print*,'and recompile doing ninja'
 if(linear_tc)then
  three_body_h_tc = .False. 
  touch three_body_h_tc
  grad_squared = .False. 
  touch grad_squared 
 endif
 if(read_tc_ints)then
  call read_fcidump_1_tc
 endif

 call test_dressing_diag
end

subroutine get_dressed_matrix(u0,h_dressed,idress)
 use bitmasks
 BEGIN_DOC
! You enter with u0, a good guess to the right eigenvector of the TC Hamiltonian
!
! You get out with a dressed symmetric matrix taking the effect of (Htilde - H) |u0>
 END_DOC
 implicit none
 integer, intent(in) :: idress
 double precision, intent(in) :: u0(N_det)
 double precision, intent(out):: h_dressed(N_det,N_det)
 double precision, allocatable :: delta_u0(:), delta_mat(:,:)
 double precision :: a
 integer :: i
 a = 1.d0
 allocate(delta_u0(N_det),delta_mat(N_det,N_det))
 delta_mat = htilde_matrix_elmt - h_matrix_all_dets ! Delta = Htilde - H
 !!!!!!!!!!!!! Computing the dressing vector 
 delta_u0 = 0.d0
 call h_non_hermite(delta_u0,u0(idress),delta_mat,a,1,N_det)  ! delta_u0 = Delta |u0> 
 delta_u0 *= 1.d0/u0(idress)
 !!!!!!!!!!!!! Computing the dressing matrix 
 h_dressed = h_matrix_all_dets
 h_dressed(idress,idress) += delta_u0(idress) 
 
 do i = 2, N_det
  h_dressed(idress,idress) -= delta_u0(i)/u0(idress) * u0(i)
  h_dressed(idress,i) += delta_u0(i)
  h_dressed(i,idress) += delta_u0(i)
 enddo
end


subroutine test_dressing_diag
 implicit none
 double precision :: e0
 double precision, allocatable :: u0(:,:),u1(:), h_dressed(:,:)
 double precision, allocatable :: eigvalues(:),eigvectors(:,:)
 double precision :: a,res,tmp1,tmp2,ebefore,thr
 integer :: N_st_diag
 double precision, allocatable :: H_jj(:)
 logical :: converged
 integer :: i,j,nitermax,idress
 thr = threshold_davidson

 !!! You assume that the first determinant is dominant 
 idress = 1

 a = 1.d0

 allocate(u1(N_det), h_dressed(N_det, N_det) )
 nitermax = 10
 print*,'Starting with a guess from the usual H'
 !! Create a guess
  
  print*,''
  print*,'********'
  res = 1.d+10
  j = 0
  ebefore = 0.d0

 if(full_tc_h_solver)then
  allocate(u0(N_det,1)) 
  do i = 1, N_det
!   u0(i,1) = psi_coef(i,1)
   u0(i,1) = reigvec_trans(i,1)
  enddo
  allocate(eigvalues(N_det), eigvectors(N_det,N_det))
  do while(res .gt. thr) 
   j += 1
   print*,''
   print*,'Iteration j ',j
   call get_dressed_matrix(u0,h_dressed,idress)
   
   call lapack_diagd(eigvalues,eigvectors,h_dressed,N_det,N_det)
   e0 = eigvalues(1)
   if(j.gt.1)then
    print*,'e0 = ',e0,dabs(ebefore-e0)
   else
    print*,'e0 = ',e0
   endif
   u0(:,1) = eigvectors(:,1) 
   u1 = 0.d0
   call h_non_hermite(u1,u0(1,1),htilde_matrix_elmt,a,1,N_det)   
   u1(:) -= e0 * u0(:,1)
   res = 0.d0
   do i = 1, N_det
    res += u1(i)*u1(i)
   enddo
   res = dsqrt(res)
   print*,'Norm of the residual vector ', res 
   ebefore = e0
  enddo
  print*,''
  print*,'Comparison between eigenvalues'
  print*,'*****'
  i = 1
  print*,'Ground state '
  print*,eigvalues(i),eigval_trans(i),dabs(eigvalues(i)-eigval_trans(i))
  print*,'*****'
  print*,'Excited states '
  do i = 2, N_det
   write(*,'(I3,X,3(F16.10,X))')i,eigvalues(i),eigval_trans(i),dabs(eigvalues(i)-eigval_trans(i))
  enddo
  print*,''
 else
  allocate(H_jj(N_det))
  do i = 1, N_det
   H_jj(i) = h_matrix_all_dets(i,i)
  enddo
  N_st_diag = N_states_diag
  allocate(u0(N_det,N_states_diag)) 
  u0 = 0.d0
  do i = 1, N_det
   u0(i,1) = psi_coef(i,1)
!   u0(i,1) = reigvec_trans(i,1)
  enddo
  j = 0
  do while(res .gt. thr) 
   j += 1
   print*,''
   print*,'Iteration j ',j
   call get_dressed_matrix(u0(1,1),h_dressed,idress)
   
   call davidson_general(u0,H_jj,e0,N_det,N_det,1,N_st_diag,converged,h_dressed)
 
   if(j.gt.1)then
    print*,'e0 = ',e0,dabs(ebefore-e0)
   else
    print*,'e0 = ',e0
   endif
   u1 = 0.d0
   call h_non_hermite(u1,u0(1,1),htilde_matrix_elmt,a,1,N_det)   
   u1(:) -= e0 * u0(:,1)
   res = 0.d0
   do i = 1, N_det
    res += u1(i)*u1(i)
   enddo
   res = dsqrt(res)
   print*,'Norm of the residual vector ', res 
   ebefore = e0
  enddo
  print*,'Comparison between eigenvalues'
  print*,'*****'
  print*,'Ground state '
  print*,e0,eigval_trans(1),dabs(e0-eigval_trans(1))
 endif
 print*,''
 print*,'End of iterations '

 double precision :: scal
 print*,'Comparision between Ground state Eigenvectors '
 print*,'Det,       Iterative,        Exact'
 tmp1 = u0(1,1)/dabs(u0(1,1))
 tmp2 = reigvec_trans(1,1)/dabs(reigvec_trans(1,1))
 scal = 0.d0
 do i = 1, N_det
  scal += u0(i,1) * reigvec_trans(i,1)
  write(*,'(I4,X,3(F16.10,X))')i,u0(i,1)/tmp1,reigvec_trans(i,1)/tmp2
 enddo
 scal = dsqrt(dabs(scal))
 print*,''
 print*,'Scalar product between the two vectors '
 print*,scal
 print*,''
! print*,'Residual vector'
! do i = 1, N_det
!  res += u1(i)*u1(i)
!  print*,u1(i)
! enddo
end
