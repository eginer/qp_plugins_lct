
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
 delta_u0 *= 1.d0/u0(1)
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
 double precision, allocatable :: u0(:),u1(:), h_dressed(:,:)
 double precision, allocatable :: eigvalues(:),eigvectors(:,:)
 double precision :: a,res
 integer :: i,j,nitermax,idress

 !!! You assume that the first determinant is dominant 
 idress = 1


 a = 1.d0
 allocate(u0(N_det), u1(N_det), h_dressed(N_det, N_det))
 allocate(eigvalues(N_det), eigvectors(N_det,N_det))
 nitermax = 10
 !! Create a guess
 do i = 1, N_det
!  u0(i) = psi_coef(i,1)
  u0(i) = reigvec_trans(i,1)
 enddo
  
  print*,''
 do j = 1, nitermax
  print*,''
  print*,'Iteration j ',j
  call get_dressed_matrix(u0,h_dressed,idress)
  call lapack_diagd(eigvalues,eigvectors,h_dressed,N_det,N_det)
  print*,'Temporary eigenvalues '
  do i = 1, N_det 
   print*,'i,',i,eigvalues(i)
  enddo
  e0 = eigvalues(1)
  print*,''
  print*,'e0 = ',e0
  u0(:) = eigvectors(:,1) 
  print*,'New eigenvector '
  do i = 1, N_det
   print*,i,u0(i)
  enddo
  print*,''
  print*,''
  print*,''
  u1 = 0.d0
  call h_non_hermite(u1,u0(1),htilde_matrix_elmt,a,1,N_det)   
  u1 -= e0 * u0
  print*,'Residual vector'
  res = 0.d0
  do i = 1, N_det
   res += u1(i)*u1(i)
   print*,u1(i)
  enddo
  print*,'Norm of the residual vector ', dsqrt(res)
 enddo
  print*,''
end
