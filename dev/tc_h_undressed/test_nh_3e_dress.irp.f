program test_nh_3e_dress
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

subroutine test_dressing_diag
 implicit none
 double precision, allocatable :: u0(:),h_dressed(:,:)
 double precision, allocatable :: eigvalues(:),eigvectors(:,:)
 double precision :: e0
 integer :: idress,i
 idress = 1
 allocate(u0(N_det),h_dressed(N_det,N_det))
 allocate(eigvalues(N_det), eigvectors(N_det,N_det))
 u0 = 0.d0
 u0(1:N_det) = psi_coef(1:N_det,1)
 call get_dressed_matrix_nh_3e(u0,h_dressed,idress)

   call lapack_diagd(eigvalues,eigvectors,h_dressed,N_det,N_det)
   do i = 1, N_det
    print*,i,eigvalues(i)
   enddo
   e0 = eigvalues(1)
!   if(j.gt.1)then
!    print*,'e0 = ',e0,dabs(ebefore-e0)
!   else
    print*,'e0 = ',e0
!   endif
  i = 1
  print*,'Ground state '
  print*,eigvalues(i),eigval_trans(i),dabs(eigvalues(i)-eigval_trans(i))
  print*,'*****'

end

subroutine get_dressed_matrix_nh_3e(u0,h_dressed,idress)
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
! delta_mat = nh_3e_matrix_elmt! Delta = Htilde - H
! delta_mat = htilde_matrix_elmt - h_matrix_all_dets ! Delta = Htilde - H
 !!!!!!!!!!!!! Computing the dressing vector 
 delta_u0 = 0.d0
 call h_non_hermite(delta_u0,u0(idress),delta_mat,a,1,N_det)  ! delta_u0 = Delta |u0> 
 do i = 1, N_det
  print*,i,delta_u0(i)
 enddo

 delta_u0 *= 1.d0/u0(idress)
 !!!!!!!!!!!!! Computing the dressing matrix 
! h_dressed =  htilde_matrix_elmt - nh_3e_matrix_elmt
! print*,'Htilde(i,i) = ',htilde_matrix_elmt(N_det,N_det)
! print*,'K+L   (i,i) = ',nh_3e_matrix_elmt(N_det,N_det)
! h_dressed = h_matrix_all_dets

 h_dressed(idress,idress) += delta_u0(idress) 
 
 do i = 2, N_det
  h_dressed(idress,idress) -= delta_u0(i)/u0(idress) * u0(i)
  h_dressed(idress,i) += delta_u0(i)
  h_dressed(i,idress) += delta_u0(i)
 enddo
! print*,'Matrix'
! do i = 1, N_det
!  write(*,'(100(F10.6,X))'),h_dressed(i,:)
! enddo
end

