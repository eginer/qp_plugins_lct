program diag_dress_iter
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
  read_wf = .True.
  touch read_wf 
  call routine

end


subroutine routine
 implicit none
 integer :: istate,j
 double precision, allocatable :: u_in(:,:),H_jj(:),Dress_jj(:),Dressing_vec(:,:),energies(:)
 double precision, allocatable :: htilde_psi(:)
 double precision :: residual,u_dot_v
 integer :: idx_dress,sze,N_st,N_st_diag
 logical :: converged
 external hcalc_template
 sze = N_det
 idx_dress = 1
 N_st = 1
 N_st_diag = N_states_diag
 allocate(H_jj(sze), Dress_jj(sze), Dressing_vec(sze,N_st))
 allocate(u_in(sze,N_st_diag), energies(N_st_diag))
 allocate(htilde_psi(sze))
 !!! Create a Guess 
 u_in = 0.d0
 do istate = 1, N_st
  do j = 1, sze
   u_in(j,istate) = psi_coef(j,istate)
  enddo
 enddo
 do istate = N_st + 1, N_st_diag
  u_in(istate,istate) = 1.d0
 enddo
 !!! H matrix diagonal elements and nul diagonal dressing 
 do j = 1, N_det
!  H_jj(i) = diag_htilde(i)
  H_jj(j) = H_matrix_all_dets(j,j)
  Dress_jj(j) = 0.d0
 enddo

 residual = 1.d0
 do while (residual.gt.threshold_davidson)
  call set_dress_vec_s(u_in,psi_det,sze,N_st,Dressing_vec)
  call dav_double_dressed(u_in,H_jj,Dress_jj,Dressing_vec,idx_dress,energies,sze,N_st,N_st_diag,converged,hcalc_template)
  call htilde_psi_no_store_no_provide(psi_det,u_in,sze,htilde_psi)
  htilde_psi(1:sze) += -energies(1) * u_in(1:sze,1)
  residual = u_dot_v(u_in(1,1),htilde_psi,sze)
  residual = dabs(residual)
  print*,'energies = ',energies
  print*,'residual = ',residual
 enddo

end
