 BEGIN_PROVIDER [double precision, eigval_tc, (N_states)]
&BEGIN_PROVIDER [double precision, reigvec_tc, (N_det,N_states)]
&BEGIN_PROVIDER [double precision, leigvec_tc, (N_det,N_states)]
  implicit none
  integer :: i,idx_dress,j
  logical :: converged
 BEGIN_DOC
 ! eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian 
 END_DOC
 if(full_tc_h_solver)then
  double precision, allocatable :: reigvec_tc_tmp(:,:),leigvec_tc_tmp(:,:),eigval_tmp(:)
  allocate(reigvec_tc_tmp(N_det,N_det),leigvec_tc_tmp(N_det,N_det),eigval_tmp(N_det))
  integer :: n_real_tc_eigval
  call non_hrmt_real_diag(N_det,htilde_matrix_elmt,reigvec_tc_tmp,leigvec_tc_tmp,n_real_tc_eigval,eigval_tmp)
  do i = 1, N_states
   eigval_tc(i) = eigval_tmp(i)
   do j = 1, N_det
    reigvec_tc(j,i) = reigvec_tc_tmp(j,i)
    leigvec_tc(j,i) = leigvec_tc_tmp(j,i)
   enddo
  enddo
 else
  idx_dress = 1
  allocate(reigvec_tc_tmp(N_det,N_states_diag))
  reigvec_tc_tmp = 0.d0
  do i = 1, N_states
   do j = 1, N_det
    reigvec_tc_tmp(j,i) = psi_coef(j,i)
   enddo
  enddo
  call iterative_davidson_tc(psi_det,reigvec_tc_tmp,N_det,N_states,N_states_diag,idx_dress,eigval_tc,converged)
  do i = 1, N_states
   do j = 1, N_det
    reigvec_tc(j,i) = reigvec_tc_tmp(j,i) 
   enddo
  enddo
 endif
END_PROVIDER 

subroutine write_left_right
 implicit none
 double precision, allocatable :: reigvec_tc_tmp(:,:),leigvec_tc_tmp(:,:)
 allocate(leigvec_tc_tmp(N_det,N_states),reigvec_tc_tmp(N_det,N_states))
 integer :: i,j
 do i = 1, N_states
  do j = 1, N_det
   leigvec_tc_tmp(j,i) = leigvec_tc(j,i)
   reigvec_tc_tmp(j,i) = reigvec_tc(j,i)
  enddo
 enddo
 call ezfio_set_tc_h_clean_reigvec_tc(reigvec_tc_tmp)
 call ezfio_set_tc_h_clean_leigvec_tc(leigvec_tc_tmp)
end

subroutine iterative_davidson_tc(psidet,psicoef,sze,N_st,N_st_diag,idx_dress,energies_out,converged)
 implicit none
 double precision, intent(inout)   :: psicoef(sze,N_st_diag),energies_out(N_st)
 integer(bit_kind), intent(in)     :: psidet(N_int,2,sze)
 integer, intent(in)               :: sze,N_st,idx_dress,N_st_diag
 logical, intent(out)              :: converged
 use bitmasks
 integer :: istate,j,i
 double precision, allocatable :: H_jj(:),Dress_jj(:),Dressing_vec(:,:),energies(:)
 double precision, allocatable :: htilde_psi(:)
 double precision :: residual,u_dot_v,e_expect,e_before, delta_e
 external hcalc_template
 allocate(H_jj(sze), Dress_jj(sze), Dressing_vec(sze,N_st))
 allocate(energies(N_st_diag),htilde_psi(sze))
 !!! Create a Guess 
 do istate = N_st + 1, N_st_diag
  psicoef(istate,istate) = 1.d0
 enddo
 !!! H matrix diagonal elements and nul diagonal dressing 
 do j = 1, N_det
  H_jj(j) = H_matrix_all_dets(j,j)
  Dress_jj(j) = 0.d0
 enddo

 residual = 1.d0
 j = 0
 print*,'Iterative Meta Davidson diagonalization with vector dressing '
 e_before = 1.d0 
 delta_e = 1.d0
 do while (residual.gt.threshold_davidson.or.dabs(delta_e).gt.1.d-6)
  j += 1
  print*,'iteration = ',j
  !!! Compute the dressing vector 
  call set_dress_vec_s(psicoef(1,1),psidet,sze,N_st,Dressing_vec)
  !!! Call the Davidson for dressing vector 
  call dav_double_dressed(psicoef,H_jj,Dress_jj,Dressing_vec,idx_dress,energies,sze,N_st,N_st_diag,converged,hcalc_template)

  !!! Compute htilde_psi = Htilde | psicoef >
  call htilde_psi_no_store_no_provide(psidet,psicoef,sze,htilde_psi)
  !!! Compute the expectation value < psicoef | Htilde | psicoef >
  e_expect = u_dot_v(htilde_psi,psicoef(1,1),sze)
  if(j>1)then
   delta_e = energies(1) - e_before 
  endif
  e_before = energies(1)
  !!! Compute the residual < psicoef | (Htilde - E) | psicoef >
  htilde_psi(1:sze) += -energies(1) * psicoef(1:sze,1)
  residual = u_dot_v(psicoef(1,1),htilde_psi,sze)
  residual = dabs(residual)
  print*, ' energies = ',energies
  if(j>1)then
   print*,'  Delta E = ', delta_e
  endif
  print*, ' residual = ',residual
 enddo
 print*,'Converged TC energy = ',energies(1)
 do i = 1, N_st
  energies_out(i) = energies(i)
 enddo
end
