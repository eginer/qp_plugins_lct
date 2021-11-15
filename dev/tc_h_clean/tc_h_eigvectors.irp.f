 BEGIN_PROVIDER [double precision, eigval_right_tc, (N_states)]
&BEGIN_PROVIDER [double precision, eigval_left_tc, (N_states)]
&BEGIN_PROVIDER [double precision, reigvec_tc, (N_det,N_states)]
&BEGIN_PROVIDER [double precision, leigvec_tc, (N_det,N_states)]
  implicit none
  integer :: i,idx_dress,j
  logical :: converged, dagger
 BEGIN_DOC
 ! eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian 
 END_DOC
 if(full_tc_h_solver.or.N_det.le.N_det_max_full)then
   double precision, allocatable :: reigvec_tc_tmp(:,:),leigvec_tc_tmp(:,:),eigval_right_tmp(:)
   allocate(reigvec_tc_tmp(N_det,N_det),leigvec_tc_tmp(N_det,N_det),eigval_right_tmp(N_det))
   integer :: n_real_tc_eigval_right
   call non_hrmt_real_diag(N_det,htilde_matrix_elmt,reigvec_tc_tmp,leigvec_tc_tmp,n_real_tc_eigval_right,eigval_right_tmp)
   do i = 1, N_states
    eigval_right_tc(i) = eigval_right_tmp(i)
    do j = 1, N_det
     reigvec_tc(j,i) = reigvec_tc_tmp(j,i)
     leigvec_tc(j,i) = leigvec_tc_tmp(j,i)
    enddo
   enddo
 else
   allocate(reigvec_tc_tmp(N_det,N_states_diag))
   idx_dress = 1
   dagger = .False. ! to compute the RIGHT eigenvector 
   reigvec_tc_tmp = 0.d0
   do i = 1, N_states
    do j = 1, N_det
     reigvec_tc_tmp(j,i) = psi_coef(j,i)
    enddo
   enddo
   call iterative_davidson_tc(psi_det,reigvec_tc_tmp,N_det,N_states,N_states_diag,idx_dress,dagger, & 
                              eigval_right_tc,converged)
   do i = 1, N_states
    do j = 1, N_det
     reigvec_tc(j,i) = reigvec_tc_tmp(j,i) 
    enddo
   enddo
   if(.not.converged)then
    print*,'Stopping ... Davidson did not converge ...'
    stop
   endif
   if(comp_left_eigv)then
     idx_dress = 1
     dagger = .True. ! to compute the LEFT  eigenvector 
     reigvec_tc_tmp = 0.d0
     do i = 1, N_states
      do j = 1, N_det
       reigvec_tc_tmp(j,i) = psi_coef(j,i)
      enddo
     enddo
     call iterative_davidson_tc(psi_det,reigvec_tc_tmp,N_det,N_states,N_states_diag,idx_dress,dagger, & 
                                eigval_left_tc,converged)
     do i = 1, N_states
      do j = 1, N_det
       leigvec_tc(j,i) = reigvec_tc_tmp(j,i) 
      enddo
     enddo
   endif
 endif
 print*,'eigval_right_tc = ',eigval_right_tc(1)
 print*,'N_det tc        = ',N_det
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

subroutine iterative_davidson_tc(psidet,psicoef,ndet,N_st,N_st_diag,idx_dress,dagger,energies_out,converged)
 implicit none
 BEGIN_DOC
! iterative Davidson technique to obtain the lowests states of the TC Hamiltonian
!
! psidet      : list of Slater determinants in input
!             
! psicoef     : guess of the coefficients 
!             
! ndet         : number of determinants in input
!             
! N_st        : Number of state targetted 
!             
! N_st_diag   : Maximum number of davidson states to diagonalize
!             
! idx_dress   : index of the Slater determinant to realize the vector dressing 
!             
! dagger      : if .True., then computes the eigenstate of the (Htilde)^DAGGER operator 
!
! energies_out: lowest N_st energies 
!
! converged   : if .False. it means that it did not converge 
 END_DOC

 use bitmasks
 double precision, intent(inout)   :: psicoef(ndet,N_st_diag),energies_out(N_st)
 integer(bit_kind), intent(in)     :: psidet(N_int,2,ndet)
 integer, intent(in)               :: ndet,N_st,idx_dress,N_st_diag
 logical, intent(in)            :: dagger
 logical, intent(out)              :: converged
 integer :: istate,j,i
 double precision, allocatable :: H_jj(:),Dress_jj(:),Dressing_vec(:,:),energies(:)
 double precision, allocatable :: htilde_psi(:)
 double precision :: residual,u_dot_v,e_expect,e_before, delta_e
 external H_u_0_nstates_openmp

 integer :: n_iter_max 
 n_iter_max = 1000
 allocate(H_jj(ndet), Dress_jj(ndet), Dressing_vec(ndet,N_st))
 allocate(energies(N_st_diag),htilde_psi(ndet))
 !!! Create a Guess 
 do istate = N_st + 1, N_st_diag
  psicoef(istate,istate) = 1.d0
 enddo
 !!! H matrix diagonal elements and nul diagonal dressing 
 do j = 1, N_det
  H_jj(j) = H_matrix_diag_all_dets(j)
  Dress_jj(j) = 0.d0
 enddo

 residual = 1.d0
 j = 0
 print*,'Iterative Meta Davidson diagonalization with vector dressing '
 print*,'************************************************************'
 if(dagger)then
  print*,'Computing the LEFT eigenvector '
 else 
  print*,'Computing the RIGHT eigenvector '
 endif
 print*,''
 e_before = 1.d0 
 delta_e = 1.d0
 converged = .False.
 do while (.not.converged)
  j += 1
  print*,'iteration = ',j
  !!! Compute the dressing vector 
  call get_dressing_tc_for_dav(psicoef(1,1),psidet,ndet,N_st,dagger,Dressing_vec)
  !!! Call the Davidson for dressing vector 
  call dav_double_dressed(psicoef,H_jj,Dress_jj,Dressing_vec,idx_dress,energies,ndet,N_st,N_st_diag,converged,H_u_0_nstates_openmp)

  if(.not.dagger)then
   !!! Compute htilde_psi = Htilde | psicoef >
   call get_htilde_psi(psidet,psicoef,ndet,htilde_psi)
  else
   call get_htilde_dagger_psi(psidet,psicoef,ndet,htilde_psi)
  endif
  !!! Compute the expectation value < psicoef | Htilde | psicoef >
  e_expect = u_dot_v(htilde_psi,psicoef(1,1),ndet)
  if(j>1)then
   delta_e = energies(1) - e_before 
  endif
  e_before = energies(1)
  !!! Compute the residual < psicoef | (Htilde - E) | psicoef >
  htilde_psi(1:ndet) += -energies(1) * psicoef(1:ndet,1)
  residual = u_dot_v(psicoef(1,1),htilde_psi,ndet)
  residual = dabs(residual)
  print*, ' energies = ',energies(1)
  if(j>1)then
   print*,'  Delta E = ', delta_e
  endif
  print*, ' residual = ',residual
  converged = residual.lt.threshold_davidson.or.dabs(delta_e).lt.1.d-5
  if(j>n_iter_max)then
   converged = .False.
   exit
  endif
 enddo
 if(.not.converged)then
  print*,'iterative_davidson_tc did not converge after ',n_iter_max,'iterations ...'
 endif
 print*,'Converged TC energy = ',energies(1)
 do i = 1, N_st
  energies_out(i) = energies(i)
 enddo
end

 BEGIN_PROVIDER [ double precision, h_mono_comp_right_tc]
&BEGIN_PROVIDER [ double precision, h_eff_comp_right_tc]
&BEGIN_PROVIDER [ double precision, h_deriv_comp_right_tc]
&BEGIN_PROVIDER [ double precision, h_three_comp_right_tc]
&BEGIN_PROVIDER [ double precision, h_tot_comp_right_tc]
 implicit none
 use bitmasks
 double precision, allocatable :: psicoef(:)
 integer(bit_kind),allocatable :: psidet(:,:,:)
 allocate(psicoef(N_det),psidet(N_int,2,N_det))
 integer :: i
 do i = 1, N_det
  psidet(1:N_int,1:2,i) = psi_det(1:N_int,1:2,i)
  psicoef(i) = psi_coef(i,1)
 enddo
 call get_e_components_htilde(psidet,psicoef,n_det,h_mono_comp_right_tc,h_eff_comp_right_tc,& 
                              h_deriv_comp_right_tc,h_three_comp_right_tc,h_tot_comp_right_tc)

END_PROVIDER 
