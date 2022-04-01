  use bitmasks ! you need to include the bitmasks_module.f90 features
 BEGIN_PROVIDER [ double precision, singles_hf_mat_elem]
 implicit none
 integer :: i,degree
 double precision :: hmono,heff,hderiv,hthree,htilde_ij
 singles_hf_mat_elem = 0.d0
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree == 1)then
   call htilde_mu_mat(psi_det(1,1,i),HF_bitmask,N_int,hmono,heff,hderiv,hthree,htilde_ij)
   singles_hf_mat_elem += dabs(htilde_ij)
  endif
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, h_tilde_expect_right]
&BEGIN_PROVIDER [ double precision, h_tilde_expect_psi_coef]
&BEGIN_PROVIDER [ double precision, h_tilde_dagger_expect_right]
&BEGIN_PROVIDER [ double precision, h_tilde_dagger_expect_psi_coef]
  implicit none
  integer :: i,j
  double precision :: hmono,heff,hderiv,hthree,htilde_ij
  h_tilde_expect_right = 0.d0
  h_tilde_dagger_expect_right = 0.d0
  h_tilde_expect_psi_coef = 0.d0
  h_tilde_dagger_expect_psi_coef = 0.d0
  do i = 1, N_det
   do j = 1, N_det
    call htilde_mu_mat(psi_det(1,1,i),psi_det(1,1,j),N_int,hmono,heff,hderiv,hthree,htilde_ij)
    h_tilde_expect_right += htilde_ij * reigvec_tc(i,1) * reigvec_tc(j,1)
    h_tilde_expect_psi_coef += htilde_ij * psi_coef(i,1) * psi_coef(j,1)
    call htilde_mu_mat(psi_det(1,1,j),psi_det(1,1,i),N_int,hmono,heff,hderiv,hthree,htilde_ij)
    h_tilde_dagger_expect_right += htilde_ij * reigvec_tc(i,1) * reigvec_tc(j,1)
    h_tilde_dagger_expect_psi_coef += htilde_ij * psi_coef(i,1) * psi_coef(j,1)
   enddo
  enddo
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, e_from_right_eigv]
 implicit none
 integer :: i,degree
 double precision :: hmono,heff,hderiv,hthree,htilde_ij
 e_from_right_eigv = 0.d0
 do i = 1, N_det
  call htilde_mu_mat(HF_bitmask,psi_det(1,1,i),N_int,hmono,heff,hderiv,hthree,htilde_ij)
  e_from_right_eigv += htilde_ij * reigvec_tc(i,1)/reigvec_tc(1,1)
  call htilde_mu_mat(psi_det(1,1,i),HF_bitmask,N_int,hmono,heff,hderiv,hthree,htilde_ij)
 enddo

 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_pt2_tc]
&BEGIN_PROVIDER [ double precision, e_pt2_tc_single]
&BEGIN_PROVIDER [ double precision, e_pt2_tc_double]
 implicit none 
 integer :: i,degree
 double precision :: hmono,heff,hderiv,hthree,htilde_ij,coef_pt1,e_i0,delta_e
 e_pt2_tc = 0.d0
 e_pt2_tc_single = 0.d0
 e_pt2_tc_double = 0.d0
 do i = 1, N_det
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree == 1 .or. degree == 2)then
   call htilde_mu_mat(psi_det(1,1,i),HF_bitmask,N_int,hmono,heff,hderiv,hthree,htilde_ij)
   call htilde_mu_mat(psi_det(1,1,i),psi_det(1,1,i),N_int,hmono,heff,hderiv,hthree,e_i0)
   delta_e = e_tilde_00 - e_i0
   coef_pt1 = htilde_ij / delta_e
   call htilde_mu_mat(HF_bitmask,psi_det(1,1,i),N_int,hmono,heff,hderiv,hthree,htilde_ij)
   e_pt2_tc += coef_pt1 * htilde_ij
   if(degree == 1)then
    e_pt2_tc_single += coef_pt1 * htilde_ij
   else 
    e_pt2_tc_double += coef_pt1 * htilde_ij
   endif
  endif
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e_tilde_00]
 implicit none
 double precision :: hmono,heff,hderiv,hthree,htilde_ij
 call htilde_mu_mat(HF_bitmask,HF_bitmask,N_int,hmono,heff,hderiv,hthree,e_tilde_00)
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, norm_ground_left_right]
 implicit none
 integer :: i
 norm_ground_left_right = 0.d0
 do i = 1, N_det
  norm_ground_left_right += reigvec_tc(i,1) * leigvec_tc(i,1)
 enddo

 if(dabs(norm_ground_left_right).lt.0.01d0)then
  print *,'Two small norm !!'
  do i = 1,N_det
   print *,i,reigvec_tc(i,1) , leigvec_tc(i,1)
  enddo
!  call routine_save_right
!  stop
 endif

 END_PROVIDER

 BEGIN_PROVIDER [ double precision, norm_ground_right]
 implicit none
 integer :: i
 norm_ground_right = 0.d0
 do i = 1, N_det
  norm_ground_right += reigvec_tc(i,1) * reigvec_tc(i,1)
 enddo
 END_PROVIDER

 BEGIN_PROVIDER [ double precision, norm_ground_left]
 implicit none
 integer :: i
 norm_ground_left = 0.d0
 do i = 1, N_det
  norm_ground_left += leigvec_tc(i,1) * leigvec_tc(i,1)
 enddo
 if(dabs(norm_ground_left).lt.0.1d0)then
  print *,'Warning ! norm_ground_left too small!'
  print *,'norm_ground_left = ',norm_ground_left
 endif
 END_PROVIDER

 BEGIN_PROVIDER [ integer, index_HF_psi_det]
 implicit none
 integer :: i,degree
 do i = 1, N_det
   call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
   if(degree == 0)then
    index_HF_psi_det = i
    exit
   endif
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, eigval_right_tc, (N_states)]
&BEGIN_PROVIDER [double precision, eigval_left_tc, (N_states)]
&BEGIN_PROVIDER [double precision, reigvec_tc, (N_det,N_states)]
&BEGIN_PROVIDER [double precision, leigvec_tc, (N_det,N_states)]

  BEGIN_DOC
  ! eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian 
  END_DOC

  implicit none
  integer                       :: i, idx_dress, j
  logical                       :: converged, dagger
  integer                       :: n_real_tc_eigval_right,igood_r,igood_l
  double precision, allocatable :: reigvec_tc_tmp(:,:),leigvec_tc_tmp(:,:),eigval_right_tmp(:)

  PROVIDE N_det N_int


  if(full_tc_h_solver.or.N_det.le.N_det_max_full)then

    allocate(reigvec_tc_tmp(N_det,N_det),leigvec_tc_tmp(N_det,N_det),eigval_right_tmp(N_det))

    call non_hrmt_real_diag(N_det,htilde_matrix_elmt,reigvec_tc_tmp,leigvec_tc_tmp,n_real_tc_eigval_right,eigval_right_tmp)
    double precision, allocatable :: coef_hf_r(:),coef_hf_l(:)
    integer, allocatable :: iorder(:)
    allocate(coef_hf_r(N_det),coef_hf_l(N_det),iorder(N_det))
    do i = 1,N_det
     iorder(i) = i
     coef_hf_r(i) = -dabs(reigvec_tc_tmp(index_HF_psi_det,i))
    enddo
    call dsort(coef_hf_r,iorder,N_det)
    igood_r = iorder(1)
    do i = 1,N_det
     iorder(i) = i
     coef_hf_l(i) = -dabs(leigvec_tc_tmp(index_HF_psi_det,i))
    enddo
    call dsort(coef_hf_l,iorder,N_det)
    igood_l = iorder(1)

    if(igood_r.ne.igood_l.and.igood_r.ne.1)then
     print *,''
     print *,'Warning, the left and right eigenvectors are "not the same" '
     print *,'Warning, the ground state is not dominated by HF...'
     print *,'State with largest RIGHT coefficient of HF ',igood_r
     print *,'coef of HF in RIGHT eigenvector = ',reigvec_tc_tmp(index_HF_psi_det,igood_r)
     print *,'State with largest LEFT  coefficient of HF ',igood_l
     print *,'coef of HF in LEFT  eigenvector = ',leigvec_tc_tmp(index_HF_psi_det,igood_l)
    endif
    if(state_following)then
     print *,'Following the states with the largest coef on HF'
     print *,'igood_r,igood_l',igood_r,igood_l
     i= igood_r
     eigval_right_tc(1) = eigval_right_tmp(i)
     do j = 1, N_det
       reigvec_tc(j,1) = reigvec_tc_tmp(j,i)
     enddo
     i= igood_l
     eigval_left_tc(1)  = eigval_right_tmp(i)
     do j = 1, N_det
       leigvec_tc(j,1) = leigvec_tc_tmp(j,i)
     enddo
    else 
     do i = 1, N_states
       eigval_right_tc(i) = eigval_right_tmp(i)
       eigval_left_tc(i)  = eigval_right_tmp(i)
       do j = 1, N_det
         reigvec_tc(j,i) = reigvec_tc_tmp(j,i)
         leigvec_tc(j,i) = leigvec_tc_tmp(j,i)
       enddo
     enddo
    endif

  else

    allocate( reigvec_tc_tmp(N_det,N_states_diag) )

    idx_dress = 1
    dagger = .False. ! to compute the RIGHT eigenvector 
    reigvec_tc_tmp = 0.d0
    do i = 1, N_states
      do j = 1, N_det
        reigvec_tc_tmp(j,i) = psi_coef(j,i)
      enddo
    enddo

    call iterative_davidson_tc(psi_det, reigvec_tc_tmp, N_det, N_int, N_states, N_states_diag, idx_dress, dagger, eigval_right_tc, converged)
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
        do j = 1, min(N_det,n_det_max_full)
!        print *,'psi_left_guess(j,i)',psi_left_guess(j,i)
        reigvec_tc_tmp(j,i) = psi_left_guess(j,i)
!         reigvec_tc_tmp(j,i) = reigvec_tc(j,i)
        enddo
      enddo
      call iterative_davidson_tc(psi_det, reigvec_tc_tmp, N_det, N_int, N_states, N_states_diag, idx_dress, dagger, eigval_left_tc, converged)
      do i = 1, N_states
        do j = 1, N_det
          leigvec_tc(j,i) = reigvec_tc_tmp(j,i) 
        enddo
      enddo
    else 
      do i = 1, N_states
        do j = 1, N_det
         leigvec_tc(j,i) = reigvec_tc(j,i) 
        enddo
      enddo
    endif

  endif
  !!!! Normalization of the scalar product of the left/right eigenvectors
  double precision  :: accu, tmp 
  do i = 1, N_states
   !!!! Normalization of right eigenvectors |Phi>
   accu = 0.d0
   do j = 1, N_det
    accu += reigvec_tc(j,i) * reigvec_tc(j,i)
   enddo
   accu = 1.d0/dsqrt(accu)
   do j = 1, N_det
    reigvec_tc(j,i) *= accu 
   enddo
   tmp = reigvec_tc(1,i) / dabs(reigvec_tc(1,i))
   do j = 1, N_det
    reigvec_tc(j,i) *= tmp
   enddo
   !!!! Adaptation of the norm of the left eigenvector such that <chi|Phi> = 1
   accu = 0.d0
   do j = 1, N_det
    accu += leigvec_tc(j,i) * reigvec_tc(j,i)
   enddo
   if(accu.gt.0.d0)then
    accu = 1.d0/dsqrt(accu)
   else
    accu = 1.d0/dsqrt(-accu)
   endif
   tmp = (leigvec_tc(1,i) * reigvec_tc(1,i) )/dabs(leigvec_tc(1,i) * reigvec_tc(1,i))
   do j = 1, N_det
    leigvec_tc(j,i) *= accu * tmp
    reigvec_tc(j,i) *= accu 
   enddo
   print*,'leigvec_tc(1,i),reigvec_tc(1,i) = ',leigvec_tc(1,i),reigvec_tc(1,i)
   accu = 0.d0
   do j = 1, N_det
    accu += leigvec_tc(j,i) * reigvec_tc(j,i)
   enddo
   print*,'norm l/r = ',accu
  enddo

END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigval_right_tc_nucl_rep, (N_states)]
&BEGIN_PROVIDER [double precision, eigval_left_tc_nucl_rep, (N_states)]
 implicit none
integer :: i
 do i = 1,N_states
  eigval_right_tc_nucl_rep(i) += eigval_right_tc(i) + nuclear_repulsion
  eigval_left_tc_nucl_rep(i)  += eigval_left_tc(i) + nuclear_repulsion
 enddo
END_PROVIDER 


subroutine write_left_right()

  implicit none
  integer :: i,j
  double precision, allocatable :: reigvec_tc_tmp(:,:),leigvec_tc_tmp(:,:)

  allocate(leigvec_tc_tmp(N_det,N_states),reigvec_tc_tmp(N_det,N_states))
  do i = 1, N_states
   do j = 1, N_det
    leigvec_tc_tmp(j,i) = leigvec_tc(j,i)
    reigvec_tc_tmp(j,i) = reigvec_tc(j,i)
   enddo
  enddo
  call ezfio_set_tc_h_clean_reigvec_tc(reigvec_tc_tmp)
  call ezfio_set_tc_h_clean_leigvec_tc(leigvec_tc_tmp)
  deallocate( leigvec_tc_tmp , reigvec_tc_tmp )

end subroutine write_left_right



subroutine iterative_davidson_tc(psidet, psicoef, ndet, Nint, N_st, N_st_diag, idx_dress, dagger, energies_out, converged)

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

 implicit none
 integer,           intent(in)   :: Nint, ndet, N_st, idx_dress, N_st_diag
 logical,           intent(in)   :: dagger
 integer(bit_kind), intent(in)   :: psidet(Nint,2,ndet)
 double precision, intent(inout) :: psicoef(ndet,N_st_diag),energies_out(N_st)
 logical, intent(out)            :: converged
 integer                         :: istate, j, i, n_iter_max 
 double precision, allocatable   :: H_jj(:), Dress_jj(:), Dressing_vec(:,:), energies(:)
 double precision, allocatable   :: htilde_psi(:)
 double precision                :: residual,u_dot_v,e_expect,e_before, delta_e
 external H_u_0_nstates_openmp

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
  call get_dressing_tc_for_dav(psicoef(1,1), psidet, ndet, Nint, N_st, dagger, Dressing_vec)
  !!! Call the Davidson for dressing vector 
  call dav_double_dressed(psicoef,H_jj,Dress_jj,Dressing_vec,idx_dress,energies,ndet,N_st,N_st_diag,converged,H_u_0_nstates_openmp)

  if(.not.dagger)then
   !!! Compute htilde_psi = Htilde | psicoef >
   call get_htilde_psi(psidet, psicoef, ndet, Nint, htilde_psi)
  else
   call get_htilde_dagger_psi(psidet, psicoef, ndet, Nint, htilde_psi)
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
  converged = residual.lt.threshold_davidson.or.dabs(delta_e).lt.thresh_it_dav
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

  use bitmasks

  implicit none
  double precision, allocatable :: psicoef(:)
  integer(bit_kind),allocatable :: psidet(:,:,:)
  integer :: i

  PROVIDE N_int N_det

  allocate( psicoef(N_det), psidet(N_int,2,N_det) )
  do i = 1, N_det
    psidet(1:N_int,1:2,i) = psi_det(1:N_int,1:2,i)
    psicoef(i) = psi_coef(i,1)
  enddo

  call get_e_components_htilde(psidet, psicoef, n_det, N_int, h_mono_comp_right_tc, h_eff_comp_right_tc,& 
                              h_deriv_comp_right_tc, h_three_comp_right_tc, h_tot_comp_right_tc)

END_PROVIDER 

BEGIN_PROVIDER [ double precision, psi_left_guess, (n_det_max_full,N_states)]
 implicit none
 print *,'providing psi_left_guess '
 psi_left_guess = 0.d0
 integer :: i
 do i = 1, min(n_det_max_full,N_states)
  psi_left_guess(i,i) = 1.d0
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_r, (N_int,2,N_det) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_r_order, (N_det) ]
 implicit none
 integer :: i
 double precision, allocatable :: psicoef(:)
 allocate(psicoef(N_det))
 do i = 1,N_det
  psicoef(i) = -dsqrt(dabs(reigvec_tc(i,1)*leigvec_tc(i,1)))
  psi_det_sorted_r_order(i) = i
 enddo
 call dsort(psicoef,psi_det_sorted_r_order,N_det)
 do i = 1, N_det
  psi_det_sorted_r(1:N_int,1:2,i) = psi_det(1:N_int,1:2,psi_det_sorted_r_order(i))
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, reigvec_tc_sorted_r, (N_det,N_states)]
&BEGIN_PROVIDER [ double precision, leigvec_tc_sorted_r, (N_det,N_states)]
 implicit none 
 integer :: i,j
 do i = 1, N_det
  reigvec_tc_sorted_r(i,:) = reigvec_tc(psi_det_sorted_r_order(i),:)
  leigvec_tc_sorted_r(i,:) = leigvec_tc(psi_det_sorted_r_order(i),:)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ integer, N_det_thresh_psi_r]
 implicit none 
 integer :: i 
 N_det_thresh_psi_r = 0
 do i = 1, N_det
  if(dabs(reigvec_tc_sorted_r(i,1)).gt.thresh_psi_r)then
   N_det_thresh_psi_r += 1
  endif
 enddo
 print*,'N_det_thresh_psi_r = ',N_det_thresh_psi_r
 END_PROVIDER 

