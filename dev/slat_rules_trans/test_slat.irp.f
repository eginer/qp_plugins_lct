program pouet
 implicit none
 read_wf = .True.
 touch read_wf
  call routine
  call test_pert
  call test_eigv
 integer :: ith
 ith = 1
! call test_left_right_eigenvalues(ith)
! call test_overlap_matrix
! provide htilde_matrix_elmt
!  provide mo_two_e_eff_dr12_pot_array
! provide ao_two_e_eff_dr12_pot_array
end

subroutine routine
 use bitmasks
 integer(bit_kind) :: key_i(N_int,2), key_j(N_int,2)
 integer :: i,j,degree
 double precision :: hij,s2,hmono,herf,heff,hderiv,htot
 double precision :: accu
 accu = 0.d0
 key_i(:,:) = psi_det(:,:,1)
 call diag_htilde_mat(key_i,hmono,herf,heff,hderiv,htot)
 print*,''
 print*,'H      eigenvalue'
 i = 1
 print*,'E0       = ',CI_energy(i)
 print*,'E0       = ',CI_electronic_energy(i)
 print*,''
 print*,'Htilde eigenvalue'
 i = 1
  print*,'E0_tilde = ',eigval_trans(i) + nuclear_repulsion
end

subroutine test_eigv
 implicit none
 integer :: i,j
 do i = 1, N_det
  print*,reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1)),leigvec_trans(i,1)/dsqrt(leigvec_trans_norm(1)),psi_coef(i,1)
 enddo
 double precision :: accu1,accu2,e
 accu1 = 0.d0
 accu2 = 0.d0
 do i = 1, N_det
  accu2 += reigvec_trans(i,1) * reigvec_trans(i,1)
  do j = 1, N_det
   accu1 += htilde_matrix_elmt(j,i) * reigvec_trans(i,1) * reigvec_trans(j,1)
  enddo
 enddo
 accu2 = dsqrt(accu2)
 e = accu1/accu2
 print*,'right eigenvector'
 print*,'e                = ',e
 accu1 = 0.d0
 accu2 = 0.d0
 do i = 1, N_det
  accu2 += leigvec_trans(i,1) * leigvec_trans(i,1)
  do j = 1, N_det
   accu1 += htilde_matrix_elmt(j,i) * leigvec_trans(i,1) * leigvec_trans(j,1)
  enddo
 enddo
 accu2 = dsqrt(accu2)
 e = accu1/accu2
 print*,'left eigenvector'
 print*,'e                = ',e
 print*,'eigval_trans(1)  = ',eigval_trans(1)
 do i = 1, N_det
  psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
 enddo
 touch psi_coef
 call save_wavefunction
end

subroutine test_pert
 implicit none
 integer :: i,j
 double precision :: accu,hmono,herf,heff,hderiv,htot
 double precision :: accu_mono,accu_double,accu2
 double precision :: pert_mono,pert_double
 double precision :: h00,hii,htotbis,phase,h0i,hi0
 integer          :: degree,exc(0:2,2,2)
 integer          :: h1, p1, h2, p2, s1, s2
 double precision :: norm_t1,norm_t2
 provide reigvec_trans
 accu = 0.d0
 accu2= 0.d0
 accu_mono = 0.d0
 accu_double = 0.d0
 pert_mono = 0.d0
 pert_double = 0.d0
 norm_t1 = 0.d0
 norm_t2 = 0.d0
 call htilde_mat(psi_det(1,1,1),psi_det(1,1,1),hmono,herf,heff,hderiv,h00)
  print*,'Printing the connection '
  do i = 1, N_det
   ! <0|H|i>
   print*,'******************'
   print*,'i = ',i
   call htilde_mat(psi_det(1,1,1),psi_det(1,1,i),hmono,herf,heff,hderiv,h0i)
 
   ! <i|H|i>
   call htilde_mat(psi_det(1,1,i),psi_det(1,1,i),hmono,herf,heff,hderiv,hii)
 
   ! <i|H|0>
   call htilde_mat(psi_det(1,1,i),psi_det(1,1,1),hmono,herf,heff,hderiv,hi0)
 
   call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,1),degree,N_int)
   if(degree==1)then
    call get_single_excitation(psi_det(1,1,i),psi_det(1,1,1),exc,phase,N_int)
    call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
    norm_t1 += (reigvec_trans(i,1)/reigvec_trans(1,1))**2.d0
    print*,'h1,p1',h1,p1
    accu_mono += h0i * reigvec_trans(i,1)/reigvec_trans(1,1)
    pert_mono += h0i*hi0/(h00 - hii)
   else if (degree==2)then
    norm_t2 += (reigvec_trans(i,1)/reigvec_trans(1,1))**2.d0
  !call debug_det(psi_det(1,1,i),N_int)
  !call debug_det(psi_det(1,1,1),N_int)
    call get_double_excitation(psi_det(1,1,i),psi_det(1,1,1),exc,phase,N_int)
    call decode_exc(exc,2,h1,p1,h2,p2,s1,s2)
    print*,'h1,p1,h2,p2',h1,p1,h2,p2
    accu_double += h0i * reigvec_trans(i,1)/reigvec_trans(1,1)
    pert_double += h0i*hi0/(h00 - hii)
   endif
   if(degree.ne.0)then
    print*,'amplitude = ',reigvec_trans(i,1)/reigvec_trans(1,1)
    print*,'ampl pert = ',hi0/(h00 - hii)
    print*,'contrib   = ', h0i * reigvec_trans(i,1)/reigvec_trans(1,1)
    print*,'cont pert = ', h0i * hi0/(h00 - hii)
    print*,'h0i       = ',h0i
    print*,'hi0       = ',hi0
    print*,'Delta E   = ',hii - h00
   endif
   ! correct sum = \sum_i <0|H|i> * t_i
   accu += h0i * reigvec_trans(i,1)/reigvec_trans(1,1)
   ! incorrect sum = \sum_i <i|H|0> * t_i
   accu2+= hi0 * reigvec_trans(i,1)/reigvec_trans(1,1)
  enddo
  print*,''
  print*,''
  print*,'accu wrong      = ',accu2
  print*,'accu            = ',accu
  print*,'eigval_trans    = ',eigval_trans(1)
  print*,'<d0|tile{H}|d0> = ',h00 
  print*,'<d0|   H   |d0> = ',ref_bitmask_energy
  print*,'E corr          = ',eigval_trans(1) - h00
  print*,'accu_mono       = ',accu_mono
  print*,'accu_double     = ',accu_double
  print*,'norm_t1         = ',norm_t1
  print*,'norm_t2         = ',norm_t2
  print*,'pert_mono       = ',pert_mono
  print*,'pert_double     = ',pert_double
  !!call print_mos
  do i = 1, N_det
   psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
  enddo
  touch psi_coef
  call save_wavefunction

end

subroutine print_mos
 implicit none
 integer :: i,nx
 double precision :: dx,xmax,r(3),mos_array(mo_num),r1(3),mu_in ,r12
 double precision :: d_d_r12_11,d_d_r12_12,d_d_r12_13,d_d_r12_16 
 double precision :: short_range
 mu_in = mu_erf 
 xmax = 5.d0
 nx = 1000
 dx = xmax / dble(nx)
 r = 0.d0
 r(1) = 0.d0
 r1 = 0.d0
 r1(1) = 1.5d0
 do i = 1, nx
  call give_all_mos_at_r(r,mos_array)
! call non_hermit_term(r1,r,1,1,mu_in,d_d_r12_11)
! call non_hermit_term(r1,r,1,2,mu_in,d_d_r12_12)
! call non_hermit_term(r1,r,1,3,mu_in,d_d_r12_13)
! call non_hermit_term(r1,r,1,6,mu_in,d_d_r12_16)
  r12 = dabs(r(1)-r1(1))
  short_range = (1.d0 - derf(mu_in * mu_erf))
  write(33,'(100(F16.10,X))')r(1),mos_array(1), mos_array(2), mos_array(3), mos_array(6)
  write(34,'(100(F16.10,X))')r12,d_d_r12_11,d_d_r12_12,d_d_r12_13,d_d_r12_16
  write(35,'(100(F16.10,X))')r12,d_d_r12_11*mos_array(1)*mos_array(2), &
                                 d_d_r12_12*mos_array(1)*mos_array(1), &
                                 d_d_r12_13*mos_array(1)*mos_array(1), &
                                 d_d_r12_16*mos_array(1)*mos_array(1)
  r(1) += dx 
 enddo

end

!subroutine non_hermit_term(r1,r2,i,j,mu_in,d_d_r12)
!implicit none
!integer, intent(in) :: i,j
!double precision, intent(in) :: r1(3), r2(3) , mu_in
!double precision, intent(out):: d_d_r12
!double precision :: mos_array_r1(mo_num),mos_grad_array_r1(3,mo_num)
!double precision :: mos_array_r2(mo_num),mos_grad_array_r2(3,mo_num)
!double precision :: r12(3), dist_r12, dist_vec(3),poly(3)
!double precision :: erf_mu_r12,derf_mu_x,poly_tot(3)
!integer :: k
!call give_all_mos_and_grad_at_r(r1,mos_array_r1,mos_grad_array_r1)
!call give_all_mos_and_grad_at_r(r2,mos_array_r2,mos_grad_array_r2)
!dist_r12 = 0.d0
!do k = 1, 3
! r12(k) = r1(k) - r2(k) 
! dist_r12 += r12(k)*r12(k)
!enddo
!dist_r12 = dsqrt(dist_r12)
!dist_vec(1) = dsqrt(r12(2)*r12(2) + r12(3)*r12(3))
!dist_vec(2) = dsqrt(r12(1)*r12(1) + r12(3)*r12(3))
!dist_vec(3) = dsqrt(r12(1)*r12(1) + r12(2)*r12(2))
!erf_mu_r12 = derf_mu_x(mu_in,dist_r12)
!call inv_r_times_poly(r12, dist_r12, dist_vec, poly)
!! poly_tot(1) = (1 - erf(mu * r12))/(2 * r12) (x1 - x2)
!do k = 1, 3
! poly_tot(k) = 0.5d0 * (poly(k) - erf_mu_r12 * r12(k) )
!enddo
!d_d_r12 = 0.d0
!do k = 1, 3
! d_d_r12 += poly_tot(k) * (mos_grad_array_r1(k,i) - mos_grad_array_r2(k,j) )
!enddo
!end

subroutine test_grad_mo
 implicit none
 integer :: ipoint,i,j
 double precision :: r1(3), accu(3),weight1
 double precision :: mos_array(mo_num),mos_grad_array(3,mo_num)
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  call give_all_mos_and_grad_at_r(r1,mos_array,mos_grad_array)
  weight1 = final_weight_at_r_vector(ipoint)
  do j = 1, mo_num
   do i = 1, 3
    accu(i) += weight1 * dabs(mos_grad_in_r_array_tranp(i,j,ipoint) - mos_grad_array(i,j))
   enddo
  enddo
 enddo
 print*,'accu = '
 print*,accu
end
