
subroutine print_energy
 implicit none
 use bitmasks
 integer :: i
 print*,'***************************************'
 
 if(.not.ten_no_jastrow)then
  print*,'Using a mu-based Jastrow'
  print*,'mu_erf = ',mu_erf
 else
  print*,'Using Ten-No Jastrow'
 endif
 print*,'H      eigenvalue'
 i = 1
 print*,'E0       = ',CI_energy(i)
 print*,''
 print*,'Htilde eigenvalue'
 i = 1
 print*,'E0_tilde = ',eigval_trans(i) + nuclear_repulsion
 print*,''
 print*,'Delta E  = ',(eigval_trans(i) + nuclear_repulsion) - CI_energy(i)
end

subroutine print_eigv
 implicit none
 integer :: i,j
 print*,'******************'
 print*,'******************'
 print*,'******************'
 print*,'Overlap betwwen psi_coef and the right eigenvector = ',overlap_psi_det_r_eigevec
 print*,'Eigenvectors '
 print*,'******************'
 print*,'Right eigenvector             Left eigenvector               psi_coef'
 do i = 1, N_det
  write(*,'(I5,X,100(F16.10,X))')i,reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1)),leigvec_trans(i,1)/dsqrt(leigvec_trans_norm(1)),psi_coef(i,1)
 enddo
 double precision :: accu1,accu2,e
 accu1 = 0.d0
 accu2 = 0.d0
 do i = 1, N_det
  accu2 += htilde_matrix_elmt(1,i) * reigvec_trans(i,1) / reigvec_trans(1,1)
  do j = 1, N_det
   accu1 += htilde_matrix_elmt(j,i) * reigvec_trans(i,1) * reigvec_trans(j,1)
  enddo
 enddo
 e = accu1/reigvec_trans_norm(1)
 print*,'******************'
 print*,'Energy computed with the RIGHT eigenvector'
 print*,'-----'
 print*,'expectation value= ',e
 print*,'projection on HF = ',accu2
 print*,''
 print*,'reigvec_trans_norm ',reigvec_trans_norm(1)
 print*,''

 accu1 = 0.d0
 accu2 = 0.d0
 do i = 1, N_det
  accu2 += htilde_matrix_elmt(i,1) * leigvec_trans(i,1) / leigvec_trans(1,1)
  do j = 1, N_det
   accu1 += htilde_matrix_elmt(j,i) * leigvec_trans(i,1) * leigvec_trans(j,1)
  enddo
 enddo
 e = accu1/reigvec_trans_norm(1)
 print*,'******************'
 print*,'Energy computed with the LEFT  eigenvector'
 print*,'-----'
 print*,'expectation value= ',e
 print*,'projection on HF = ',accu2
 print*,''
 print*,'reigvec_trans_norm ',leigvec_trans_norm(1)
end

subroutine print_overlap_left_right
 implicit none
 integer :: i
 print*,'***************'
 print*,'***************'
 print*,'OVERLAP PRINTING'
 print*,'***************'
 print*,'Printing the overlap between the left and right eigenvectors'
 do i = 1, n_tc_ovlp_print
  write(*,'(I4,X,100(F9.5,X))')i,left_right_overlap(i,1:min(100,n_tc_ovlp_print))
 enddo

 print*,'               '
 print*,'***************'
 print*,'Printing the overlap between the right and right eigenvectors'
 do i = 1, n_tc_ovlp_print
  write(*,'(I4,X,100(F9.5,X))')i,right_right_overlap(i,1:min(100,n_tc_ovlp_print))
 enddo

 print*,'               '
 print*,'***************'
 print*,'Printing the overlap between the left and left eigenvectors'
 do i = 1, n_tc_ovlp_print
  write(*,'(I4,X,100(F9.5,X))')i,left_left_overlap(i,1:min(100,n_tc_ovlp_print))
 enddo


end

subroutine print_pert
 implicit none
 integer :: i,j
 double precision :: accu,hmono,herf,heff,hderiv,htot
 double precision :: accu_mono,accu_double,accu2
 double precision :: pert_mono,pert_double
 double precision :: h00,hii,htotbis,phase,h0i,hi0,hthree
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
 call htilde_mat(psi_det(1,1,1),psi_det(1,1,1),hmono,herf,heff,hderiv,hthree,h00)
 print*,''
 print*,'*************************************'
 print*,'*************************************'
 print*,'*************************************'
 print*,''
 print*,'Detailed analysis by projection on HF'
 print*,''
 print*,''
  print*,'Printing the connection '
  do i = 1, N_det
   ! <0|H|i>
   print*,'******************'
   print*,'i = ',i
   call htilde_mat(psi_det(1,1,1),psi_det(1,1,i),hmono,herf,heff,hderiv,hthree,h0i)
 
   ! <i|H|i>
   call htilde_mat(psi_det(1,1,i),psi_det(1,1,i),hmono,herf,heff,hderiv,hthree,hii)
 
   ! <i|H|0>
   call htilde_mat(psi_det(1,1,i),psi_det(1,1,1),hmono,herf,heff,hderiv,hthree,hi0)
 
   call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,1),degree,N_int)
   print*,'degree = ',degree
   if(degree==1)then
    call get_single_excitation(psi_det(1,1,i),psi_det(1,1,1),exc,phase,N_int)
    call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
    norm_t1 += (reigvec_trans(i,1)/reigvec_trans(1,1))**2.d0
    print*,'h1,p1',h1,p1
    accu_mono += h0i * reigvec_trans(i,1)/reigvec_trans(1,1)
    pert_mono += h0i*hi0/(h00 - hii)
   else if (degree==2)then
    norm_t2 += (reigvec_trans(i,1)/reigvec_trans(1,1))**2.d0
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
  print*,'projection on HF= ',accu
  print*,'eigval_trans    = ',eigval_trans(1)
  print*,''
  print*,'<d0|tile{H}|d0> = ',h00 
  print*,'<d0|   H   |d0> = ',ref_bitmask_energy
  print*,''
  print*,'E corr          = ',eigval_trans(1) - h00
  print*,''
  print*,'Analysis of energy by components '
  print*,'accu_mono       = ',accu_mono
  print*,'accu_double     = ',accu_double
  print*,'Analysis of energy by perturbation'
  print*,'pert_mono       = ',pert_mono
  print*,'pert_double     = ',pert_double
  print*,'Norm of amplitudes by components '
  print*,'norm_t1         = ',norm_t1
  print*,'norm_t2         = ',norm_t2
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

subroutine non_hermit_term(r1,r2,i,j,mu_in,d_d_r12)
 implicit none
 integer, intent(in) :: i,j
 double precision, intent(in) :: r1(3), r2(3) , mu_in
 double precision, intent(out):: d_d_r12
 double precision :: mos_array_r1(mo_num),mos_grad_array_r1(3,mo_num)
 double precision :: mos_array_r2(mo_num),mos_grad_array_r2(3,mo_num)
 double precision :: r12(3), dist_r12, dist_vec(3),poly(3)
 double precision :: erf_mu_r12,derf_mu_x,poly_tot(3)
 integer :: k
 call give_all_mos_and_grad_at_r(r1,mos_array_r1,mos_grad_array_r1)
 call give_all_mos_and_grad_at_r(r2,mos_array_r2,mos_grad_array_r2)
 dist_r12 = 0.d0
 do k = 1, 3
  r12(k) = r1(k) - r2(k) 
  dist_r12 += r12(k)*r12(k)
 enddo
 dist_r12 = dsqrt(dist_r12)
 dist_vec(1) = dsqrt(r12(2)*r12(2) + r12(3)*r12(3))
 dist_vec(2) = dsqrt(r12(1)*r12(1) + r12(3)*r12(3))
 dist_vec(3) = dsqrt(r12(1)*r12(1) + r12(2)*r12(2))
 erf_mu_r12 = derf_mu_x(mu_in,dist_r12)
 call inv_r_times_poly(r12, dist_r12, dist_vec, poly)
 ! poly_tot(1) = (1 - erf(mu * r12))/(2 * r12) (x1 - x2)
 do k = 1, 3
  poly_tot(k) = 0.5d0 * (poly(k) - erf_mu_r12 * r12(k) )
 enddo
 d_d_r12 = 0.d0
 do k = 1, 3
  d_d_r12 += poly_tot(k) * (mos_grad_array_r1(k,i) - mos_grad_array_r2(k,j) )
 enddo
end

subroutine plot_on_top_left_right
 implicit none
 integer :: i
 integer                        :: i_unit_output,getUnitAndOpen
 character*(128)                :: output
 PROVIDE ezfio_filename
 
 provide reigvec_trans_norm
 print*,'Psi '
 output=trim(ezfio_filename)//'.on_top_psi'
 i_unit_output = getUnitAndOpen(output,'w')
 call plot_on_top(i_unit_output)

 print*,'Right eigenvector'
 output=trim(ezfio_filename)//'.on_top_right'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, N_det
  psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
 enddo
 touch psi_coef
 call plot_on_top(i_unit_output)


 print*,'Left eigenvector'
 output=trim(ezfio_filename)//'.on_top_left'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, N_det
  psi_coef(i,1) = leigvec_trans(i,1)/dsqrt(leigvec_trans_norm(1))
 enddo
 touch psi_coef
 call plot_on_top(i_unit_output)

end

subroutine plot_on_top(i_unit)
 implicit none
 integer, intent(in) :: i_unit
 integer :: i,nx,istate 
 double precision :: x,xmax,dx,r(3)
 double precision :: on_top_in_r,dm_a,dm_b
 istate = 1
 nx = 1000
 xmax = 3.d0
 dx = xmax / dble(nx)
 r = 0.d0
 do i = 1, nx
  call give_on_top_in_r_one_state(r,istate,on_top_in_r)
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  write(i_unit,'(10(F16.10,X))')r(1), on_top_in_r , dm_a+dm_b
  r(1) += dx 
 enddo
end
