subroutine test_one_e_int
 implicit none
 double precision :: mu_in, D_center(3)
 integer :: i,j,ipoint,m
 double precision :: r1(3), weight,accu_xyz(3),accu_dxyz(3),r12
 double precision :: num_int,exact,err_abs,accu_abs,err_relat,accu_relat
 double precision :: aos_array_r1(ao_num),aos_grad_array_r1(3,ao_num),derf_mu_x,xyz_ints(3),dxyz_ints(3)
 double precision :: coulomb,accu_coulomb,phi_j_erf_mu_r_phi
 mu_in = mu_erf
 D_center = 0.d0
 D_center(1) = 1.2d0
 D_center(3) = -0.1d0

 do i = 1, ao_num
  do j = 1, ao_num
   call phi_j_erf_mu_r_xyz_phi(i,j,mu_in, D_center, xyz_ints)
   call phi_j_erf_mu_r_dxyz_phi(i,j,mu_in, D_center, dxyz_ints)
   coulomb = phi_j_erf_mu_r_phi(i,j,mu_in, D_center)
   accu_coulomb = 0.d0
   accu_xyz = 0.d0
   accu_dxyz= 0.d0
   do ipoint = 1, n_points_final_grid
    r1(1) = final_grid_points(1,ipoint)
    r1(2) = final_grid_points(2,ipoint)
    r1(3) = final_grid_points(3,ipoint)
    r12   = (r1(1) - D_center(1))**2.d0 + (r1(2) - D_center(2))**2.d0 + (r1(3) - D_center(3))**2.d0
    r12   = dsqrt(r12)
    call give_all_aos_and_grad_at_r(r1,aos_array_r1,aos_grad_array_r1)
    weight = final_weight_at_r_vector(ipoint)
    accu_coulomb += weight * aos_array_r1(j) * aos_array_r1(i) * derf_mu_x(mu_in,r12)
    do m = 1, 3
     accu_xyz(m)  += weight * aos_array_r1(j) * aos_array_r1(i) * r1(m) * derf_mu_x(mu_in,r12)
     accu_dxyz(m) += weight * aos_array_r1(j) * aos_grad_array_r1(m,i)  * derf_mu_x(mu_in,r12)
    enddo
   enddo

   num_int = accu_coulomb
   exact   = coulomb
   err_abs = dabs(num_int - exact)
   err_relat = 0.d0
   if(dabs(num_int).gt.1.d-09)then
    err_relat = err_abs/dabs(num_int)
   endif
   if(err_relat .gt. 1.d-3)then
    print*,'AHAHAH'
    print*,'Coulomb'
    print*,'i,j',i,j
    print*,'num_int,exact    ',num_int,exact
    print*,'err_abs,err_relat',err_abs,err_relat
!    stop
   endif

   do m = 1, 3
    num_int = accu_xyz(m)
    exact   = xyz_ints(m)
    err_abs = dabs(num_int - exact)
    err_relat = 0.d0
    if(dabs(num_int).gt.1.d-09)then
     err_relat = err_abs/dabs(num_int)
    endif
    if(err_relat .gt. 1.d-3)then
     print*,'AHAHAH'
     print*,'xyz '
     print*,'i,j',i,j
     print*,'num_int,exact    ',num_int,exact
     print*,'err_abs,err_relat',err_abs,err_relat
!     stop
    endif

    num_int = accu_dxyz(m)
    exact   = dxyz_ints(m)
    err_abs = dabs(num_int - exact)
    err_relat = 0.d0
    if(dabs(num_int).gt.1.d-09)then
     err_relat = err_abs/dabs(num_int)
    endif
    if(err_relat .gt. 1.d-3)then
     print*,'AHAHAH'
     print*,'dxyz '
     print*,'i,j',i,j
     print*,'num_int,exact    ',num_int,exact
     print*,'err_abs,err_relat',err_abs,err_relat
!     stop
    endif
   enddo
  enddo
 enddo

end
