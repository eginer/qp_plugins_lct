program test_psi
  implicit none
  read_wf = .True.
  touch read_wf 
!  call routine
  call routine_psi
end

 subroutine routine
  implicit none
  include 'constants.include.F'
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.hole'
  i_unit_output = getUnitAndOpen(output,'w')
  double precision :: rmax,rmin, dr,r(3),geminal,geminal_tot
  integer :: i,nx
  rmin = 0.d0
  rmax = 4.D0
  nx = 10000
  dr = (rmax - rmin)/dble(nx)
  r = 0.d0
  r(1) = rmin
  integer :: i_occ,j_occ
  i_occ = 1
  j_occ = 1
 double precision :: accu_1, accu_2,dm_a,dm_b,mu,rho,geminal_bis
  accu_1 = 0.d0
  accu_2 = 0.d0
  do i = 1, nx
   call geminal_ij_at_r1_r2(r,r,i_occ,j_occ,geminal)
   call geminal_at_r1_r2(r,r,geminal_tot)
!   if(dabs(geminal_tot - geminal).gt.1.d-10)then
!    print*,geminal,geminal_tot,dabs(geminal_tot - geminal)
!    print*,'bug in geminal !'
!    stop
!   endif
   call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
   mu = (-1.d0/geminal)*1.d0/(2.d0*sqpi)
   rho = r(1)**2*(dm_a+dm_b)
   accu_1 += rho*dr 
   accu_2 += rho * mu*dr
   double precision :: mu_mf,dm
   call mu_of_r_mean_field(r,mu_mf, dm)
   write(i_unit_output,'(100(F16.10,X))')r(1),geminal,mu,mu_mf,rho, dm_a+dm_b
   r(1) += dr
  enddo
  print*,'nelec = ',accu_1*4.d0*pi
  print*,'av mu = ',accu_2/accu_1
end


 subroutine routine_psi
  implicit none
  include 'constants.include.F'
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.hole'
  i_unit_output = getUnitAndOpen(output,'w')
  double precision :: rmax,rmin, dr,r1(3),r2(3),geminal,geminal_tot
  double precision :: unit_vec(3),norm
  integer :: i_occ,j_occ
  integer :: i,nx
  double precision :: accu_1, accu_2,mu,rho,geminal_bis
  double precision :: dm_a1,dm_b1,dm_a2,dm_b2,get_norm,r12(3)
  double precision :: grad_dm_a2(3), grad_dm_b2(3), aos_array(ao_num), grad_aos_array(3,ao_num)
  !!!!!!! Unit vector
  unit_vec(1) = 1.d0
  unit_vec(2) = 0.d0
  unit_vec(3) = 0.d0
  norm = get_norm(unit_vec)
  unit_vec *= 1.d0/norm
  !!!!!!!! R1 vector 
  r1 = 0.d0
  r1(1) = 1.0d0
  !!!!!!!! 
  rmin = -2.D0
  rmax = 6.D0
  nx = 10000
  dr = (rmax - rmin)/dble(nx)
  r2 = r1
  r2 += rmin*unit_vec
  i_occ = 1
  j_occ = 1
  accu_1 = 0.d0
  accu_2 = 0.d0
  do i = 1, nx
   call geminal_at_r1_r2(r1,r2,geminal_tot)
   call dm_dft_alpha_beta_at_r(r1,dm_a1,dm_b1)
   call dm_dft_alpha_beta_at_r(r2,dm_a2,dm_b2)
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r& 
   (r2,dm_a2,dm_b2, grad_dm_a2, grad_dm_b2, aos_array, grad_aos_array)
   double precision :: grad_norm
   grad_norm = get_norm(grad_dm_a2)
   r12 = r1 - r2
   norm = get_norm(r12)
   write(i_unit_output,'(100(F16.10,X))')norm,r2(1),geminal_tot,grad_norm,dm_a1,dm_a2 
   r2 += unit_vec * dr
  enddo
end


double precision function get_norm(vec)
 implicit none
 double precision, intent(in) :: vec(3)
 get_norm = vec(1)**2 + vec(2)**2 + vec(3)**2
 get_norm = dsqrt(get_norm)
end
