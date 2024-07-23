program test_psi
  implicit none
  read_wf = .True.
  touch read_wf 
  call routine
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
  double precision :: rmax,rmin, dr,r(3),geminal
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
 double precision :: accu_1, accu_2,dm_a,dm_b,mu,rho
  accu_1 = 0.d0
  accu_2 = 0.d0
  do i = 1, nx
   call geminal_at_r1_r2(r,r,i_occ,j_occ,geminal)
   call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
   mu = (-1.d0/geminal)*1.d0/(2.d0*sqpi)
   rho = r(1)**2*(dm_a+dm_b)
   accu_1 += rho*dr 
   accu_2 += rho * mu*dr
   write(i_unit_output,'(100(F16.10,X))')r(1),geminal,mu,rho, dm_a+dm_b
   r(1) += dr
  enddo
  print*,'nelec = ',accu_1*4.d0*pi
  print*,'av mu = ',accu_2/accu_1
end
