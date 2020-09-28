program ao_one_e_eff_pot
  implicit none
  call plot_e_vs_mu
! call test_quick
end
subroutine test_quick
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  integer :: i,j 
  provide one_e_j_eigval_trans
  print*,'Eigenvalues '
  do i = 1, mo_num
   print*,'i',one_e_j_eigval_trans(i)
  enddo

  double precision :: accu
  accu = 0.d0
  do i = 1, ao_num
   do j = 1, ao_num
    accu += dabs( ao_erf_integrals_n_e(j,i) - ao_erf_mu_r_inv_r(j,i))
   enddo
  enddo
  print*,'accu = ',accu
  call plot_eff_pot 
end


subroutine plot_eff_pot
 implicit none
 include 'utils/constants.include.F'
 double precision :: r_max,dx,r(3),local_v_ne
 integer :: i,nx
 double precision :: rho, mos_array(mo_num),mu_in
 r_max = 3.d0
 nx = 1000
 dx = r_max/dble(nx)
 r = 0.d0
 mu_in = local_v_ne(r) * 0.5d0 * sqpi
 print*,'mu_in = ',mu_in
 mu_in = (local_v_ne(r) + 0.5d0) * sqpi / 3.d0
 print*,'mu_in = ',mu_in
 do i = 1, nx
  r(1) += dx
  call give_all_mos_at_r(r,mos_array)
  rho = mos_array(1) * mos_array(1)
  write(33,*),r(1),local_v_ne(r),-1.d0/r(1),rho
 enddo
end

subroutine plot_e_vs_mu
 implicit none
 double precision :: dmu,mu_max,mu_in
 integer :: i,nmu
 nmu = 1000
 mu_max = 50.d0
 dmu = mu_max/dble(nmu)
 mu_one_e_j = dmu
 do i = 1, nmu
  write(*,*)mu_one_e_j,one_e_j_eigval_trans(1)  
  write(33,*)mu_one_e_j,one_e_j_eigval_trans(1)  
  mu_one_e_j += dmu
  touch mu_one_e_j
 enddo

end
