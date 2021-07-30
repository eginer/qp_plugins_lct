subroutine fit_erfc_mu(mu, expo_fit, coef_fit)
 implicit none
 double precision, intent(in) :: mu
 double precision, intent(out):: expo_fit(n_max_fit_slat)
 double precision, intent(out):: coef_fit(n_max_fit_slat)
 BEGIN_DOC
! (1 - erf(mu*x)) = \sum_i coef_gauss_1_erf_x(i) * exp(-expo_gauss_1_erf_x(i) * x^2)
!
! This is based on a fit of (1 - erf(mu*x)) by exp(-alpha * x) exp(-beta*mu^2x^2)
!
! and the slater function exp(-alpha * x) is fitted with n_max_fit_slat gaussians 
 END_DOC
 integer :: i
 double precision :: expos(n_max_fit_slat),alpha,beta
 alpha = expos_slat_gauss_1_erf_x(1) * mu_erf 
 call expo_fit_slater_gam(alpha,expos)                                                                                    
 beta = expos_slat_gauss_1_erf_x(2) * mu_erf**2.d0
 do i = 1, n_max_fit_slat
  expo_fit(i) = expos(i) + beta
  coef_fit(i) = coef_fit_slat_gauss(i)
 enddo
end

double precision function gauss_fit_erfc_mu(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision,allocatable :: expo_fit(:)
 double precision,allocatable :: coef_fit(:)
 integer :: i
 allocate( expo_fit(n_max_fit_slat) , coef_fit(n_max_fit_slat))
 call fit_erfc_mu(mu, expo_fit, coef_fit)
 gauss_fit_erfc_mu = 0.d0
 do i = 1, n_max_fit_slat
  gauss_fit_erfc_mu += coef_fit(i) * dexp(-expo_fit(i) * x * x)
 enddo
end

subroutine test_fit_erfc
 implicit none
 integer :: i,nx
 double precision :: mu,fit_erfc_mu,x,xmax,dx
 double precision :: gauss_fit_erfc_mu
 mu = mu_erf
 xmax = 5.D0
 nx = 1000
 dx = xmax/dble(nx)
 x = 0.d0
 do i = 1, nx
  write(33,'(100(F16.10,X))')x,(1.d0 - derf(mu*x)),gauss_fit_erfc_mu(x,mu)
  x += dx
 enddo
end

BEGIN_PROVIDER [double precision, erfc_mu_ij_rk,( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! erfc_mu_ij_rk(i,j,r) = \int dr' phi_i(r') phi_j(r') (1 - erf(mu(r) * |r - r'|))/|r - r'|
!
! it is fitted with fit_erfc_mu
 END_DOC
 double precision :: expo_fit(n_max_fit_slat)
 double precision :: coef_fit(n_max_fit_slat),overlap_gauss_r12_ao
 double precision :: mu,r(3),int_mu,delta,wall0,wall1
 integer :: i,j,ipoint,ifit
  provide mu_of_r_for_ints
  call fit_erfc_mu(1.d0, expo_fit, coef_fit)
 call wall_time(wall0)
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta,ifit,expo_fit,coef_fit) & 
  !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,erfc_mu_ij_rk,final_grid_points,n_max_fit_slat)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
   mu = mu_of_r_for_ints(ipoint,1)
   call fit_erfc_mu(mu, expo_fit, coef_fit)
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   do i = 1, ao_num
    do j = i, ao_num
     erfc_mu_ij_rk(j,i,ipoint) = 0.d0
     do ifit = 1, n_max_fit_slat
      delta = expo_fit(ifit)
      call erf_mu_gauss_ij_ao(i,j,1.d+10, r, delta,int_mu)
      erfc_mu_ij_rk(j,i,ipoint) += coef_fit(ifit) * int_mu 
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    erfc_mu_ij_rk(j,i,ipoint)= erfc_mu_ij_rk(i,j,ipoint)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for erfc_mu_ij_rk  ',wall1 - wall0

END_PROVIDER 

subroutine test_erfc_rij
 implicit none
 integer :: ipoint,i,j
 double precision :: integral, semi_num
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     integral = -v_ij_erf_rk(j,i,ipoint)
     semi_num = erfc_mu_ij_rk(i,j,ipoint)
     if(dabs(integral - semi_num).gt.1.d-10)then
      print*,'ipoint,i,j',ipoint,i,j
      print*,integral,semi_num,dabs(integral - semi_num)
     endif
    enddo
   enddo
  enddo


end
