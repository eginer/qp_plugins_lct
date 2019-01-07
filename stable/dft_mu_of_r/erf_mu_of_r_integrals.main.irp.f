program erf_mu_of_r_integrals
 implicit none 
!call test_int_erf_bielec_ijkl_mo
 provide ao_bielec_integrals_erf_mu_of_r_in_map
 integer :: i,j,k,l
 double precision :: integral,integral_2,integral_3
 double precision :: get_ao_bielec_integral_erf_mu_of_r, erf_mu_of_r_ao
 do i =1, ao_num!1
  do j = 1, ao_num!2
   do k = 1, ao_num!1
    do l =1,  ao_num!2
     integral   = get_ao_bielec_integral_erf_mu_of_r(i,j,k,l,ao_integrals_erf_mu_of_r_map)
     integral_2 = integral
     integral_3 = erf_mu_of_r_ao(i,k,j,l) 
     if(dabs(0.5D0 * (integral + integral_2) - integral_3) .GT.1.d-10)then
      print*,'ijkl',i,j,k,l
      print*,'integral, integral_left, integral correct'
      print*,integral , integral_2, integral_3  
      print*,dabs(0.5D0 * (integral + integral_2) - integral_3)
     endif
    enddo
   enddo
  enddo
 enddo
 !call test_mu_erf
!call test_int_erf_bielec_ijkl
!call test_aos
!call bis
!call test_prim
!call test_integrals_rint
!call test_erf
end

!subroutine test_mu_erf
!implicit none
!double precision, allocatable :: integrals_mo(:,:),mos_array(:)
!double precision :: r(3),mu_in
!allocate(integrals_mo(mo_tot_num,mo_tot_num),mos_array(mo_tot_num))
!r = 0.d0
!mu_in = 0.001d0
!integer :: i_mu,n_mui,j,k,l
!double precision :: dmu,mu_tot
!mu_tot = 2.d0
!n_mu = 100
!do i_mu = 1, n_mu
! call give_all_erf_mu_of_r_kl_mo(integrals_mo,mu_in,r)
! call give_all_mos_at_r(r,mos_array)
! tmp = 0.d0
! do i = 1, elec_alpha_num
!  do j = 1, elec_beta_num
!   tmp += mos_array(j)**2 * integrals_mo(i,i)
!  enddo
! enddo
!enddo
!    

!end

subroutine test_erf
 implicit none
 integer :: i,nt
 double precision :: integral,t,dt,x,mu_in

 x = 0.5d0
 mu_in = 0.5d0


 nt = 1000000
 dt = mu_in/dble(nt)

 t=0.d0
 integral = 0.d0
 do i = 1, nt
  integral += dexp(-x**2 * t**2)
  t+=dt
 enddo
 integral *= 2.d0/dsqrt(dacos(-1.d0)) * dt
 print*,'integral,erf'
 print*,integral,erf(mu_in*x)/x

end

subroutine test_integrals_rint
 implicit none
 integer :: i,nx
 double precision :: alpha,beta,dx,dxalpha
 alpha = 10.d0
 nx = 1000000

 beta = 1.0d0
 dxalpha = alpha/dble(nx)
 dx = 1.d0/dble(nx)
 double precision :: rint,integral,integral_alpha
 double precision :: x
 integral_alpha = 0.d0
 x = 0.d0
 do i = 1, nx
  integral_alpha +=  dexp(-beta * x**2)
  x += dxalpha
 enddo
 integral_alpha*= dxalpha

 x = 0.d0
 integral = 0.d0
 do i = 1, nx
  integral +=  dexp(-beta * alpha **2 * x**2)
 enddo
 integral *= dx / alpha
 print*,'integral alpha , integral'
 print*,integral_alpha, integral
 print*,'rint'
 print*,rint(0,beta*alpha**2) * alpha
 print*,'dsqrt(pi.)/2'
 print*,dsqrt(dacos(-1.d0))/2.d0
 



end


subroutine test_prim
 implicit none
  
 integer :: j,k,l,i
 integer :: m,n
 double precision :: dist,mu_in
 double precision, allocatable :: r(:), aos_array(:), aos_integrals(:,:),aos_integrals_bis(:,:),integral_bourrin(:,:)
 double precision               :: alpha, beta, weight
 integer                        :: num_A,num_B,n_pt_in
 double precision               :: A_center(3),B_center(3),C_center(3),P_center(3)
 integer                        :: power_A(3),power_B(3)
 double precision :: primitive_value
 integer :: nu,i_u
 double precision :: du,u
 double precision :: E_ab, p,rho
 integer :: i_x,nx
 double precision :: x,dx,domain,E_ab_vec(3)
 double precision :: Int_x(3),dist_ax,dist_bx,dist_cx,dist_px
 integer :: i_ao,j_ao
 integer :: nt,i_t
 double precision :: dt,t
 double precision :: Int_Ox(3),dist_pc
 allocate(r(3), aos_array(ao_num),integral_bourrin(ao_num,ao_num),aos_integrals(ao_num,ao_num),aos_integrals_bis(ao_prim_num_max,ao_prim_num_max))
 C_center = 0.00000000d0
 mu_in = 0.01d0
 aos_integrals = 0.d0
 i_ao = 1
 j_ao = 2
 do m = 1, ao_prim_num(i_ao)
  alpha = ao_expo_ordered_transp(m,i_ao)
  do n = 1, ao_prim_num(j_ao)
   beta = ao_expo_ordered_transp(n,j_ao)
   do j = 1, nucl_num
    do k = 1, n_points_radial_grid  -1
     do l = 1, n_points_integration_angular 
      r(1) = grid_points_per_atom(1,l,k,j)
      r(2) = grid_points_per_atom(2,l,k,j)
      r(3) = grid_points_per_atom(3,l,k,j)
      weight = final_weight_functions_at_grid_points(l,k,j)
      dist = 0.d0
      do i = 1, 3
       dist += (r(i) - C_center(i))**2
      enddo
      dist = max(dist,1.d-10)
      dist = dsqrt(dist)
      aos_integrals(n,m) += primitive_value(i_ao,m,r) * primitive_value(j_ao,n,r) * weight * erf(mu_in * dist)/dist
      enddo
     enddo
    enddo
!   print*,'integrals good = ',aos_integrals(n,m)
  enddo
 enddo

!integral_bourrin = 0.d0
!integer :: nt,i_t
!double precision :: dt,t
!nt = 10000
!dt = mu_in/dble(nt)
!do m = 1, ao_prim_num(i_ao)
! alpha = ao_expo_ordered_transp(m,i_ao)
! do n = 1, ao_prim_num(j_ao)
!  beta = ao_expo_ordered_transp(n,j_ao)
!  t = 0.d0
!  do i_t = 1, nt
!   do j = 1, nucl_num
!    do k = 1, n_points_radial_grid  -1
!     do l = 1, n_points_integration_angular 
!      r(1) = grid_points_per_atom(1,l,k,j)
!      r(2) = grid_points_per_atom(2,l,k,j)
!      r(3) = grid_points_per_atom(3,l,k,j)
!      weight = final_weight_functions_at_grid_points(l,k,j)
!      dist = 0.d0
!      do i = 1, 3
!       dist += (r(i) - C_center(i))**2
!      enddo
!      dist = max(dist,1.d-10)
!      dist = dsqrt(dist)
!      integral_bourrin(n,m) += primitive_value(i_ao,m,r) * primitive_value(j_ao,n,r) * weight * dexp(-dist**2*t**2)
!      enddo
!     enddo
!    enddo
!   t += dt
!  enddo
!  integral_bourrin(n,m) *= dt * 2.d0/dsqrt(dacos(-1.d0))
! enddo
!enddo

!do m = 1, ao_prim_num(i_ao)
! do n = 1, ao_prim_num(j_ao)
!  if(dabs(aos_integrals(n,m)-integral_bourrin(n,m)).gt.1.d-6)then
!   print*,'prb !!'
!   print*,m,n
!   print*,dabs(aos_integrals(n,m)-integral_bourrin(n,m))
!   print*,aos_integrals(n,m),integral_bourrin(n,m)
!  endif
! enddo
!enddo



!!integration over U 
!integral_bourrin = 0.d0
!if(.False.)then
! print*,'Integration over U'
! nu = 100000
! du = mu_in/dble(nu)
! nx = 1000
! domain = 20.d0
! dx = domain/dble(nx)
! do m = 1, ao_prim_num(i_ao)
!  num_A = ao_nucl(i_ao)
!  power_A(1:3)= ao_power(i_ao,1:3)
!  A_center(1:3) = nucl_coord(num_A,1:3)
!  alpha = ao_expo_ordered_transp(m,i_ao)
!  do n = 1, ao_prim_num(j_ao)
!   beta = ao_expo_ordered_transp(n,j_ao)
!   num_B = ao_nucl(j_ao)
!   power_B(1:3)= ao_power(j_ao,1:3)
!   B_center(1:3) = nucl_coord(num_B,1:3)
!   p = alpha + beta 
!   rho = alpha * beta / p
!   E_ab_vec = 0.d0
!   E_ab = 1.d0
!   do i = 1, 3
!    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) / p
!    E_ab_vec(i) = dexp(-rho * (A_center(i)-B_center(i))**2)
!    E_ab *= E_ab_vec(i)
!   enddo
!   u = 0.d0
!   integral = 0.d0
!   do i_u = 1, nu
!    ! integration over x, y and z
!    do i = 1, 3
!     Int_x(i) = 0.d0
!     x = -domain * 0.5d0
!     do i_x = 1, nx
!      dist_ax =  x-A_center(i)
!      dist_bx =  x-B_center(i)
!      dist_cx =  x-C_center(i)
!      dist_px =  x-P_center(i)
! !    Int_x(i) += dist_ax**power_A(i) * dist_bx**power_B(i) * dexp(-alpha*dist_ax**2) * dexp(-beta*dist_bx**2) * dexp(-u**2*dist_cx**2)
!      Int_x(i) += dist_ax**power_A(i) * dist_bx**power_B(i) * dexp(-p * dist_px**2) * dexp(-u**2*dist_cx**2) * E_ab_vec(i)
!      x += dx
!     enddo
!     Int_x(i) *= dx
!    enddo
!    integral += Int_x(1) * Int_x(2) * Int_x(3) * 2.d0/dsqrt(dacos(-1.d0)) * du
!    u+= du
!   enddo
!   if(dabs(aos_integrals(n,m)-integral).gt.1.d-6)then
!    print*,'prb !!'
!    print*,m,n
!    print*,dabs(aos_integrals(n,m)/integral)
!    print*,aos_integrals(n,m),integral
!   endif
!  enddo
! enddo
!endif
 

!if(.False.)then
! print*,'Integration over t'
!! integration over t 
! integral_bourrin = 0.d0
! nt = 100
!!dt = 1.d0/dble(nt)
!
! nx = 10000000
! domain = 50.d0
! dx = domain/dble(nx)
! do m = 1, ao_prim_num(i_ao)
!  num_A = ao_nucl(i_ao)
!  power_A(1:3)= ao_power(i_ao,1:3)
!  A_center(1:3) = nucl_coord(num_A,1:3)
!  alpha = ao_expo_ordered_transp(m,i_ao)
!  do n = 1, ao_prim_num(j_ao)
!   beta = ao_expo_ordered_transp(n,j_ao)
!   num_B = ao_nucl(j_ao)
!   power_B(1:3)= ao_power(j_ao,1:3)
!   B_center(1:3) = nucl_coord(num_B,1:3)
!!  numerical integration over x
!   p = alpha + beta 
!   rho = alpha * beta / p
!   E_ab_vec = 0.d0
!   E_ab = 1.d0
!   dist_pc = 0.d0
!   do i = 1, 3
!    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) / p
!    E_ab_vec(i) = dexp(-rho * (A_center(i)-B_center(i))**2)
!    E_ab *= E_ab_vec(i)
!    dist_pc += (P_center(i) - C_center(i))**2
!   enddo
!   t = 0.d0
!   integral = 0.d0
!   dt = (mu_in / dsqrt(p+ mu_in**2))/dble(nt)
!   do i_t = 1, nt
!    ! integration over x, y and z
!    do i = 1, 3
!     Int_x(i) = 0.d0
!     x = -domain * 0.5d0
!     do i_x = 1, nx
!      dist_ax =  x-A_center(i)
!      dist_bx =  x-B_center(i)
!      dist_cx =  x-C_center(i)
!      dist_px =  x-P_center(i)
!      Int_Ox(i) += dexp(-p * dist_px**2) * dexp(-t**2 * rho / (1.d0 -t**2)*dist_cx**2) * E_ab_vec(i)
!      Int_x(i)  += dist_ax**power_A(i) * dist_bx**power_B(i) * dexp(-p * dist_px**2) * dexp(-t**2 * rho / (1.d0 -t**2)*dist_cx**2) * E_ab_vec(i)
!      x += dx
!     enddo
!     if(isnan(Int_Ox(i)))then
!      print*,'AAHH NAN'
!      print*,'************'
!      print*,x,dexp(-p * dist_px**2) * dexp(-t**2 * rho / (1.d0 -t**2)*dist_cx**2) * E_ab_vec(i)
!      print*,dexp(-p * dist_px**2) 
!      print*,dexp(-t**2 * rho / (1.d0 -t**2)*dist_cx**2) 
!      print*,E_ab_vec(i) 
!      print*,'************'
!      stop
!     endif
!     Int_x(i) *= dx
!     Int_Ox(i) *= dx
!    enddo
!    integral += 2.d0*dacos(-1.d0) / p * E_ab * Int_x(1)/Int_Ox(1) * Int_x(2)/Int_Ox(2) * Int_x(3)/Int_Ox(3) *  dexp(-p*t**2*dist_pc) *  dt 
!    t+= dt
!   enddo
!
!   if(dabs(aos_integrals(n,m)-integral).gt.1.d-6)then
!    print*,'prb !!'
!    print*,m,n
!    print*,dabs(aos_integrals(n,m)/integral)
!    print*,aos_integrals(n,m),integral
!   endif
!  enddo
! enddo
!endif
 
!print*,'Analytic integrals'
!do m = 1, ao_prim_num(i_ao)
! alpha = ao_expo_ordered_transp(m,i_ao)
! num_A = ao_nucl(i_ao)
! power_A(1:3)= ao_power(i_ao,1:3)
! A_center(1:3) = nucl_coord(num_A,1:3)
! do n = 1, ao_prim_num(j_ao)
!  beta = ao_expo_ordered_transp(n,j_ao)
!  num_B = ao_nucl(j_ao)
!  power_B(1:3)= ao_power(j_ao,1:3)
!  B_center(1:3) = nucl_coord(num_B,1:3)

!  p = alpha + beta 
!  rho = alpha * beta / p
!  E_ab_vec = 0.d0
!  E_ab = 1.d0

!  do i = 1, 3
!   P_center(i) = (alpha * A_center(i) + beta * B_center(i)) / p
!   E_ab_vec(i) = dexp(-rho * (A_center(i)-B_center(i))**2)
!   E_ab *= E_ab_vec(i)
!  enddo
!  integral  = 2.d0 * dacos(-1.d0) / p * E_ab * mu_in/dsqrt(p+mu_in**2)
!  if(dabs(aos_integrals(n,m)-integral).gt.1.d-6)then
!   print*,'prb !!'
!   print*,m,n
!   print*,dabs(aos_integrals(n,m)/integral)
!   print*,aos_integrals(n,m),integral
!  endif
!
! enddo
!enddo

    !call give_all_aos_at_r(r,aos_array)
    !integral = 0.d0
    !do m = 1, ao_prim_num(i_ao)
    ! integral +=  ao_coef_normalized_ordered_transp(m,i_ao) * primitive_value(i_ao,m,r) 
    !enddo
    !if(dabs(integral - aos_array(i_ao)).gt.1.d-10)then
    ! print*,'pb !! '
    ! print*,integral,aos_array(i_ao)
    !endif
    !integral = 0.d0
    !do m = 1, ao_prim_num(j_ao)
    ! integral +=  ao_coef_normalized_ordered_transp(m,j_ao) * primitive_value(j_ao,m,r) 
    !enddo
    !if(dabs(integral - aos_array(j_ao)).gt.1.d-10)then
    ! print*,'pb !! '
    ! print*,integral,aos_array(j_ao)
    !endif
! do m = 1, ao_num
!  do n = 1, ao_num
!   if(dabs(aos_integrals(n,m) + ao_nucl_elec_integral(n,m)/nucl_charge(1)).gt.1.d-6)then
!    print*,'prb !!'
!    print*,n,m
!    print*,ao_nucl_elec_integral(n,m)/nucl_charge(1),aos_integrals(n,m)
!   endif
!  enddo
! enddo

  double precision  :: NAI_pol_mult_erf
  double precision  :: NAI_pol_mult
  double precision :: integral
  print*,'NAI pol mult erf '
  n_pt_in = n_pt_max_integrals
  do l = 1, ao_prim_num(j_ao)
   num_A = ao_nucl(i_ao)
   power_A(1:3)= ao_power(i_ao,1:3)
   A_center(1:3) = nucl_coord(num_A,1:3)
   alpha = ao_expo_ordered_transp(l,i_ao)
   do m = 1, ao_prim_num(i_ao)
    beta = ao_expo_ordered_transp(m,j_ao)
    num_B = ao_nucl(j_ao)
    power_B(1:3)= ao_power(j_ao,1:3)
    B_center(1:3) = nucl_coord(num_B,1:3)
    integral =  NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
    aos_integrals_bis(m,l) = integral
   enddo
  enddo

 do m = 1, ao_prim_num(i_ao)
  do n = 1, ao_prim_num(j_ao)
  !if(dabs(aos_integrals(n,m)-aos_integrals_bis(n,m)).gt.1.d-6)then
  ! print*,'prb !!'
  ! print*,m,n
    print*,aos_integrals(n,m),aos_integrals_bis(n,m),dabs(aos_integrals(n,m)/aos_integrals_bis(n,m))
  !endif
  enddo
 enddo
 deallocate(r,aos_array)
end

subroutine bis
 implicit none
 double precision  :: NAI_pol_mult_erf
 double precision  :: NAI_pol_mult
 double precision :: integral
 double precision :: A_center(3), B_center(3), C_center(3), r(3), mu_in, alpha,beta,weight,dist
 integer ::i_ao,j_ao,l,m,n_pt_in,n,j,k,u,i
 integer :: power_A(3), power_B(3),num_A,num_B
 double precision, allocatable :: aos_array(:), aos_integrals(:,:),aos_integrals_bis(:,:)
 double precision :: dmu ,x
 double precision :: domain,dx,primitive_value,mu_max
 integer :: nx,i_x,nmu,i_mu
 allocate(aos_array(ao_num),aos_integrals(ao_num,ao_num),aos_integrals_bis(ao_prim_num_max,ao_prim_num_max))

 i_ao = 1
 j_ao = 1
 C_center = 0.d0
 C_center(1) = 0.13d0
 nx = 10
 domain = 10.d0
 mu_max = 10.d0
 nmu = 10
 dmu = mu_max/dble(nmu)
 dx = domain/dble(nx)
 do i_ao = 1, ao_num
  do j_ao = 1, ao_num
   mu_in = 0.01d0
   print*,i_ao,j_ao
   do i_mu = 1, nmu
    x = 0.d0
    do i_x = 1, nx
!    print*,'C_center '
!    print*,C_center 
     aos_integrals = 0.d0
     do m = 1, ao_prim_num(i_ao)
      alpha = ao_expo_ordered_transp(m,i_ao)
      do n = 1, ao_prim_num(j_ao)
       beta = ao_expo_ordered_transp(n,j_ao)
       do j = 1, nucl_num
        do k = 1, n_points_radial_grid  -1
         do l = 1, n_points_integration_angular 
          r(1) = grid_points_per_atom(1,l,k,j)
          r(2) = grid_points_per_atom(2,l,k,j)
          r(3) = grid_points_per_atom(3,l,k,j)
          weight = final_weight_functions_at_grid_points(l,k,j)
          dist = 0.d0
          do i = 1, 3
           dist += (r(i) - C_center(i))**2
          enddo
          dist = max(dist,1.d-10)
          dist = dsqrt(dist)
          aos_integrals(n,m) += primitive_value(i_ao,m,r) * primitive_value(j_ao,n,r) * weight * erf(mu_in * dist)/dist
          enddo
         enddo
        enddo
      enddo
     enddo
    
     n_pt_in = n_pt_max_integrals
     do m = 1, ao_prim_num(i_ao)
      alpha = ao_expo_ordered_transp(m,i_ao)
      num_A = ao_nucl(i_ao)
      power_A(1:3)= ao_power(i_ao,1:3)
      A_center(1:3) = nucl_coord(num_A,1:3)
      do n = 1, ao_prim_num(j_ao)
       beta = ao_expo_ordered_transp(n,j_ao)
       num_B = ao_nucl(j_ao)
       B_center(1:3) = nucl_coord(num_B,1:3)
       power_B(1:3)= ao_power(j_ao,1:3)
       integral =  NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
       aos_integrals_bis(n,m) = integral
      enddo
     enddo
     do m = 1, ao_prim_num(i_ao)
      do n = 1, ao_prim_num(j_ao)
       if(dabs(aos_integrals(n,m)-aos_integrals_bis(n,m)).gt.1.d-6 )then
        print*,'prb !!'
        print*,'mu_in',mu_in
        print*,'C_center = '
        print*,C_center
        print*,'i_ao,j_ao = ',i_ao,j_ao
        print*,m,n
        print*,aos_integrals(n,m),aos_integrals_bis(n,m),dabs(aos_integrals(n,m)/aos_integrals_bis(n,m))
       endif
       if(isnan(aos_integrals_bis(n,m)))then
        print*,'NANANANANANA'
        print*,'mu_in',mu_in
        print*,'C_center = '
        print*,C_center
        print*,m,n
        print*,aos_integrals(n,m),aos_integrals_bis(n,m),dabs(aos_integrals(n,m)/aos_integrals_bis(n,m))
       endif
      enddo
     enddo
     x += dx
     C_center(2) = x
   
    enddo ! i_x
    mu_in += dmu
   enddo ! i_mu
  enddo ! j_ao
 enddo ! i_ao

end



subroutine test_aos
 implicit none
 double precision  :: NAI_pol_mult_erf_ao
 double precision :: integral
 double precision :: A_center(3), B_center(3), C_center(3), r(3), mu_in, alpha,beta,weight,dist
 integer ::i_ao,j_ao,l,m,n_pt_in,n,j,k,u,i
 integer :: power_A(3), power_B(3),num_A,num_B
 double precision, allocatable :: aos_array(:)
 double precision :: dmu ,x
 double precision :: domain,dx,primitive_value,mu_max,integral_brut
 integer :: nx,i_x,nmu,i_mu
 allocate(aos_array(ao_num))

 C_center = 0.d0
 C_center(1) = 0.13d0
 nx = 3
 domain = 10.d0
 mu_max = 10.d0
 nmu = 3
 dmu = mu_max/dble(nmu)
 dx = domain/dble(nx)
 do i_ao = 20, ao_num
  do j_ao = 1, ao_num
   mu_in = 0.01d0
   print*,i_ao,j_ao
   do i_mu = 1, nmu
    x = 0.d0
    do i_x = 1, nx
     integral_brut = 0.d0
     do j = 1, nucl_num
      do k = 1, n_points_radial_grid  -1
       do l = 1, n_points_integration_angular 
        r(1) = grid_points_per_atom(1,l,k,j)
        r(2) = grid_points_per_atom(2,l,k,j)
        r(3) = grid_points_per_atom(3,l,k,j)
        weight = final_weight_functions_at_grid_points(l,k,j)
        dist = 0.d0
        do i = 1, 3
         dist += (r(i) - C_center(i))**2
        enddo
        call give_all_aos_at_r(r,aos_array)
        dist = max(dist,1.d-10)
        dist = dsqrt(dist)
        integral_brut += aos_array(i_ao) * aos_array(j_ao) * weight * erf(mu_in * dist)/dist
       enddo
      enddo
     enddo
     integral = NAI_pol_mult_erf_ao(i_ao,j_ao,mu_in,C_center)  
     if(dabs(integral - integral_brut).gt.1.d-6 )then
      print*,'prb !!'
      print*,'mu_in',mu_in
      print*,'C_center = '
      print*,C_center
      print*,'i_ao,j_ao = ',i_ao,j_ao
      print*,i_ao,j_ao
      print*,integral , integral_brut,dabs(integral / integral_brut)
     endif
     if(isnan(integral))then
      print*,'NANANANANANA'
      print*,'mu_in',mu_in
      print*,'C_center = '
      print*,C_center
      print*,i_ao,j_ao
      print*,integral , integral_brut,dabs(integral / integral_brut)
     endif
      print*,integral , integral_brut,dabs(integral / integral_brut)
     x += dx
     C_center(2) = x
   
    enddo ! i_x
    mu_in += dmu
   enddo ! i_mu
  enddo ! j_ao
 enddo ! i_ao

end

subroutine test_int_erf_bielec_ijkl
 implicit none
 integer :: i,j,k,l
 double precision :: accu_relative,accu_absolute
 double precision, allocatable :: integrals_array(:,:)
 double precision :: integral_erf,erf_mu_of_r_ao,test,ao_bielec_integral_erf,test_2,erf_mu_of_r_ao_old
 allocate(integrals_array(ao_num,ao_num))
 accu_relative = 0.d0
 accu_absolute = 0.d0
 i_count = 0
 
 integral_erf = erf_mu_of_r_ao(1,1,1,1)
 print*,''
 call wall_time(wall1)
 provide ao_bielec_integrals_erf_in_map ao_bielec_integrals_erf_mu_of_r_in_map 
 do i = 1, ao_num ! electron 1 
  do j = 1, ao_num ! electron 1 
   print*,i,j
!  call give_all_erf_mu_of_r_kl(i,j,integrals_array)
   do k = 1, ao_num ! electron 2 
    do l = 1, ao_num ! electron 2 
     double precision :: get_ao_bielec_integral_erf_mu_of_r,integral_erf_map,i_count
     
     integral_erf = ao_bielec_integral_erf(i,j,k,l)
    !integral_erf = erf_mu_of_r_ao(i,j,k,l)
  ! !test = erf_mu_of_r_ao(i,j,k,l)
!    test = integrals_array(l,k)
     test = get_ao_bielec_integral_erf_mu_of_r(i,k,j,l,ao_integrals_erf_mu_of_r_map)
     if(dabs(integral_erf + test).gt.0.d0)then
      i_count += 1.d0
      accu_absolute += dabs(integral_erf - test)!/dabs(integral_erf + test)*0.5d0
      accu_relative += dabs(integral_erf - test)/dabs(integral_erf + test)*0.5d0
     endif
     if(dabs(integral_erf - test).gt.1.d-10)then
      print*,'AHAHAH'
      print*,i,j,k,l
      print*,integral_erf,test
      print*,dabs(integral_erf - test),dabs(integral_erf - test)/dabs(integral_erf + test)*0.5d0
     endif
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall2)
 double precision :: wall1, wall2
 print*,'time = ',wall2-wall1
 accu_absolute = accu_absolute/(i_count)
 accu_relative = accu_relative/(i_count)
 print*,'average absolute error = ',accu_absolute 
 print*,'average relative error = ',accu_relative 
end

subroutine test_int_erf_bielec_ijkl_mo
 implicit none
 integer :: i,j,k,l
 integer :: m,n,p,q
 double precision :: accu,get_mo_bielec_integral_erf_mu_of_r, integral, get_ao_bielec_integral_erf_mu_of_r,integral_2

 do i = 1, mo_tot_num ! 1 
  do j = 1, mo_tot_num ! 2 
   do k = 1, mo_tot_num ! 1 
    do l = 1, mo_tot_num ! 2
     integral = get_mo_bielec_integral_erf_mu_of_r(i,j,k,l,mo_integrals_erf_mu_of_r_map)
     integral_2 = 0.d0
     do m = 1, ao_num
      do n = 1, ao_num
       do p = 1, ao_num
        do q = 1, ao_num   !                                                                          1             2                1              2
         integral_2 += get_ao_bielec_integral_erf_mu_of_r(m,n,p,q,ao_integrals_erf_mu_of_r_map) * mo_coef(m,i) * mo_coef(n,j)  * mo_coef(p,k) * mo_coef(q,l)  
        enddo
       enddo
      enddo
     enddo
      print*,i,j,k,l
      print*,integral,integral_2
     if(dabs(integral - integral_2).gt.1.d-6)then
      print*,'AHAHAH'
      print*,i,j,k,l
      print*,integral,integral_2
      print*,dabs(integral- integral_2),dabs(integral- integral_2)/dabs(integral+integral_2)*0.5d0
     endif
    enddo
   enddo
  enddo 
 enddo

end
