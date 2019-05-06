subroutine n_order_partial_derivative_ao_polynomial_part(r,ao_num_local,order_der,pol_part_derivative)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
!            : ao_num_local : number of the considered ao 
!
!            : order_der : order of the partial derivative for: order_der(1) ==> x, order_der(2) ==> y and order_der(3) ==> z respectively
!
! output     : pol_part_derivative = d^(n_x)/dx^(n_x)*d^(n_y)/dy^(n_y)*d^(n_z)/dz^(n_z)[pol_ao(ao_num_local)] at r, where n_x,n_y,n_z are respectively order_der(1),order_der(2),order_der(3)
!
 END_DOC
 double precision, intent(in) :: r(3)
 integer         , intent(in) :: ao_num_local
 integer         , intent(in) :: order_der(3)
 double precision, intent(out) :: pol_part_derivative 

 integer :: power_ao(3)
 integer :: num_ao
 double precision :: center_ao(3)
 double precision :: dxn,dyn,dzn
 double precision :: fact

 num_ao = ao_nucl(ao_num_local)
 center_ao(1:3) = nucl_coord(num_ao,1:3) 
 power_ao(1:3) = ao_power(ao_num_local,1:3)

 if(power_ao(1) .ge. order_der(1))then
  dxn = (fact(power_ao(1))/fact(power_ao(1)-order_der(1)))*(r(1) -center_ao(1))**(power_ao(1)-order_der(1))
 else
  dxn = 0.d0
 endif
 
 if(power_ao(2) .ge. order_der(2))then
  dyn = (fact(power_ao(2))/fact(power_ao(2)-order_der(2)))*(r(2) -center_ao(2))**(power_ao(2)-order_der(2))
 else
  dyn = 0.d0
 endif

 if(power_ao(3) .ge. order_der(3))then
  dzn = (fact(power_ao(3))/fact(power_ao(3)-order_der(3)))*(r(3) -center_ao(3))**(power_ao(3)-order_der(3))
 else
  dzn = 0.d0
 endif

 pol_part_derivative = dxn* dyn *dzn




end


subroutine n_order_partial_derivative_ao_exponential_part(r,ao_num_local,order_der,exp_part_derivative)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
!            : ao_num_local : number of the considered ao 
!
!            : order_der : order of the partial derivative for: order_der(1) ==> x, order_der(2) ==> y and order_der(3) ==> z respectively
!
! output     : exp_part_derivative = d^(n_x)/dx^(n_x)*d^(n_y)/dy^(n_y)*d^(n_z)/dz^(n_z)[exp_ao(ao_num_local)] at r, where n_x,n_y,n_z are respectively order_der(1),order_der(2),order_der(3)
!
 END_DOC
 double precision, intent(in) :: r(3)
 integer         , intent(in) :: ao_num_local
 integer         , intent(in) :: order_der(3)
 double precision, intent(out) :: exp_part_derivative 

 integer :: a,b,c,m,i,k,l
 integer :: num_ao
 double precision :: center_ao(3)
 double precision :: binom_func 
 double precision :: dx,dy,dz,r2
 double precision :: beta 
 double precision, allocatable :: d_h(:,:,:),d_x_f(:),d_y_f(:),d_z_f(:)
 allocate(d_h(0:order_der(1),0:order_der(2),0:order_der(3)),d_x_f(0:order_der(1)+2),d_y_f(0:order_der(2)+2),d_z_f(0:order_der(3)+2))

 num_ao = ao_nucl(ao_num_local)
 center_ao(1:3) = nucl_coord(num_ao,1:3) 

 dx = (r(1) - center_ao(1))
 dy = (r(2) - center_ao(2))
 dz = (r(3) - center_ao(3))
 r2 = dx*dx + dy*dy + dz*dz

 exp_part_derivative = 0.d0

!print*,' order_der_loc_dasn la sub =',order_der

 do m=1,ao_prim_num(ao_num_local)
  beta = ao_expo_ordered_transp(m,ao_num_local)
  d_h = 0.d0
  d_h(0,0,0)= dexp(-beta*r2)

  d_x_f = 0.d0
  d_y_f = 0.d0 
  d_z_f = 0.d0

  d_x_f(1) = (-2.d0 *beta*dx)
  d_y_f(1) = (-2.d0 *beta*dy)
  d_z_f(1) = (-2.d0 *beta*dz)

  d_x_f(2) = (-2.d0 *beta)
  d_y_f(2) = (-2.d0 *beta)
  d_z_f(2) = (-2.d0 *beta)

  if(order_der(1) .gt. 0 .and. order_der(2) .gt. 0 .and. order_der(3) .gt. 0)then
  !d_h(2,1,1)= (-2.d0 *beta*dx) * dexp(-beta*r2)
  !d_h(1,2,1)= (-2.d0 *beta*dy) * dexp(-beta*r2)
  !d_h(1,1,2)= (-2.d0 *beta*dz) * dexp(-beta*r2)
  !d_h(1,2,2)= (-2.d0 *beta*dy)*(-2.d0 *beta*dz) * dexp(-beta*r2)
  !d_h(2,2,1)= (-2.d0 *beta*dx)*(-2.d0 *beta*dy) * dexp(-beta*r2)
  !d_h(2,1,2)= (-2.d0 *beta*dx)*(-2.d0 *beta*dz) * dexp(-beta*r2)
  !d_h(2,2,2)= (-2.d0 *beta*dx)* (-2.d0 *beta*dy) * (-2.d0 *beta*dz) * dexp(-beta*r2)
   do c=1,order_der(3)
    do b=1,order_der(2)
     do a=1,order_der(1)
      do i=0,c-1
       do l=0,b-1
        do k=0,a-1
         d_h(a,b,c) += binom_func(a-1,k)*binom_func(b-1,l)*binom_func(c-1,i)*d_x_f(k+1)*d_y_f(l+1)*d_z_f(i+1)*d_h(a-1-k,b-1-l,c-1-i)
        enddo
       enddo
      enddo
     enddo 
    enddo
   enddo
  else if (order_der(1) .gt. 0 .and. order_der(2) .gt. 0 .and. order_der(3) .eq. 0)then
   d_h(1,0,0)= d_x_f(1) * d_h(0,0,0)
   d_h(0,1,0)= d_y_f(1) * d_h(0,0,0)
  !d_h(1,1,0)= d_x_f(1) * d_y_f(1) * d_h(0,0,0) 
  !d_h(2,1,0)= d_y_f(1)*(d_x_f(1)*d_h(1,0,0)+ d_x_f(2)*d_h(0,0,0)) 
  !d_h(2,1,0)= d_y_f(1)*((d_x_f(1)**2.d0)*d_h(0,0,0)+ d_x_f(2)*d_h(0,0,0)) 
   do b=1,order_der(2)
    do a=1,order_der(1)
     do l=0,b-1
      do k=0,a-1
       d_h(a,b,0) += binom_func(a-1,k)*binom_func(b-1,l)*d_x_f(k+1)*d_y_f(l+1)*d_h(a-1-k,b-1-l,0)
       !print*,'d_h(',a,b,0,')   =', d_h(a,b,0)
      enddo
     enddo
    enddo 
   enddo
  else if (order_der(1) .gt. 0 .and. order_der(2) .eq. 0 .and. order_der(3) .gt. 0)then
   d_h(1,0,0)= d_x_f(1) * d_h(0,0,0)
   d_h(0,0,1)= d_z_f(1) * d_h(0,0,0)
  !d_h(2,1,1)= (-2.d0 *beta*dx) * dexp(-beta*r2)
  !d_h(1,1,2)= (-2.d0 *beta*dz) * dexp(-beta*r2)
  !d_h(2,1,2)= (-2.d0 *beta*dx)*(-2.d0 *beta*dz) * dexp(-beta*r2)
   do c=1,order_der(3)
    do a=1,order_der(1)
     do i=0,c-1
      do k=0,a-1
       d_h(a,0,c) += binom_func(a-1,k)*binom_func(c-1,i)*d_x_f(k+1)*d_z_f(i+1)*d_h(a-1-k,0,c-1-i)
      enddo
     enddo
    enddo 
   enddo
  else if (order_der(1) .eq. 0 .and. order_der(2) .gt. 0 .and. order_der(3) .gt.0) then
   d_h(0,1,0)= d_y_f(1) * d_h(0,0,0)
   d_h(0,0,1)= d_z_f(1) * d_h(0,0,0)
  !d_h(1,2,1)= (-2.d0 *beta*dy) * dexp(-beta*r2)
  !d_h(1,1,2)= (-2.d0 *beta*dz) * dexp(-beta*r2)
  !d_h(1,2,2)= (-2.d0 *beta*dy)*(-2.d0 *beta*dz) * dexp(-beta*r2)
   do c=1,order_der(3)
    do b=1,order_der(2)
     do i=0,c-1
      do l=0,b-1
       d_h(0,b,c) += binom_func(b-1,l)*binom_func(c-1,i)*d_y_f(l+1)*d_z_f(i+1)*d_h(0,b-1-l,c-1-i)
      enddo
     enddo
    enddo
   enddo
  else if (order_der(1) .gt. 0 .and. order_der(2) .eq. 0 .and. order_der(3) .eq. 0)then
  !d_h(2,1,1)= (-2.d0 *beta*dx) * dexp(-beta*r2)
 ! print*,'je passe ici'
   do a=1,order_der(1)
    do k=0,a-1
     d_h(a,0,0) += binom_func(a-1,k)*d_x_f(k+1)*d_h(a-1-k,0,0)
    enddo
   enddo 
  else if (order_der(1) .eq. 0 .and. order_der(2) .gt. 0 .and. order_der(3) .eq. 0)then
  !d_h(1,2,1)= (-2.d0 *beta*dy) * dexp(-beta*r2)
   do b=1,order_der(2)
    do l=0,b-1
     d_h(0,b,0) += binom_func(b-1,l)*d_y_f(l+1)*d_h(0,b-1-l,0)
    enddo
   enddo
  else if (order_der(1) .eq. 0 .and. order_der(2) .eq. 0 .and. order_der(3) .gt. 0)then
  !d_h(1,1,2)= (-2.d0 *beta*dz) * dexp(-beta*r2)
   do c=1,order_der(3)
    do i=0,c-1
     d_h(0,0,c) += binom_func(c-1,i)*d_z_f(i+1)*d_h(0,0,c-1-i)
    enddo
   enddo
  endif
 
  exp_part_derivative += ao_coef_normalized_ordered_transp(m,ao_num_local) * d_h(order_der(1),order_der(2),order_der(3)) 
 !print*,'exp_part_derivative  =',exp_part_derivative


 enddo

end


subroutine n_order_partial_derivative_ao(r,ao_num_local,order_der,ao_part_derivative)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
!            : ao_num_local : number of the considered ao 
!
!            : order_der : order of the partial derivative for: order_der(1) ==> x, order_der(2) ==> y and order_der(3) ==> z respectively
!
! output     : exp_part_derivative = d^(n_x)/dx^(n_x)*d^(n_y)/dy^(n_y)*d^(n_z)/dz^(n_z)[ao(ao_num_local)] at r, where n_x,n_y,n_z are respectively order_der(1),order_der(2),order_der(3)
!
 END_DOC
 double precision, intent(in) :: r(3)
 integer         , intent(in) :: ao_num_local
 integer         , intent(in) :: order_der(3)
 double precision, intent(out) :: ao_part_derivative 

 integer :: k,l,m
 double precision :: binom_func
 double precision :: pol_part_derivative,exp_part_derivative
 integer :: order_der_loc_pol(3)
 integer :: order_der_loc_exp(3)

 ao_part_derivative = 0.d0

 do k=0, order_der(1)
  do l=0, order_der(2)
   do m=0, order_der(3)

    order_der_loc_pol(1)= order_der(1)-k
    order_der_loc_pol(2)= order_der(2)-l
    order_der_loc_pol(3)= order_der(3)-m

    order_der_loc_exp(1)= k 
    order_der_loc_exp(2)= l
    order_der_loc_exp(3)= m

   !print*,' order_der_loc =',order_der
   !print*,' order_der_loc_exp =',order_der_loc_exp
     
    call n_order_partial_derivative_ao_polynomial_part(r,ao_num_local,order_der_loc_pol,pol_part_derivative)
    call n_order_partial_derivative_ao_exponential_part(r,ao_num_local,order_der_loc_exp,exp_part_derivative) 
   !print*,'*******'
   !print*,'order_der_loc_pol',order_der_loc_pol(1),order_der_loc_pol(2),order_der_loc_pol(3)
   !print*,'order_der_loc_exp',order_der_loc_exp(1),order_der_loc_exp(2),order_der_loc_exp(3)
   !print*,'pol_part_derivative',pol_part_derivative
   !print*,'exp_part_derivative',exp_part_derivative
   !print*,'binome tot',binom_func(order_der(3),m)* binom_func(order_der(1),k) * binom_func(order_der(2),l) 
    ao_part_derivative += binom_func(order_der(1),k) * binom_func(order_der(2),l) * binom_func(order_der(3),m) * pol_part_derivative * exp_part_derivative
   !print*,'ao_part_derivative  =',ao_part_derivative
   enddo
  enddo
 enddo

end


subroutine n_order_partial_derivative_ao_product(r,ao_num_local_1,ao_num_local_2,order_der,ao_part_derivative)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
!            : ao_num_local_1 : number of the first ao 
!
!            : ao_num_local_2 : number of the second ao 
!
!            : order_der : order of the partial derivative for: order_der(1) ==> x, order_der(2) ==> y and order_der(3) ==> z respectively
!
! output     : exp_part_derivative = d^(n_x)/dx^(n_x)*d^(n_y)/dy^(n_y)*d^(n_z)/dz^(n_z)[ao(ao_num_local_1)ao(ao_num_local_2)] at r, where n_x,n_y,n_z are respectively order_der(1),order_der(2),order_der(3)
!
 END_DOC
 double precision, intent(in) :: r(3)
 integer         , intent(in) :: ao_num_local_1
 integer         , intent(in) :: ao_num_local_2
 integer         , intent(in) :: order_der(3)
 double precision, intent(out) :: ao_part_derivative 

 integer :: k,l,m
 double precision :: binom_func
 double precision :: ao_part_derivative_1,ao_part_derivative_2 
 integer :: order_der_loc_1(3)
 integer :: order_der_loc_2(3)

 ao_part_derivative = 0.d0

 do k=0, order_der(1)
  do l=0, order_der(2)
   do m=0, order_der(3)

    order_der_loc_1(1)= order_der(1)-k
    order_der_loc_1(2)= order_der(2)-l
    order_der_loc_1(3)= order_der(3)-m

    order_der_loc_2(1)= k 
    order_der_loc_2(2)= l
    order_der_loc_2(3)= m

   !print*,' order_der_loc =',order_der
   !print*,' order_der_loc_exp =',order_der_loc_exp
 
    call n_order_partial_derivative_ao(r,ao_num_local_1,order_der_loc_1,ao_part_derivative_1) 
    call n_order_partial_derivative_ao(r,ao_num_local_2,order_der_loc_2,ao_part_derivative_2) 
   !print*,'*******'
   !print*,'order_der_loc_pol',order_der_loc_pol(1),order_der_loc_pol(2),order_der_loc_pol(3)
   !print*,'order_der_loc_exp',order_der_loc_exp(1),order_der_loc_exp(2),order_der_loc_exp(3)
   !print*,'pol_part_derivative',pol_part_derivative
   !print*,'exp_part_derivative',exp_part_derivative
   !print*,'binome tot',binom_func(order_der(3),m)* binom_func(order_der(1),k) * binom_func(order_der(2),l) 
    ao_part_derivative += binom_func(order_der(1),k) * binom_func(order_der(2),l) * binom_func(order_der(3),m) * ao_part_derivative_1* ao_part_derivative_2 
   !print*,'ao_part_derivative  =',ao_part_derivative
   enddo
  enddo
 enddo

end




subroutine n_order_delta_ao_product(r,ao_num_local_1,ao_num_local_2,n_delta,n_order_delta)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
!            : ao_num_local_1 : number of the first ao 
!
!            : ao_num_local_2 : number of the second ao 
!
!            : n_delta : number of application of the laplacian operator on the ao product 
!
! output     : n_order_delta = Delta^{n_delta}[ao(ao_num_local_1)ao(ao_num_local_2)] at r, where Delta is the laplacian operator
!
 END_DOC
 double precision, intent(in) :: r(3)
 integer         , intent(in) :: ao_num_local_1
 integer         , intent(in) :: ao_num_local_2
 integer         , intent(in) :: n_delta 
 double precision, intent(out) :: n_order_delta 

 integer :: k,l,m
 double precision :: binom_func
 double precision :: ao_part_derivative
 integer :: order_der(3)

 n_order_delta = 0.d0

 do k=0, n_delta
  do l=0, n_delta - k

    order_der(1)= 2*k 
    order_der(2)= 2*l
    order_der(3)= 2*(n_delta - k - l)

    call n_order_partial_derivative_ao_product(r,ao_num_local_1,ao_num_local_2,order_der,ao_part_derivative) 

    n_order_delta += binom_func(n_delta,k) * binom_func(n_delta-k,l) * ao_part_derivative
  enddo
 enddo

end



 subroutine give_delta_n_at_r(r,n_delta,delta_n_at_r)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! input      : n_delta => number of application of the laplacian operator on the ao product 
!
! output     : delta_n_at_r ==> Delta^{n_delta}[ao(ao_num_local_1)ao(ao_num_local_2)] at r, where Delta is the laplacian operator 
!
!            : 
 END_DOC
 double precision, intent(in) :: r(3)
 integer, intent(in) :: n_delta
 double precision, intent(out) :: delta_n_at_r(ao_num,ao_num) 

 integer :: i,j
 double precision :: n_order_delta 
 
 do i = 1,ao_num
  do j = 1,ao_num

   call n_order_delta_ao_product(r,j,i,n_delta,n_order_delta)
   delta_n_at_r(j,i) = n_order_delta 

  enddo
 enddo

 end


 BEGIN_PROVIDER[double precision, n_sphe_av_exp,(n_points_final_grid,0:order_derivative_ontop,n_states)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: k,l,n,istate,r
 double precision :: n_order_delta,r1(3)

 n_sphe_av_exp = 0.d0 

 do istate = 1, n_states
  do n = 0, order_derivative_ontop      
   do r = 1, n_points_final_grid
    r1(1) = final_grid_points(1,r)
    r1(2) = final_grid_points(2,r)
    r1(3) = final_grid_points(3,r)
    do l = 1, ao_num       
     do k = 1, ao_num        
      call n_order_delta_ao_product(r1,l,k,n,n_order_delta)
      n_sphe_av_exp(r,n,istate) += (1.d0/(2.d0*dble(n)+1.d0))*((one_e_dm_alpha_ao_for_dft(k,l,istate)+one_e_dm_beta_ao_for_dft(k,l,istate)) * n_order_delta)
     enddo
    enddo
   enddo
  enddo
 enddo  


 END_PROVIDER
