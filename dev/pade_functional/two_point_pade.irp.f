
 BEGIN_PROVIDER[double precision, two_point_pade_h_a_3,(0:3,n_states)]
&BEGIN_PROVIDER[double precision, two_point_pade_h_b_3,(1:3,n_states)]
 implicit none
 integer :: istate,i
 double precision :: c1(0:3,n_states),c2(0:2,n_states)

 c2 = 0.d0
 c1 = 0.d0
 two_point_pade_h_a_3 = 0.d0
 two_point_pade_h_b_3 = 0.d0

 do istate = 1, n_states

  do i = 0,1
   c2(2*i,istate)=  (2.d0 *dble(i)+1.d0) * sr_hatree_energy_small_mu_coef(i,istate) 
  enddo

  do i = 0,1
   c1(2*i+1,istate)= -(2.d0 *dble(i)+2.d0) * coeff_H_exp(i,istate) 
  enddo

  do i = 0,2
   print*,"c2(",i,") = ", c2(i,istate)
  enddo
  print*,"************************"

  do i = 0,3 
   print*,"c1(",i,") = ", c1(i,istate)
  enddo
  
  two_point_pade_h_a_3(0,istate)=c2(0,istate) 
  
  two_point_pade_h_b_3(2,istate)=-c2(2,istate)/c2(0,istate)
  two_point_pade_h_b_3(3,istate)=c2(0,istate)/c1(3,istate)

 enddo
 END_PROVIDER


 subroutine give_two_p_pade_dh_3(mu,two_pade_dh_3)
 implicit none
 double precision, intent(in)  :: mu
 double precision, intent(out) :: two_pade_dh_3(n_states)
 double precision :: a(0:3,n_states),b(1:3,n_states)
 double precision :: num,denom
 integer :: istate,i
 a=two_point_pade_h_a_3
 b=two_point_pade_h_b_3 

 do istate = 1, n_states
  num = 0.d0
  denom= 0.d0

  do i= 0,3
   num += a(i,istate)*mu**(dble(i))
  enddo

  do i=1,3
   denom += b(i,istate)*mu**(dble(i))
  enddo

 two_pade_dh_3(istate) = num/(1+denom)
 enddo  
 end

 BEGIN_PROVIDER[double precision, two_point_pade_h_a_4,(0:4,n_states)]
&BEGIN_PROVIDER[double precision, two_point_pade_h_b_4,(1:4,n_states)]
 implicit none
 integer :: istate,i
 double precision :: c1(0:4,n_states),c2(0:3,n_states)

 c2 = 0.d0
 c1 = 0.d0
 two_point_pade_h_a_4 = 0.d0
 two_point_pade_h_b_4 = 0.d0

 do istate = 1, n_states

  do i = 0,1
   c2(2*i,istate)=  (2.d0 *dble(i)+1.d0) * sr_hatree_energy_small_mu_coef(i,istate) 
  enddo

  do i = 0,1
   c1(2*i+1,istate)= -(2.d0 *dble(i)+2.d0) * coeff_H_exp(i,istate) 
  enddo

  do i = 0,3
   print*,"c2(",i,") = ", c2(i,istate)
  enddo
  print*,"************************"

  do i = 0,4
   print*,"c1(",i,") = ", c1(i,istate)
  enddo
  
  two_point_pade_h_a_4(0,istate)=c2(0,istate) 
  two_point_pade_h_a_4(1,istate)=-c2(0,istate)**(3.d0)/(c1(3,istate)*c2(2,istate)) 
  
  two_point_pade_h_b_4(1,istate)=-c2(0,istate)**(2.d0)/(c1(3,istate)*c2(2,istate))
  two_point_pade_h_b_4(2,istate)=-c2(2,istate)/c2(0,istate)
  two_point_pade_h_b_4(3,istate)=c2(0,istate)/c1(3,istate)
  two_point_pade_h_b_4(4,istate)=-c2(0,istate)**(3.d0)/(c1(3,istate)**(2.d0)*c2(2,istate))

 enddo
 END_PROVIDER


 subroutine give_two_p_pade_dh_4(mu,two_pade_dh_4)
 implicit none
 double precision, intent(in)  :: mu
 double precision, intent(out) :: two_pade_dh_4(n_states)
 double precision :: a(0:4,n_states),b(1:4,n_states)
 double precision :: num,denom
 integer :: istate,i
 a=two_point_pade_h_a_4
 b=two_point_pade_h_b_4 

 do istate = 1, n_states
  num = 0.d0
  denom= 0.d0

  do i= 0,4
   num += a(i,istate)*mu**(dble(i))
  enddo

  do i=1,4
   denom += b(i,istate)*mu**(dble(i))
  enddo

 two_pade_dh_4(istate) = num/(1+denom)
 enddo  
 end


 BEGIN_PROVIDER[double precision, two_point_pade_h_a_5,(0:5,n_states)]
&BEGIN_PROVIDER[double precision, two_point_pade_h_b_5,(1:5,n_states)]
 implicit none
 integer :: istate,i
 double precision :: c1(0:5,n_states),c2(0:4,n_states)
 double precision :: denom

 c2 = 0.d0
 c1 = 0.d0
 two_point_pade_h_a_5 = 0.d0
 two_point_pade_h_b_5 = 0.d0

 do istate = 1, n_states

  do i = 0,2
   c2(2*i,istate)=  (2.d0 *dble(i)+1.d0) * sr_hatree_energy_small_mu_coef(i,istate) 
  enddo

  do i = 0,2
   c1(2*i+1,istate)= -(2.d0 *dble(i)+2.d0) * coeff_H_exp(i,istate) 
  enddo

  print*,"************************"
  do i = 0,4
   print*,"c2(",i,") = ", c2(i,istate)
  enddo
  print*,"************************"
  do i = 0,5
   print*,"c1(",i,") = ", c1(i,istate)
  enddo
  
  denom=c1(5,istate)*(c2(0,istate)**4.d0)+(c1(3,istate)**3.d0)*c2(2,istate)**2.d0


  two_point_pade_h_a_5(0,istate)=c2(0,istate) 
  two_point_pade_h_a_5(1,istate)=-((c1(3,istate)*c2(0,istate)**2.d0*(c1(3,istate)*c2(0,istate)*c2(2,istate)-c1(5,istate)*c2(2,istate)**2.d0+c1(5,istate)*c2(0,istate)*c2(4,istate)))/(denom))

  two_point_pade_h_a_5(2,istate)=((c1(3,istate)*c2(0,istate)**5.d0+c1(3,istate)**3.d0*c2(2,istate)**3.d0-c1(3,istate)**3.d0*c2(0,istate)*c2(2,istate)*c2(4,istate))/(denom))


  two_point_pade_h_b_5(1,istate)=-((c1(3,istate)*c2(0,istate)*(c1(3,istate)*c2(0,istate)*c2(2,istate)-c1(5,istate)*c2(2,istate)**2.d0+c1(5,istate)*c2(0,istate)*c2(4,istate)))/(denom))

  two_point_pade_h_b_5(2,istate)= ((c1(3,istate)*c2(0,istate)**4.d0-c1(5,istate)*c2(0,istate)**3.d0*c2(2,istate)-c1(3,istate)**3.d0*c2(2,istate)*c2(4,istate))/(denom))

  two_point_pade_h_b_5(3,istate)=((c1(3,istate)*c2(2,istate)*(c1(3,istate)*c2(0,istate)*c2(2,istate)-c1(5,istate)*c2(2,istate)**2.d0+c1(5,istate)*c2(0,istate)*c2(4,istate)))/(denom))

  two_point_pade_h_b_5(4,istate)=((-c1(3,istate)*c2(0,istate)**3.d0*c2(2,istate)+c1(5,istate)*c2(0,istate)**2.d0*c2(2,istate)**2.d0-c1(5,istate)*c2(0,istate)**3.d0*c2(4,istate))/(denom))

  two_point_pade_h_b_5(5,istate)=((c2(0,istate)**5.d0+c1(3,istate)**2.d0*c2(2,istate)**3.d0-c1(3,istate)**2.d0*c2(0,istate)*c2(2,istate)*c2(4,istate))/(denom))

 enddo
 END_PROVIDER

 subroutine give_two_p_pade_dh_5(mu,two_pade_dh_5)
 implicit none
 double precision, intent(in)  :: mu
 double precision, intent(out) :: two_pade_dh_5(n_states)
 double precision :: a(0:5,n_states),b(1:5,n_states)
 double precision :: num,denom
 integer :: istate,i
 a=two_point_pade_h_a_5
 b=two_point_pade_h_b_5

 do istate = 1, n_states
  num = 0.d0
  denom= 0.d0

  do i= 0,5
   num += a(i,istate)*mu**(dble(i))
  enddo

  do i=1,5
   denom += b(i,istate)*mu**(dble(i))
  enddo

 two_pade_dh_5(istate) = num/(1+denom)
 enddo  
 end


 BEGIN_PROVIDER[double precision, two_point_pade_h_a_6,(0:6,n_states)]
&BEGIN_PROVIDER[double precision, two_point_pade_h_b_6,(1:6,n_states)]
 implicit none
 integer :: istate,i
 double precision :: c1(0:6,n_states),c2(0:5,n_states)
 double precision :: denom

 c2 = 0.d0
 c1 = 0.d0
 two_point_pade_h_a_6 = 0.d0
 two_point_pade_h_b_6 = 0.d0

 do istate = 1, n_states

  do i = 0,2
   c2(2*i,istate)=  (2.d0 *dble(i)+1.d0) * sr_hatree_energy_small_mu_coef(i,istate) 
  enddo

  do i = 0,2
   c1(2*i+1,istate)= -(2.d0 *dble(i)+2.d0) * coeff_H_exp(i,istate) 
  enddo

  print*,"************************"
  do i = 0,5
   print*,"c2(",i,") = ", c2(i,istate)
  enddo
  print*,"************************"
  do i = 0,6
   print*,"c1(",i,") = ", c1(i,istate)
  enddo
  
  denom= c1(3,istate)**2.d0*c2(0,istate)**4.d0-2.d0*c1(3,istate)*c1(5,istate)*c2(0,istate)**3.d0*c2(2,istate)+c1(5,istate)**2.d0*c2(0,istate)**2.d0*c2(2,istate)**2.d0-c1(5,istate)**2.d0*c2(0,istate)**3.d0*c2(4,istate)-c1(3,istate)**4.d0*c2(2,istate)*c2(4,istate)


  two_point_pade_h_a_6(0,istate)=c2(0,istate) 

  two_point_pade_h_a_6(1,istate)=((c2(0,istate)*(c1(5,istate)*c2(0,istate)**5.d0+2.d0*c1(3,istate)**3.d0*c2(0,istate)*c2(2,istate)**2.d0-c1(3,istate)**2.d0*c1(5,istate)*c2(2,istate)**3.d0 - c1(3,istate)**3.d0*c2(0,istate)**2.d0*c2(4,istate)+c1(3,istate)**2.d0*c1(5,istate)*c2(0,istate)*c2(2,istate)*c2(4,istate)))/(denom))

  two_point_pade_h_a_6(2,istate)=((c1(3,istate)**2.d0*c2(0,istate)**4.d0*c2(2,istate)-c1(3,istate)*c1(5,istate)*c2(0,istate)**3.d0*c2(2,istate)**2.d0+c1(3,istate)*c1(5,istate)*c2(0,istate)**4.d0*c2(4,istate)+c1(3,istate)**4.d0*c2(2,istate)**2.d0*c2(4,istate)-c1(3,istate)**4.d0*c2(0,istate)*c2(4,istate)**2.d0)/(-denom))

  two_point_pade_h_a_6(3,istate)=((c1(3,istate)*(-c2(0,istate)**6.d0-2.d0*c1(3,istate)**2.d0*c2(0,istate)*c2(2,istate)**3.d0+c1(3,istate)*c1(5,istate)*c2(2,istate)**4.d0+2.d0*c1(3,istate)**2.d0*c2(0,istate)**2.d0*c2(2,istate)*c2(4,istate)-2.d0*c1(3,istate)*c1(5,istate)*c2(0,istate)*c2(2,istate)**2.d0*c2(4,istate)+c1(3,istate)*c1(5,istate)*c2(0,istate)**2.d0*c2(4,istate)**2.d0))/(-denom))


  two_point_pade_h_b_6(1,istate)= (( -c1(5,istate)*c2(0,istate)**5.d0-2.d0*c1(3,istate)**3.d0*c2(0,istate)*c2(2,istate)**2.d0+c1(3,istate)**2.d0*c1(5,istate)*c2(2,istate)**3.d0+c1(3,istate)**3.d0*c2(0,istate)**2.d0*c2(4,istate)-c1(3,istate)**2.d0*c1(5,istate)*c2(0,istate)*c2(2,istate)*c2(4,istate))/(-denom))

  two_point_pade_h_b_6(2,istate)=((2.d0*c1(3,istate)**2.d0*c2(0,istate)**3.d0*c2(2,istate)-3.d0*c1(3,istate)*c1(5,istate)*c2(0,istate)**2.d0*c2(2,istate)**2.d0+c1(5,istate)**2.d0*c2(0,istate)*c2(2,istate)**3.d0+c1(3,istate)*c1(5,istate)*c2(0,istate)**3.d0*c2(4,istate)-c1(5,istate)**2.d0*c2(0,istate)**2.d0*c2(2,istate)*c2(4,istate)-c1(3,istate)**4.d0*c2(4,istate)**2.d0)/(-denom))

  two_point_pade_h_b_6(3,istate)=((-c1(3,istate)*c2(0,istate)**5.d0+c1(5,istate)*c2(0,istate)**4.d0*c2(2,istate)+c1(3,istate)**3.d0*c2(0,istate)*c2(2,istate)*c2(4,istate)-c1(3,istate)**2.d0*c1(5,istate)*c2(2,istate)**2.d0*c2(4,istate)+c1(3,istate)**2.d0*c1(5,istate)*c2(0,istate)*c2(4,istate)**2.d0)/(-denom))

  two_point_pade_h_b_6(4,istate)= ((-2.d0*c1(3,istate)**2.d0*c2(0,istate)**2.d0*c2(2,istate)**2.d0+3.d0*c1(3,istate)*c1(5,istate)*c2(0,istate)*c2(2,istate)**3.d0-c1(5,istate)**2.d0*c2(2,istate)**4.d0+c1(3,istate)**2.d0*c2(0,istate)**3.d0*c2(4,istate)-3.d0*c1(3,istate)*c1(5,istate)*c2(0,istate)**2.d0*c2(2,istate)*c2(4,istate)+2.d0*c1(5,istate)**2.d0*c2(0,istate)*c2(2,istate)**2.d0*c2(4,istate)-c1(5,istate)**2.d0*c2(0,istate)**2.d0*c2(2,istate)**2.d0)/(-denom))

  two_point_pade_h_b_6(5,istate)= ((c1(3,istate)*c2(0,istate)**4.d0*c2(2,istate)-c1(5,istate)*c2(0,istate)**3.d0*c2(2,istate)**2.d0+c1(5,istate)*c2(0,istate)**4.d0*c2(4,istate)+c1(3,istate)**3.d0*c2(2,istate)**2.d0*c2(4,istate)-c1(3,istate)**3.d0*c2(0,istate)*c2(4,istate)**2.d0)/(-denom))

  two_point_pade_h_b_6(6,istate)=((-c2(0,istate)**6.d0-2.d0*c1(3,istate)**2.d0*c2(0,istate)*c2(2,istate)**3.d0+c1(3,istate)*c1(5,istate)*c2(2,istate)**4.d0+2.d0*c1(3,istate)**2.d0*c2(0,istate)**2.d0*c2(2,istate)*c2(4,istate)-2.d0*c1(3,istate)*c1(5,istate)*c2(0,istate)*c2(2,istate)**2.d0*c2(4,istate)+c1(3,istate)*c1(5,istate)*c2(0,istate)**2.d0*c2(4,istate)**2.d0)/(-denom))

 enddo
 END_PROVIDER

 subroutine give_two_p_pade_dh_6(mu,two_pade_dh_6)
 implicit none
 double precision, intent(in)  :: mu
 double precision, intent(out) :: two_pade_dh_6(n_states)
 double precision :: a(0:6,n_states),b(1:6,n_states)
 double precision :: num,denom
 integer :: istate,i
 a=two_point_pade_h_a_6
 b=two_point_pade_h_b_6

 do istate = 1, n_states
  num = 0.d0
  denom= 0.d0

  do i= 0,6
   num += a(i,istate)*mu**(dble(i))
  enddo

  do i=1,6
   denom += b(i,istate)*mu**(dble(i))
  enddo

 two_pade_dh_6(istate) = num/(1+denom)
 enddo  
 end


 BEGIN_PROVIDER[double precision, two_point_pade_h_a_7,(0:7,n_states)]
&BEGIN_PROVIDER[double precision, two_point_pade_h_b_7,(1:7,n_states)]
 implicit none
 integer :: istate,i
 double precision :: C1(0:7),C2(0:6)
 double precision :: denom

 C2 = 0.d0
 C1 = 0.d0
 two_point_pade_h_a_7 = 0.d0
 two_point_pade_h_b_7 = 0.d0

 do istate = 1, n_states

  do i = 0,3
   C2(2*i)=  (2.d0 *dble(i)+1.d0) * sr_hatree_energy_small_mu_coef(i,istate) 
  enddo

  do i = 0,3
   C1(2*i+1)= -(2.d0 *dble(i)+2.d0) * coeff_H_exp(i,istate) 
  enddo

  print*,"************************"
  do i = 0,6
   print*,"c2(",i,") = ", C2(i)
  enddo
  print*,"************************"
  do i = 0,7
   print*,"c1(",i,") = ", C1(i)
  enddo
  
 include 'pade_7.inc'

 enddo
 END_PROVIDER

 subroutine give_two_p_pade_dh_7(mu,two_pade_dh_7)
 implicit none
 double precision, intent(in)  :: mu
 double precision, intent(out) :: two_pade_dh_7(n_states)
 double precision :: a(0:7,n_states),b(1:7,n_states)
 double precision :: num,denom
 integer :: istate,i
 a=two_point_pade_h_a_7
 b=two_point_pade_h_b_7

 do istate = 1, n_states
  num = 0.d0
  denom= 0.d0

  do i= 0,7
   num += a(i,istate)*mu**(dble(i))
  enddo

  do i=1,7
   denom += b(i,istate)*mu**(dble(i))
  enddo

 two_pade_dh_7(istate) = num/(1+denom)
 enddo  
 end

 
 BEGIN_PROVIDER[double precision, simpson_int_test,(n_states)]
&BEGIN_PROVIDER[double precision, trapez_int_test,(n_states)]
 implicit none
 integer :: istate,i,n
 double precision :: a,b,h,exact,part1,part2,approx,part_trapez,trapez 
 b=10000.d0
 a = 0.5d0
 n= 10000000

 h = (b - a)/dble(n)
 part1=0.d0
 part2=0.d0

 do i=1,n/2-1
  part1 += (a+2.d0*dble(i)*h)**(-2.d0)
 enddo

 do i=1,n/2
  part2 += (a+(2.d0*dble(i)-1.d0)*h)**(-2.d0)
 enddo

 approx = (h/3.d0)*((a)**(-2.d0)+(b)**(-2.d0)+2.d0*part1 + 4.d0*part2 )

 exact= 1.d0/a - 1/b

 part_trapez = 0.d0
 do i=1,n-1
  part_trapez += (a+(dble(i))*h)**(-2.d0)
 enddo
 
 trapez= h*(((a)**(-2.d0)+(b)**(-2.d0))*0.5d0 + part_trapez)

 double precision :: integral, error
 double precision, allocatable :: listx(:), listy(:)
 allocate(listx(n),listy(n))

  do i=1,n
   listx(i) = a+dble(i)*h
   listy(i) = (a+(dble(i))*h)**(-2.d0) 
  enddo

  call DsplIntegr (listx, listy, n , a, b, integral, error)


 print*,"simpson        =",approx
 print*,"exact          =",exact 
 print*,"Trapez         =",trapez
 print*,"JULIEN         =",integral
 print*,"error simpson  =",approx-exact
 print*,"error Trapez   =",trapez-exact
 print*,"error Julien   =",integral -exact
 print*,"b-2     =",(b)**(-2.d0)


 END_PROVIDER



!BEGIN_PROVIDER[double precision, two_p_pade_h_3,(n_states)]
!BEGIN_PROVIDER[double precision, two_p_pade_h_4,(n_states)]
!BEGIN_PROVIDER[double precision, two_p_pade_h_5,(n_states)]
!BEGIN_PROVIDER[double precision, two_p_pade_h_6,(n_states)]
!BEGIN_PROVIDER[double precision, two_p_pade_h_7,(n_states)]
!implicit none
!integer :: istate,i,n
!double precision :: a,b,h,part1,part2,mu,mu_c,f_b,f_a,two_pade_dh_3(n_states)
!double precision :: part1_4,part2_4,f_b_4,f_a_4,two_pade_dh_4(n_states)
!double precision :: part1_5,part2_5,f_b_5,f_a_5,two_pade_dh_5(n_states)
!double precision :: part1_6,part2_6,f_b_6,f_a_6,two_pade_dh_6(n_states)
!double precision :: part1_7,part2_7,f_b_7,f_a_7,two_pade_dh_7(n_states)
!b=10000.d0
!a = mu_erf
!n= 10000000 

!b=1.d0
!a =sqrt(mu_erf**2.d0/(1.d0+mu_erf**2.d0))
!n= 10000000

!h = (b - a)/dble(n)

!do istate = 1, n_states
!!part1 = 0.d0
!!part2 = 0.d0

! part1_4 = 0.d0
! part2_4 = 0.d0

!!part1_5 = 0.d0
!!part2_5 = 0.d0

!!part1_6 = 0.d0
!!part2_6 = 0.d0

!!part1_7 = 0.d0
!!part2_7 = 0.d0

! do i=1,n/2-1
!  mu = a+(2.d0*dble(i))*h
!  mu_c =sqrt(mu**2.d0/(1+mu**2.d0))
!! call give_two_p_pade_dh_3(mu,two_pade_dh_3)
!  call give_two_p_pade_dh_4(mu_c,two_pade_dh_4)
!! call give_two_p_pade_dh_5(mu,two_pade_dh_5)
!! call give_two_p_pade_dh_6(mu,two_pade_dh_6)
!! call give_two_p_pade_dh_7(mu,two_pade_dh_7)
!! part1 += two_pade_dh_3(istate) 
!  part1_4 += two_pade_dh_4(istate)*(1.d0-mu_c**2.d0)**(-3.d0/2.d0) 
!! part1_5 += two_pade_dh_5(istate) 
!! part1_6 += two_pade_dh_6(istate) 
!! part1_7 += two_pade_dh_7(istate) 
! enddo
!
! do i=1,n/2
!  mu=a+(2.d0*dble(i)-1.d0)*h
!  mu_c =sqrt(mu**2.d0/(1+mu**2.d0))
!! call give_two_p_pade_dh_3(mu,two_pade_dh_3)
!  call give_two_p_pade_dh_4(mu_c,two_pade_dh_4)
!! call give_two_p_pade_dh_5(mu,two_pade_dh_5)
!! call give_two_p_pade_dh_6(mu,two_pade_dh_6)
!! call give_two_p_pade_dh_7(mu,two_pade_dh_7)
!! part2 += two_pade_dh_3(istate) 
!  part2_4 += two_pade_dh_4(istate)*(1.d0-mu_c**2.d0)**(-3.d0/2.d0)
!! part2_5 += two_pade_dh_5(istate) 
!! part2_6 += two_pade_dh_6(istate) 
!! part2_7 += two_pade_dh_7(istate) 
! enddo
! 
!  mu=a
!  mu_c =sqrt(mu**2.d0/(1+mu**2.d0))
!! call give_two_p_pade_dh_3(mu,two_pade_dh_3)
!  call give_two_p_pade_dh_4(mu_c,two_pade_dh_4)
!! call give_two_p_pade_dh_5(mu,two_pade_dh_5)
!! call give_two_p_pade_dh_6(mu,two_pade_dh_6)
!! call give_two_p_pade_dh_7(mu,two_pade_dh_7)
!! f_a=two_pade_dh_3(istate)
!  f_a_4=two_pade_dh_4(istate)*(1.d0-mu_c**2.d0)**(-3.d0/2.d0)
!! f_a_5=two_pade_dh_5(istate)
!! f_a_6=two_pade_dh_6(istate)
!! f_a_7=two_pade_dh_7(istate)

!  mu=b
!  mu_c =sqrt(mu**2.d0/(1+mu**2.d0))
!! call give_two_p_pade_dh_3(mu,two_pade_dh_3)
!  call give_two_p_pade_dh_4(mu_c,two_pade_dh_4)
!! call give_two_p_pade_dh_5(mu,two_pade_dh_5)
!! call give_two_p_pade_dh_6(mu,two_pade_dh_6)
!! call give_two_p_pade_dh_7(mu,two_pade_dh_7)
!! f_b=two_pade_dh_3(istate)
!  f_b_4=two_pade_dh_4(istate)*(1.d0-mu_c**2.d0)**(-3.d0/2.d0)
!! f_b_5=two_pade_dh_5(istate)
!! f_b_6=two_pade_dh_6(istate)
!! f_b_7=two_pade_dh_7(istate)

!
!!two_p_pade_h_3(istate) = -(h/3.d0)*(f_a+f_b+2.d0*part1 + 4.d0*part2 )
! two_p_pade_h_4(istate) = -(h/3.d0)*(f_a_4+f_b_4+2.d0*part1_4 + 4.d0*part2_4 )
!!two_p_pade_h_5(istate) = -(h/3.d0)*(f_a_5+f_b_5+2.d0*part1_5 + 4.d0*part2_5 )
!!two_p_pade_h_6(istate) = -(h/3.d0)*(f_a_6+f_b_6+2.d0*part1_6 + 4.d0*part2_6 )
!!two_p_pade_h_7(istate) = -(h/3.d0)*(f_a_7+f_b_7+2.d0*part1_7 + 4.d0*part2_7 )
!
!!print*,"two_p_pade_h_3 =",two_p_pade_h_3(istate)
!!print*,"f_b            =",f_b
!enddo
!END_PROVIDER



 BEGIN_PROVIDER[double precision, two_p_pade_h_4,(n_states)]
 implicit none
 integer :: istate,i,n
 double precision :: a,b,h,part1,part2,mu,mu_c,f_b,f_a,two_pade_dh_3(n_states)
 double precision :: part1_4,part2_4,f_b_4,f_a_4,two_pade_dh_4(n_states)
 double precision :: integral, error


 b=0.99999999d0
 a =sqrt(mu_erf**2.d0/(1.d0+mu_erf**2.d0))
!b=10000.d0
!a =mu_erf
 n= 10000000
 double precision, allocatable :: listx(:), listy(:)
 allocate(listx(n),listy(n))

 h = (b - a)/dble(n)

 do istate = 1, n_states

 !do i=1,n
 ! mu = a+dble(i)*h
 ! call give_two_p_pade_dh_4(mu,two_pade_dh_4)
 ! listx(i) = mu
 ! listy(i) = two_pade_dh_4(istate)
 !enddo

  do i=1,n
   mu = a+dble(i)*h
   mu_c =sqrt(mu**2.d0/(1-mu**2.d0))
   call give_two_p_pade_dh_4(mu_c,two_pade_dh_4)
   listx(i) = mu
   listy(i) = two_pade_dh_4(istate)*(1.d0-mu**2.d0)**(-3.d0/2.d0)
  enddo


  call DsplIntegr (listx, listy, n , a, b, integral, error)
 
  two_p_pade_h_4(istate) = -integral 


 enddo
 END_PROVIDER
