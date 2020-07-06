 BEGIN_PROVIDER [double precision, coef_xyz_ao, (2,3,ao_num)]
&BEGIN_PROVIDER [integer, power_xyz_ao, (2,3,ao_num) ]
 implicit none
 BEGIN_DOC
! coefficient for the basis function :: (x * phi_i(r), y * phi_i(r), * z_phi(r))
!
! x * (x - A_x)^a_x = A_x (x - A_x)^a_x + 1 * (x - A_x)^{a_x+1}
 END_DOC
 integer :: i,j,k,num_ao,power_ao(1:3)
 double precision :: center_ao(1:3)
 do i = 1, ao_num
  power_ao(1:3)= ao_power(i,1:3) 
  num_ao = ao_nucl(i)
  center_ao(1:3) = nucl_coord(num_ao,1:3)
  do j = 1, 3
   coef_xyz_ao(1,j,i) = center_ao(j) ! A_x (x - A_x)^a_x
   power_xyz_ao(1,j,i)= power_ao(j)
   coef_xyz_ao(2,j,i) = 1.d0         ! 1 * (x - A_x)^a_{x+1}
   power_xyz_ao(2,j,i)= power_ao(j) + 1
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, ao_coef_ord_grad_transp, (2,3,ao_prim_num_max,ao_num) ]
&BEGIN_PROVIDER [ integer, power_ord_grad_transp, (2,3,ao_num) ]
  implicit none
  BEGIN_DOC
  ! grad AO in terms of polynoms and coefficients 
  END_DOC
  integer                        :: i,j,power_ao(3), m
  do j=1, ao_num
    power_ao(1:3)= ao_power(j,1:3) 
    do i=1, ao_prim_num_max
     do m = 1, 3
      ao_coef_ord_grad_transp(1,m,i,j) = ao_coef_normalized_ordered(j,i) * dble(power_ao(m)) ! a_x * c_i 
      power_ord_grad_transp(1,m,j) = power_ao(m) - 1
      ao_coef_ord_grad_transp(2,m,i,j) = -2.d0 * ao_coef_normalized_ordered(j,i) * ao_expo_ordered_transp(i,j) ! -2 * c_i * alpha_i 
      power_ord_grad_transp(2,m,j) = power_ao(m) + 1
     enddo
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_coef_ord_xyz_grad_transp, (4,3,ao_prim_num_max,ao_num) ]
&BEGIN_PROVIDER [ integer, power_ord_xyz_grad_transp, (4,3,ao_num) ]
  implicit none
  BEGIN_DOC
  ! x * d/dx of an AO in terms of polynoms and coefficients 
  END_DOC
  integer                        :: i,j,power_ao(3), m,num_ao
  double precision :: center_ao(1:3)
  do j=1, ao_num
   power_ao(1:3)= ao_power(j,1:3) 
   num_ao = ao_nucl(j)
   center_ao(1:3) = nucl_coord(num_ao,1:3)
   do i=1, ao_prim_num_max
    do m = 1, 3
     power_ord_xyz_grad_transp(1,m,j)   = power_ao(m) - 1
     ao_coef_ord_xyz_grad_transp(1,m,i,j) = dble(power_ao(m)) * ao_coef_normalized_ordered(j,i) * center_ao(m)
     power_ord_xyz_grad_transp(2,m,j)   = power_ao(m)
     ao_coef_ord_xyz_grad_transp(2,m,i,j) = dble(power_ao(m)) * ao_coef_normalized_ordered(j,i) 
     power_ord_xyz_grad_transp(3,m,j)   = power_ao(m) + 1
     ao_coef_ord_xyz_grad_transp(3,m,i,j) = -2.d0 * ao_coef_normalized_ordered(j,i) * ao_expo_ordered_transp(i,j) * center_ao(m)
     power_ord_xyz_grad_transp(4,m,j)   = power_ao(m) + 2
     ao_coef_ord_xyz_grad_transp(4,m,i,j) = -2.d0 * ao_coef_normalized_ordered(j,i) * ao_expo_ordered_transp(i,j) 
    enddo
   enddo
  enddo

END_PROVIDER

subroutine xyz_grad_phi_ao(r,i_ao,xyz_grad_phi)
 implicit none
 integer, intent(in) :: i_ao
 double precision, intent(in) :: r(3)
 double precision, intent(out):: xyz_grad_phi(3) ! x * d/dx phi i, y * d/dy phi_i, z * d/dz phi_
 double precision :: center_ao(3),beta
 double precision :: accu(3,4),dr(3),r2,pol_usual(3)
 integer :: m,power_ao(3),num_ao,j_prim
 power_ao(1:3)= ao_power(i_ao,1:3) 
 num_ao = ao_nucl(i_ao)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dr(1) = (r(1) - center_ao(1))
 dr(2) = (r(2) - center_ao(2))
 dr(3) = (r(3) - center_ao(3))
 r2 = 0.d0
 do m = 1, 3
  r2 += dr(m)*dr(m)
 enddo
 ! computes the gaussian part 
 accu = 0.d0
 do j_prim =1,ao_prim_num(i_ao)
   beta = ao_expo_ordered_transp(j_prim,i_ao)
   if(dabs(beta*r2).gt.50.d0)cycle
   do m = 1, 3
    accu(m,1) += ao_coef_ord_xyz_grad_transp(1,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,2) += ao_coef_ord_xyz_grad_transp(2,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,3) += ao_coef_ord_xyz_grad_transp(3,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,4) += ao_coef_ord_xyz_grad_transp(4,m,j_prim,i_ao) * dexp(-beta*r2) 
   enddo
 enddo
 ! computes the polynom part
 pol_usual = 0.d0
 pol_usual(1) = dr(2)**dble(power_ao(2)) * dr(3)**dble(power_ao(3)) 
 pol_usual(2) = dr(1)**dble(power_ao(1)) * dr(3)**dble(power_ao(3)) 
 pol_usual(3) = dr(1)**dble(power_ao(1)) * dr(2)**dble(power_ao(2)) 

 xyz_grad_phi = 0.d0
 do m = 1, 3
  xyz_grad_phi(m) += accu(m,2) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(2,m,i_ao))
  xyz_grad_phi(m) += accu(m,3) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(3,m,i_ao))
  xyz_grad_phi(m) += accu(m,4) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(4,m,i_ao))
  if(power_ord_xyz_grad_transp(1,m,i_ao).lt.0)cycle
  xyz_grad_phi(m) += accu(m,1) * pol_usual(m) * dr(m)**dble(power_ord_xyz_grad_transp(1,m,i_ao))
 enddo
end

subroutine grad_phi_ao(r,i_ao,grad_xyz_phi)
 implicit none
 integer, intent(in) :: i_ao
 double precision, intent(in) :: r(3)
 double precision, intent(out):: grad_xyz_phi(3) ! x * phi i, y * phi_i, z * phi_
 double precision :: center_ao(3),beta
 double precision :: accu(3,2),dr(3),r2,pol_usual(3)
 integer :: m,power_ao(3),num_ao,j_prim
 power_ao(1:3)= ao_power(i_ao,1:3) 
 num_ao = ao_nucl(i_ao)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dr(1) = (r(1) - center_ao(1))
 dr(2) = (r(2) - center_ao(2))
 dr(3) = (r(3) - center_ao(3))
 r2 = 0.d0
 do m = 1, 3
  r2 += dr(m)*dr(m)
 enddo
 ! computes the gaussian part 
 accu = 0.d0
 do j_prim =1,ao_prim_num(i_ao)
   beta = ao_expo_ordered_transp(j_prim,i_ao)
   if(dabs(beta*r2).gt.50.d0)cycle
   do m = 1, 3
    accu(m,1) += ao_coef_ord_grad_transp(1,m,j_prim,i_ao) * dexp(-beta*r2) 
    accu(m,2) += ao_coef_ord_grad_transp(2,m,j_prim,i_ao) * dexp(-beta*r2) 
   enddo
 enddo
 ! computes the polynom part
 pol_usual = 0.d0
 pol_usual(1) = dr(2)**dble(power_ao(2)) * dr(3)**dble(power_ao(3)) 
 pol_usual(2) = dr(1)**dble(power_ao(1)) * dr(3)**dble(power_ao(3)) 
 pol_usual(3) = dr(1)**dble(power_ao(1)) * dr(2)**dble(power_ao(2)) 
 do m = 1, 3
  grad_xyz_phi(m)  = accu(m,2) * pol_usual(m) * dr(m)**dble(power_ord_grad_transp(2,m,i_ao))
  if(power_ao(m)==0)cycle
  grad_xyz_phi(m) += accu(m,1) * pol_usual(m) * dr(m)**dble(power_ord_grad_transp(1,m,i_ao))
 enddo
end

subroutine xyz_phi_ao(r,i_ao,xyz_phi)
 implicit none
 integer, intent(in) :: i_ao
 double precision, intent(in) :: r(3)
 double precision, intent(out):: xyz_phi(3) ! x * phi i, y * phi_i, z * phi_i
 double precision :: center_ao(3),beta
 double precision :: accu,dr(3),r2,pol_usual(3)
 integer :: m,power_ao(3),num_ao
 power_ao(1:3)= ao_power(i_ao,1:3) 
 num_ao = ao_nucl(i_ao)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dr(1) = (r(1) - center_ao(1))
 dr(2) = (r(2) - center_ao(2))
 dr(3) = (r(3) - center_ao(3))
 r2 = 0.d0
 do m = 1, 3
  r2 += dr(m)*dr(m)
 enddo
 ! computes the gaussian part 
 accu = 0.d0
 do m=1,ao_prim_num(i_ao)
   beta = ao_expo_ordered_transp(m,i_ao)
   if(dabs(beta*r2).gt.50.d0)cycle
   accu += ao_coef_normalized_ordered_transp(m,i_ao) * dexp(-beta*r2)
 enddo
 ! computes the polynom part
 pol_usual = 0.d0
 pol_usual(1) = dr(2)**dble(power_ao(2)) * dr(3)**dble(power_ao(3)) 
 pol_usual(2) = dr(1)**dble(power_ao(1)) * dr(3)**dble(power_ao(3)) 
 pol_usual(3) = dr(1)**dble(power_ao(1)) * dr(2)**dble(power_ao(2)) 
 do m = 1, 3
  xyz_phi(m) = accu * pol_usual(m) * dr(m)**(dble(power_ao(m))) * ( coef_xyz_ao(1,m,i_ao) + coef_xyz_ao(2,m,i_ao) * dr(m) )
 enddo
end
