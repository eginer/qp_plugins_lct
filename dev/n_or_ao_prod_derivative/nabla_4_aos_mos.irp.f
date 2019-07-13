
subroutine give_all_aos_and_fourth_at_r(r,aos_array,aos_grad_array,aos_lapl_array,aos_3rd_direct_array,aos_4th_direct_array)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output     : aos_array(i) = ao(i) evaluated at r
!
!            : aos_grad_array(1,i) = gradient X of the ao(i) evaluated at r
!
!            : aos_lapl_array(1,i) =  X componetnt laplcian of the ao(i) evaluated at r (d^2ao(i)/dx^2)
!
!            : aos_3rd_direct_array(1,i) = Third partial derivative of the ao(i) evaluated at r by X (d^3ao(i)/dx^3)
!
!            : aos_3rd_direct_array(1,i) = Fourth partial derivative of the ao(i) evaluated at r by X (d^4ao(i)/dx^4)
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: aos_array(ao_num)
 double precision, intent(out) :: aos_grad_array(3,ao_num)
 double precision, intent(out) :: aos_lapl_array(3,ao_num)
 double precision, intent(out) :: aos_3rd_direct_array(3,ao_num)
 double precision, intent(out) :: aos_4th_direct_array(3,ao_num)

 integer :: power_ao(3)
 integer :: i,j,k,l,m
 double precision :: dx,dy,dz,r2
 double precision ::      dx2,dy2,dz2
 double precision ::      dx1,dy1,dz1
 double precision ::      dx3,dy3,dz3
 double precision ::      dx4,dy4,dz4
 double precision ::      dx5,dy5,dz5
 double precision ::      tx1,tx2,tx3,tx4
 double precision ::      ty1,ty2,ty3,ty4
 double precision ::      tz1,tz2,tz3,tz4
 double precision ::      fx1,fx2,fx3,fx4,fx5
 double precision ::      fy1,fy2,fy3,fy4,fy5
 double precision ::      fz1,fz2,fz3,fz4,fz5
 double precision :: center_ao(3)
 double precision :: beta,accu_1,accu_2,accu_3,accu_4,accu_5,contrib
 do i = 1, nucl_num
  center_ao(1:3) = nucl_coord(i,1:3)
  dx = (r(1) - center_ao(1))
  dy = (r(2) - center_ao(2))
  dz = (r(3) - center_ao(3))
  r2 = dx*dx + dy*dy + dz*dz
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i) ! index of the ao in the ordered format
   aos_array(k) = 0.d0
   aos_grad_array(1,k) = 0.d0
   aos_grad_array(2,k) = 0.d0
   aos_grad_array(3,k) = 0.d0

   aos_lapl_array(1,k) = 0.d0
   aos_lapl_array(2,k) = 0.d0
   aos_lapl_array(3,k) = 0.d0

   aos_3rd_direct_array(1,k) = 0.d0
   aos_3rd_direct_array(2,k) = 0.d0
   aos_3rd_direct_array(3,k) = 0.d0

   aos_4th_direct_array(1,k) = 0.d0
   aos_4th_direct_array(2,k) = 0.d0
   aos_4th_direct_array(3,k) = 0.d0
  
   power_ao(1:3)= ao_power_ordered_transp_per_nucl(1:3,j,i)
   dx2 = dx**power_ao(1)
   dy2 = dy**power_ao(2)
   dz2 = dz**power_ao(3)
   if(power_ao(1) .ne. 0)then
    dx1 = dble(power_ao(1)) * dx**(power_ao(1)-1)
   else
    dx1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(1) .ge. 2)then
    dx3 = dble(power_ao(1)) * dble((power_ao(1)-1))  * dx**(power_ao(1)-2)
   else
    dx3 = 0.d0
   endif
   dx4 = dble((2 * power_ao(1) + 1))  * dx**(power_ao(1)) 

   dx5 = dx**(power_ao(1)+2)

 ! For the third order 
   if(power_ao(1) .ge. 3)then
    tx1 = dble(power_ao(1)) * dble((power_ao(1)-1))*dble((power_ao(1)-2))  * dx**(power_ao(1)-3)
   else
    tx1 = 0.d0
   endif

   if(power_ao(1) .ge. 1)then
    tx2 = dble(power_ao(1)) * 3.d0 * dble(power_ao(1)) * dx**(power_ao(1)-1)
   else
    tx2 = 0.d0
   endif
   tx3 =  (3.d0*dble(power_ao(1)) + 3.d0 ) * dx**(power_ao(1)+1) 
   tx4 = dx**(power_ao(1)+3)


 ! For the fourth order 
   if(power_ao(1) .ge. 4)then
    fx1 = dble(power_ao(1)) * dble((power_ao(1)-1))*dble((power_ao(1)-2))*dble((power_ao(1)-3)) * dx**(power_ao(1)-4)
   else
    fx1 = 0.d0
   endif

   if(power_ao(1) .ge. 2)then
    fx2 = dble(power_ao(1)) * (dble(power_ao(1))- 1d0) * (4d0 * dble(power_ao(1))- 2d0) * dx**(power_ao(1)-2)
   else
    fx2 = 0.d0
   endif
   fx3 =  (12.d0*dble(power_ao(1))**2d0 + 4.d0*(3d0*dble(power_ao(1))+3d0)*(dble(power_ao(1))+1d0)) * dx**(power_ao(1)) 
   fx4 = (4d0*dble(power_ao(1))+6d0 ) * dx**(power_ao(1)+2)
   fx5 =  dx**(power_ao(1)+4)   

   !!!!!!!!!!!!!!! For Y!!!!!!!!!!!!!!!!!!!

   if(power_ao(2) .ne. 0)then
    dy1 = dble(power_ao(2)) * dy**(power_ao(2)-1)
   else
    dy1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(2) .ge. 2)then
    dy3 = dble(power_ao(2)) * dble((power_ao(2)-1))  * dy**(power_ao(2)-2)
   else
    dy3 = 0.d0
   endif

   dy4 = dble((2 * power_ao(2) + 1))  * dy**(power_ao(2)) 
   dy5 = dy**(power_ao(2)+2)


 ! For the third order 
   if(power_ao(2) .ge. 3)then
    ty1 = dble(power_ao(2)) * dble((power_ao(2)-1))*dble((power_ao(2)-2))  * dy**(power_ao(2)-3)
   else
    ty1 = 0.d0
   endif

   if(power_ao(2) .ge. 1)then
    ty2 = dble(power_ao(2)) * 3.d0 * dble(power_ao(2)) * dy**(power_ao(2)-1)
   else
    ty2 = 0.d0
   endif
   ty3 =  (3.d0*dble(power_ao(2)) + 3.d0 ) * dy**(power_ao(2)+1) 
   ty4 = dy**(power_ao(2)+3)
 
 ! For the fourth order 
   if(power_ao(2) .ge. 4)then
    fy1 = dble(power_ao(2)) * dble((power_ao(2)-1))*dble((power_ao(2)-2))*dble((power_ao(2)-3)) * dy**(power_ao(2)-4)
   else
    fy1 = 0.d0
   endif

   if(power_ao(2) .ge. 2)then
    fy2 = dble(power_ao(2)) * (dble(power_ao(2))- 1d0) * (4d0 * dble(power_ao(2))- 2d0) * dy**(power_ao(2)-2)
   else
    fy2 = 0.d0
   endif
   fy3 =  (12.d0*dble(power_ao(2))**2d0 + 4.d0*(3d0*dble(power_ao(2))+3d0)*(dble(power_ao(2))+1d0)) * dy**(power_ao(2)) 
   fy4 = (4d0*dble(power_ao(2))+6d0 ) * dy**(power_ao(2)+2)
   fy5 =  dy**(power_ao(2)+4)   

   !!!!!!!!!!!!!!! For Z!!!!!!!!!!!!!!!!!!!


  if(power_ao(3) .ne. 0)then
    dz1 = dble(power_ao(3)) * dz**(power_ao(3)-1)
   else
    dz1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(3) .ge. 2)then
    dz3 = dble(power_ao(3)) * dble((power_ao(3)-1))  * dz**(power_ao(3)-2)
   else
    dz3 = 0.d0
   endif

   dz4 = dble((2 * power_ao(3) + 1))  * dz**(power_ao(3)) 
   dz5 = dz**(power_ao(3)+2)

 ! For the third order 
   if(power_ao(3) .ge. 3)then
    tz1 = dble(power_ao(3)) * dble((power_ao(3)-1))*dble((power_ao(3)-2))  * dz**(power_ao(3)-3)
   else
    tz1 = 0.d0
   endif

   if(power_ao(3) .ge. 1)then
    tz2 = dble(power_ao(3)) * 3.d0 * dble(power_ao(3)) * dz**(power_ao(3)-1)
   else
    tz2 = 0.d0
   endif
   tz3 =  (3.d0*dble(power_ao(3)) + 3.d0 ) * dz**(power_ao(3)+1) 
   tz4 = dz**(power_ao(3)+3)

 ! For the fourth order 
   if(power_ao(3) .ge. 4)then
    fz1 = dble(power_ao(3)) * dble((power_ao(3)-1))*dble((power_ao(3)-2))*dble((power_ao(3)-3)) * dz**(power_ao(3)-4)
   else
    fz1 = 0.d0
   endif

   if(power_ao(3) .ge. 2)then
    fz2 = dble(power_ao(3)) * (dble(power_ao(3))- 1d0) * (4d0 * dble(power_ao(3))- 2d0) * dz**(power_ao(3)-2)
   else
    fz2 = 0.d0
   endif
   fz3 =  (12.d0*dble(power_ao(3))**2d0 + 4.d0*(3d0*dble(power_ao(3))+3d0)*(dble(power_ao(3))+1d0)) * dz**(power_ao(3)) 
   fz4 = (4d0*dble(power_ao(3))+6d0 ) * dz**(power_ao(3)+2)
   fz5 =  dz**(power_ao(3)+4)   


   accu_1 = 0.d0
   accu_2 = 0.d0
   accu_3 = 0.d0
   accu_4 = 0.d0
   accu_5 = 0.d0
   do l = 1,ao_prim_num(k)
    beta = ao_expo_ordered_transp_per_nucl(l,j,i)
    contrib = ao_coef_normalized_ordered_transp_per_nucl(l,j,i) * dexp(-beta*r2)
    accu_1 += contrib
    accu_2 += contrib * beta
    accu_3 += contrib * beta**2
    accu_4 += contrib * beta**3
    accu_5 += contrib * beta**4
   enddo


   aos_array(k) = accu_1 * dx2 * dy2 * dz2

   aos_grad_array(1,k) = accu_1 * dx1  * dy2 * dz2- 2.d0 * dx2 * dx  * dy2 * dz2 * accu_2
   aos_grad_array(2,k) = accu_1 * dx2  * dy1 * dz2- 2.d0 * dx2 * dy2 * dy  * dz2 * accu_2
   aos_grad_array(3,k) = accu_1 * dx2  * dy2 * dz1- 2.d0 * dx2 * dy2 * dz2 * dz  * accu_2

   aos_lapl_array(1,k) = accu_1 * dx3  * dy2 * dz2- 2.d0 * dx4 * dy2 * dz2* accu_2 +4.d0 * dx5 *dy2 * dz2* accu_3
   aos_lapl_array(2,k) = accu_1 * dx2  * dy3 * dz2- 2.d0 * dx2 * dy4 * dz2* accu_2 +4.d0 * dx2 *dy5 * dz2* accu_3
   aos_lapl_array(3,k) = accu_1 * dx2  * dy2 * dz3- 2.d0 * dx2 * dy2 * dz4* accu_2 +4.d0 * dx2 *dy2 * dz5* accu_3

   aos_3rd_direct_array(1,k) = dy2 * dz2 *( tx1 * accu_1 - 2.d0 * tx2 * accu_2 + 4.d0 * tx3 * accu_3 -  8.d0 * tx4 * accu_4 ) 
   aos_3rd_direct_array(2,k) = dx2 * dz2 *( ty1 * accu_1 - 2.d0 * ty2 * accu_2 + 4.d0 * ty3 * accu_3 -  8.d0 * ty4 * accu_4 ) 
   aos_3rd_direct_array(3,k) = dx2 * dy2 *( tz1 * accu_1 - 2.d0 * tz2 * accu_2 + 4.d0 * tz3 * accu_3 -  8.d0 * tz4 * accu_4 ) 

   aos_4th_direct_array(1,k) = dy2 * dz2 *( fx1 * accu_1 - 2.d0 * fx2 * accu_2 + fx3 * accu_3 - 8.d0 * fx4 * accu_4 +  16.d0 * fx5 * accu_5)
   aos_4th_direct_array(2,k) = dx2 * dz2 *( fy1 * accu_1 - 2.d0 * fy2 * accu_2 + fy3 * accu_3 - 8.d0 * fy4 * accu_4 +  16.d0 * fy5 * accu_5) 
   aos_4th_direct_array(3,k) = dx2 * dy2 *( fz1 * accu_1 - 2.d0 * fz2 * accu_2 + fz3 * accu_3 - 8.d0 * fz4 * accu_4 +  16.d0 * fz5 * accu_5)  
  enddo
 enddo
end



 subroutine give_all_aos_and_fourth_order_cross_terms_at_r(r,aos_array,aos_grad_array,aos_2nd_cross_array,aos_3rd_cross_array,aos_4th_cross_array)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output     : aos_array(i) = ao(i) evaluated at r
!
!            : aos_grad_array(1,i) = gradient X of the ao(i) evaluated at r
!
!            : aos_2nd_cross_array(1,i) =  Second order partial derivative of the ao(i) evaluated at r by X and Y (d^2ao(i)/dydx) 
!
!            : aos_2nd_cross_array(2,i) =  Second order partial derivative of the ao(i) evaluated at r by Y and Z (d^2ao(i)/dydz) 
!
!            : aos_2nd_cross_array(3,i) =  Second order partial derivative of the ao(i) evaluated at r by X and Z (d^2ao(i)/dxdz)
!
!            : aos_3rd_cross_array(1,i) = Third partial derivative of the ao(i) evaluated at r by X and Y (d^3ao(i)/dydx^2)
!
!            : aos_3rd_cross_array(2,i) = Third partial derivative of the ao(i) evaluated at r by Y and Z (d^3ao(i)/dzdy^2)
!
!            : aos_3rd_cross_array(3,i) = Third partial derivative of the ao(i) evaluated at r by Z and X (d^3ao(i)/dxdz^2)

!            : aos_3rd_cross_array(4,i) = Third partial derivative of the ao(i) evaluated at r by X and Y (d^3ao(i)/dxdy^2)
!
!            : aos_3rd_cross_array(5,i) = Third partial derivative of the ao(i) evaluated at r by Y and Z (d^3ao(i)/dydz^2)
!
!            : aos_3rd_cross_array(6,i) = Third partial derivative of the ao(i) evaluated at r by Z and X (d^3ao(i)/dzdx^2)
!
!            : aos_4th_cross_array(1,i) = Fourth partial derivative of the ao(i) evaluated at r by X and Y (d^4ao(i)/dx^2dy^2)
!
!            : aos_4th_cross_array(2,i) = Fourth partial derivative of the ao(i) evaluated at r by Y and Z (d^4ao(i)/dy^2dz^2)
!
!            : aos_4th_cross_array(3,i) = Fourth partial derivative of the ao(i) evaluated at r by Z and X (d^4ao(i)/dx^2dz^2)
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: aos_array(ao_num)
 double precision, intent(out) :: aos_grad_array(3,ao_num)
 double precision, intent(out) :: aos_2nd_cross_array(3,ao_num)
 double precision, intent(out) :: aos_3rd_cross_array(6,ao_num)
 double precision, intent(out) :: aos_4th_cross_array(3,ao_num)

 integer :: power_ao(3)
 integer :: i,j,k,l,m
 double precision :: dx,dy,dz,r2
 double precision ::      dx2,dy2,dz2
 double precision ::      dx1,dy1,dz1
 double precision ::      dx3,dy3,dz3
 double precision ::      dx4,dy4,dz4
 double precision ::      dx5,dy5,dz5
 double precision ::      tx1,tx2,tx3,tx4
 double precision ::      ty1,ty2,ty3,ty4
 double precision ::      tz1,tz2,tz3,tz4
 double precision ::      fx1,fx2,fx3,fx4,fx5
 double precision :: center_ao(3)
 double precision :: beta,accu_1,accu_2,accu_3,accu_4,accu_5,contrib

 aos_array = 0.d0
 aos_grad_array = 0.d0    
 aos_2nd_cross_array = 0.d0
 aos_3rd_cross_array = 0.d0
 aos_4th_cross_array = 0.d0

 do i = 1, nucl_num
  center_ao(1:3) = nucl_coord(i,1:3)
  dx = (r(1) - center_ao(1))
  dy = (r(2) - center_ao(2))
  dz = (r(3) - center_ao(3))
  r2 = dx*dx + dy*dy + dz*dz
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i) ! index of the ao in the ordered format

   power_ao(1:3)= ao_power_ordered_transp_per_nucl(1:3,j,i)
   dx2 = dx**power_ao(1)
   dy2 = dy**power_ao(2)
   dz2 = dz**power_ao(3)
   if(power_ao(1) .ne. 0)then
    dx1 = dble(power_ao(1)) * dx**(power_ao(1)-1)
   else
    dx1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(1) .ge. 2)then
    dx3 = dble(power_ao(1)) * dble((power_ao(1)-1))  * dx**(power_ao(1)-2)
   else
    dx3 = 0.d0
   endif
   dx4 = dble((2 * power_ao(1) + 1))  * dx**(power_ao(1)) 

   dx5 = dx**(power_ao(1)+2)

 ! For the third order 
   if(power_ao(1) .ge. 3)then
    tx1 = dble(power_ao(1)) * dble((power_ao(1)-1))*dble((power_ao(1)-2))  * dx**(power_ao(1)-3)
   else
    tx1 = 0.d0
   endif

   if(power_ao(1) .ge. 1)then
    tx2 = dble(power_ao(1)) * 3.d0 * dble(power_ao(1)) * dx**(power_ao(1)-1)
   else
    tx2 = 0.d0
   endif
   tx3 =  (3.d0*dble(power_ao(1)) + 3.d0 ) * dx**(power_ao(1)+1) 
   tx4 = dx**(power_ao(1)+3)

   !!!!!!!!!!!!!!! For Y!!!!!!!!!!!!!!!!!!!

   if(power_ao(2) .ne. 0)then
    dy1 = dble(power_ao(2)) * dy**(power_ao(2)-1)
   else
    dy1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(2) .ge. 2)then
    dy3 = dble(power_ao(2)) * dble((power_ao(2)-1))  * dy**(power_ao(2)-2)
   else
    dy3 = 0.d0
   endif

   dy4 = dble((2 * power_ao(2) + 1))  * dy**(power_ao(2)) 
   dy5 = dy**(power_ao(2)+2)

   !!!!!!!!!!!!!!! For Z!!!!!!!!!!!!!!!!!!!

  if(power_ao(3) .ne. 0)then
    dz1 = dble(power_ao(3)) * dz**(power_ao(3)-1)
   else
    dz1 = 0.d0
   endif
   ! For the Laplacian
   if(power_ao(3) .ge. 2)then
    dz3 = dble(power_ao(3)) * dble((power_ao(3)-1))  * dz**(power_ao(3)-2)
   else
    dz3 = 0.d0
   endif

   dz4 = dble((2 * power_ao(3) + 1))  * dz**(power_ao(3)) 
   dz5 = dz**(power_ao(3)+2)

   accu_1 = 0.d0
   accu_2 = 0.d0
   accu_3 = 0.d0
   accu_4 = 0.d0
   accu_5 = 0.d0
   do l = 1,ao_prim_num(k)
    beta = ao_expo_ordered_transp_per_nucl(l,j,i)
    contrib = ao_coef_normalized_ordered_transp_per_nucl(l,j,i) * dexp(-beta*r2)
    accu_1 += contrib
    accu_2 += contrib * beta
    accu_3 += contrib * beta**2
    accu_4 += contrib * beta**3
    accu_5 += contrib * beta**4
   enddo


   aos_array(k) = accu_1 * dx2 * dy2 * dz2

   aos_grad_array(1,k) = accu_1 * dx1  * dy2 * dz2- 2.d0 * dx2 * dx  * dy2 * dz2 * accu_2
   aos_grad_array(2,k) = accu_1 * dx2  * dy1 * dz2- 2.d0 * dx2 * dy2 * dy  * dz2 * accu_2
   aos_grad_array(3,k) = accu_1 * dx2  * dy2 * dz1- 2.d0 * dx2 * dy2 * dz2 * dz  * accu_2

   aos_2nd_cross_array(1,k) = dz2 *( dx1*(dy1*accu_1 - 2d0* dy2 * dy*accu_2)-2.d0 * dx2 * dx* (dy1*accu_2 - 2d0* dy2 * dy*accu_3 ) )
   aos_2nd_cross_array(2,k) = dx2 *( dy1*(dz1*accu_1 - 2d0* dz2 * dz*accu_2)-2.d0 * dy2 * dy* (dz1*accu_2 - 2d0* dz2 * dz*accu_3 ) ) 
   aos_2nd_cross_array(3,k) = dy2 *( dz1*(dx1*accu_1 - 2d0* dx2 * dx*accu_2)-2.d0 * dz2 * dz* (dx1*accu_2 - 2d0* dx2 * dx*accu_3 ) ) 

   aos_3rd_cross_array(1,k) = dz2 *( dx3 * (  dy1 * accu_1 - 2d0 *  dy2 * dy * accu_2 ) - 2d0 * dx4 *( dy1 * accu_2 - 2d0 *  dy2 * dy * accu_3 ) + 4d0* dx5 *( dy1 * accu_3 - 2d0 *  dy2 * dy * accu_4 )    ) 
   aos_3rd_cross_array(2,k) = dx2 *( dy3 * (  dz1 * accu_1 - 2d0 *  dz2 * dz * accu_2 ) - 2d0 * dy4 *( dz1 * accu_2 - 2d0 *  dz2 * dz * accu_3 ) + 4d0* dy5 *( dz1 * accu_3 - 2d0 *  dz2 * dz * accu_4 )    ) 
   aos_3rd_cross_array(3,k) = dy2 *( dz3 * (  dx1 * accu_1 - 2d0 *  dx2 * dx *accu_2 ) - 2d0 * dz4 *( dx1 * accu_2 - 2d0 *  dx2 * dx * accu_3 ) + 4d0* dz5 *( dx1 * accu_3 - 2d0 *  dx2 * dx * accu_4 )) 
   aos_3rd_cross_array(4,k) = dz2 *( dy3 * (  dx1 * accu_1 - 2d0 *  dx2 * dx * accu_2 ) - 2d0 * dy4 *( dx1 * accu_2 - 2d0 *  dx2 * dx * accu_3 ) + 4d0* dy5 *( dx1 * accu_3 - 2d0 *  dx2 * dx * accu_4 )    ) 
   aos_3rd_cross_array(5,k) = dx2 *( dz3 * (  dy1 * accu_1 - 2d0 *  dy2 * dy * accu_2 ) - 2d0 * dz4 *( dy1 * accu_2 - 2d0 *  dy2 * dy * accu_3 ) + 4d0* dz5 *( dy1 * accu_3 - 2d0 *  dy2 * dy * accu_4 )    ) 
   aos_3rd_cross_array(6,k) = dy2 *( dx3 * (  dz1 * accu_1 - 2d0 *  dz2 * dz *accu_2 ) - 2d0 * dx4 *( dz1 * accu_2 - 2d0 *  dz2 * dz * accu_3 ) + 4d0* dx5 *( dz1 * accu_3 - 2d0 *  dz2 * dz * accu_4 )) 



   aos_4th_cross_array(1,k) = dz2 *( dx3 * (dy3*accu_1- 2d0 * dy4*accu_2 + 4d0 *dy5*accu_3) - 2d0 *dx4*(dy3*accu_2 - 2d0 * dy4*accu_3 + 4d0 *dy5*accu_4) + 4d0 *dx5*(dy3*accu_3 - 2d0 * dy4*accu_4 + 4d0 *dy5*accu_5) ) 
   aos_4th_cross_array(2,k) = dx2 *( dy3 * (dz3*accu_1- 2d0 * dz4*accu_2 + 4d0 *dz5*accu_3) - 2d0 *dy4*(dz3*accu_2 - 2d0 * dz4*accu_3 + 4d0 *dz5*accu_4) + 4d0 *dy5*(dz3*accu_3 - 2d0 * dz4*accu_4 + 4d0 *dz5*accu_5) )
   aos_4th_cross_array(3,k) = dy2 *( dz3 * (dx3*accu_1- 2d0 * dx4*accu_2 + 4d0 *dx5*accu_3) - 2d0 *dz4*(dx3*accu_2 - 2d0 * dx4*accu_3 + 4d0 *dx5*accu_4) + 4d0 *dz5*(dx3*accu_3 - 2d0 * dx4*accu_4 + 4d0 *dx5*accu_5) ) 
  enddo
 enddo
end



 subroutine give_nabla_2_at_r(r,nabla_2_at_r,nabla_2_tot_at_r)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output     : aos_array(i) = ao(i) evaluated at r
!
!            : blablablablablabla 
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: nabla_2_at_r(3,ao_num,ao_num) 
 double precision, intent(out) :: nabla_2_tot_at_r(ao_num,ao_num)

 integer :: i,j,m

 double precision :: aos_array(ao_num)
 double precision :: aos_grad_array(3,ao_num)
 double precision :: aos_lapl_array(3,ao_num)

 call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)

 do j = 1,ao_num
  do i = 1,ao_num
   do m = 1,3 

    nabla_2_at_r(m,i,j) = aos_lapl_array(m,i)*aos_array(j)+2d0*aos_grad_array(m,i)*aos_grad_array(m,j)+aos_array(i)*aos_lapl_array(m,j)

   enddo
   nabla_2_tot_at_r(i,j) = nabla_2_at_r(1,i,j) + nabla_2_at_r(2,i,j) + nabla_2_at_r(3,i,j)
  enddo
 enddo

 end



 subroutine give_nabla_4_at_r(r,nabla_4_at_r)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output     : aos_array(i) = ao(i) evaluated at r
!
!            : blablablablablabla 
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: nabla_4_at_r(ao_num,ao_num) 

 integer :: i,j
 double precision :: d2xd2y,d2yd2z,d2xd2z
 double precision :: d4_x,d4_y,d4_z

 double precision :: aos_array(ao_num)
 double precision :: aos_grad_array(3,ao_num)

 double precision :: aos_lapl_array(3,ao_num)
 double precision :: aos_3rd_direct_array(3,ao_num)
 double precision :: aos_4th_direct_array(3,ao_num)

 double precision :: aos_2nd_cross_array(3,ao_num)
 double precision :: aos_3rd_cross_array(6,ao_num)
 double precision :: aos_4th_cross_array(3,ao_num)

 call give_all_aos_and_fourth_at_r(r,aos_array,aos_grad_array,aos_lapl_array,aos_3rd_direct_array,aos_4th_direct_array)
 call give_all_aos_and_fourth_order_cross_terms_at_r(r,aos_array,aos_grad_array,aos_2nd_cross_array,aos_3rd_cross_array,aos_4th_cross_array)
 call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)
 
 do i = 1,ao_num
  do j = 1,ao_num
  
   d2xd2y = aos_4th_cross_array(1,i)*aos_array(j) + aos_lapl_array(1,i)*aos_lapl_array(2,j)+ aos_lapl_array(2,i)*aos_lapl_array(1,j)+aos_array(i)*aos_4th_cross_array(1,j) + 2d0*(aos_3rd_cross_array(1,i)*aos_grad_array(2,j) + aos_3rd_cross_array(4,i)*aos_grad_array(1,j)+ aos_grad_array(1,i)*aos_3rd_cross_array(4,j) + aos_grad_array(2,i)*aos_3rd_cross_array(1,j)+ 2d0* aos_2nd_cross_array(1,i)*aos_2nd_cross_array(1,j)) 

   d2yd2z = aos_4th_cross_array(2,i)*aos_array(j) + aos_lapl_array(2,i)*aos_lapl_array(3,j)+ aos_lapl_array(3,i)*aos_lapl_array(2,j)+aos_array(i)*aos_4th_cross_array(2,j)+ 2d0*(aos_3rd_cross_array(2,i)*aos_grad_array(3,j) + aos_3rd_cross_array(5,i)*aos_grad_array(2,j)+ aos_grad_array(2,i)*aos_3rd_cross_array(5,j)+ aos_grad_array(3,i)*aos_3rd_cross_array(2,j)+ 2d0*aos_2nd_cross_array(2,i)*aos_2nd_cross_array(2,j))

   d2xd2z = aos_4th_cross_array(3,i)*aos_array(j) + aos_lapl_array(1,i)*aos_lapl_array(3,j)+ aos_lapl_array(3,i)*aos_lapl_array(1,j)+aos_array(i)*aos_4th_cross_array(3,j)+ 2d0*(aos_3rd_cross_array(6,i)*aos_grad_array(3,j) + aos_3rd_cross_array(3,i)*aos_grad_array(1,j)+ aos_grad_array(1,i)*aos_3rd_cross_array(3,j)+ aos_grad_array(3,i)*aos_3rd_cross_array(6,j)+ 2d0* aos_2nd_cross_array(3,i)*aos_2nd_cross_array(3,j))
 

   d4_x = aos_4th_direct_array(1,i)*aos_array(j)+aos_array(i)*aos_4th_direct_array(1,j)+ 2d0*(2d0*aos_3rd_direct_array(1,i)*aos_grad_array(1,j)+2d0*aos_grad_array(1,i)*aos_3rd_direct_array(1,j)+ 3d0*aos_lapl_array(1,i)*aos_lapl_array(1,j) )
   
   d4_y = aos_4th_direct_array(2,i)*aos_array(j)+aos_array(i)*aos_4th_direct_array(2,j)+ 2d0*(2d0*aos_3rd_direct_array(2,i)*aos_grad_array(2,j)+2d0*aos_grad_array(2,i)*aos_3rd_direct_array(2,j)+ 3d0*aos_lapl_array(2,i)*aos_lapl_array(2,j) ) 

   d4_z = aos_4th_direct_array(3,i)*aos_array(j)+aos_array(i)*aos_4th_direct_array(3,j)+ 2d0*(2d0*aos_3rd_direct_array(3,i)*aos_grad_array(3,j)+2d0*aos_grad_array(3,i)*aos_3rd_direct_array(3,j)+ 3d0*aos_lapl_array(3,i)*aos_lapl_array(3,j) )

   nabla_4_at_r(i,j) = d4_x + d4_y + d4_z + 2.d0 *(d2xd2y+d2yd2z+d2xd2z) 

  enddo
 enddo

 end


 subroutine give_nabla_3_at_r(r,nabla_3_at_r)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output     : aos_array(i) = ao(i) evaluated at r
!
!            : blablablablablabla 
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: nabla_3_at_r(9,ao_num,ao_num) 

 integer :: i,j

 double precision :: aos_array(ao_num)
 double precision :: aos_grad_array(3,ao_num)

 double precision :: aos_lapl_array(3,ao_num)
 double precision :: aos_3rd_direct_array(3,ao_num)
 double precision :: aos_4th_direct_array(3,ao_num)

 double precision :: aos_2nd_cross_array(3,ao_num)
 double precision :: aos_3rd_cross_array(6,ao_num)
 double precision :: aos_4th_cross_array(3,ao_num)

 call give_all_aos_and_fourth_at_r(r,aos_array,aos_grad_array,aos_lapl_array,aos_3rd_direct_array,aos_4th_direct_array)
 call give_all_aos_and_fourth_order_cross_terms_at_r(r,aos_array,aos_grad_array,aos_2nd_cross_array,aos_3rd_cross_array,aos_4th_cross_array)
 call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)
 
 do i = 1,ao_num
  do j = 1,ao_num
  
   nabla_3_at_r(1,i,j) = aos_3rd_cross_array(1,i)*aos_array(j)+aos_lapl_array(1,i)*aos_grad_array(2,j)+2.d0*aos_2nd_cross_array(1,i)*aos_grad_array(1,j)+2.d0*aos_grad_array(1,i)*aos_2nd_cross_array(1,j)+aos_grad_array(2,i)*aos_lapl_array(1,j)+aos_array(i)*aos_3rd_cross_array(1,j)

   nabla_3_at_r(2,i,j) = aos_3rd_cross_array(2,i)*aos_array(j)+aos_lapl_array(2,i)*aos_grad_array(3,j)+2.d0*aos_2nd_cross_array(2,i)*aos_grad_array(2,j)+2.d0*aos_grad_array(2,i)*aos_2nd_cross_array(2,j)+aos_grad_array(3,i)*aos_lapl_array(2,j)+aos_array(i)*aos_3rd_cross_array(2,j)

   nabla_3_at_r(3,i,j) = aos_3rd_cross_array(3,i)*aos_array(j)+aos_lapl_array(3,i)*aos_grad_array(1,j)+2.d0*aos_2nd_cross_array(3,i)*aos_grad_array(3,j)+2.d0*aos_grad_array(3,i)*aos_2nd_cross_array(3,j)+aos_grad_array(1,i)*aos_lapl_array(3,j)+aos_array(i)*aos_3rd_cross_array(3,j)

   nabla_3_at_r(4,i,j) = aos_3rd_cross_array(4,i)*aos_array(j)+aos_lapl_array(2,i)*aos_grad_array(1,j)+2.d0*aos_2nd_cross_array(1,i)*aos_grad_array(2,j)+2.d0*aos_grad_array(2,i)*aos_2nd_cross_array(1,j)+aos_grad_array(1,i)*aos_lapl_array(2,j)+aos_array(i)*aos_3rd_cross_array(4,j)

   nabla_3_at_r(5,i,j) = aos_3rd_cross_array(5,i)*aos_array(j)+aos_lapl_array(3,i)*aos_grad_array(2,j)+2.d0*aos_2nd_cross_array(2,i)*aos_grad_array(3,j)+2.d0*aos_grad_array(3,i)*aos_2nd_cross_array(2,j)+aos_grad_array(2,i)*aos_lapl_array(3,j)+aos_array(i)*aos_3rd_cross_array(5,j)

   nabla_3_at_r(6,i,j) = aos_3rd_cross_array(6,i)*aos_array(j)+aos_lapl_array(1,i)*aos_grad_array(3,j)+2.d0*aos_2nd_cross_array(3,i)*aos_grad_array(1,j)+2.d0*aos_grad_array(1,i)*aos_2nd_cross_array(3,j)+aos_grad_array(3,i)*aos_lapl_array(1,j)+aos_array(i)*aos_3rd_cross_array(6,j)


   nabla_3_at_r(7,i,j) = aos_3rd_direct_array(1,i)*aos_array(j)+ aos_array(i)*aos_3rd_direct_array(1,j)+ 3.d0*(aos_lapl_array(1,i)*aos_grad_array(1,j)+aos_grad_array(1,i)*aos_lapl_array(1,j)) 

   nabla_3_at_r(8,i,j) = aos_3rd_direct_array(2,i)*aos_array(j)+ aos_array(i)*aos_3rd_direct_array(2,j)+ 3.d0*(aos_lapl_array(2,i)*aos_grad_array(2,j)+aos_grad_array(2,i)*aos_lapl_array(2,j)) 

   nabla_3_at_r(9,i,j) = aos_3rd_direct_array(3,i)*aos_array(j)+ aos_array(i)*aos_3rd_direct_array(3,j)+ 3.d0*(aos_lapl_array(3,i)*aos_grad_array(3,j)+aos_grad_array(3,i)*aos_lapl_array(3,j)) 
  enddo
 enddo

 end


 subroutine give_nabla_4_at_contributions(r,nabla_4_at_r_contrib)
 implicit none
 BEGIN_DOC
! input      : r(1) ==> r(1) = x, r(2) = y, r(3) = z
!
! output     : aos_array(i) = ao(i) evaluated at r
!
!            : blablablablablabla 
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: nabla_4_at_r_contrib(6,ao_num,ao_num) 

 integer :: i,j
 double precision :: d2xd2y,d2yd2z,d2xd2z
 double precision :: d4_x,d4_y,d4_z

 double precision :: aos_array(ao_num)
 double precision :: aos_grad_array(3,ao_num)

 double precision :: aos_lapl_array(3,ao_num)
 double precision :: aos_3rd_direct_array(3,ao_num)
 double precision :: aos_4th_direct_array(3,ao_num)

 double precision :: aos_2nd_cross_array(3,ao_num)
 double precision :: aos_3rd_cross_array(6,ao_num)
 double precision :: aos_4th_cross_array(3,ao_num)

 call give_all_aos_and_fourth_at_r(r,aos_array,aos_grad_array,aos_lapl_array,aos_3rd_direct_array,aos_4th_direct_array)
 call give_all_aos_and_fourth_order_cross_terms_at_r(r,aos_array,aos_grad_array,aos_2nd_cross_array,aos_3rd_cross_array,aos_4th_cross_array)
 call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)
 
 do i = 1,ao_num
  do j = 1,ao_num
  
   d2xd2y = aos_4th_cross_array(1,i)*aos_array(j) + aos_lapl_array(1,i)*aos_lapl_array(2,j)+ aos_lapl_array(2,i)*aos_lapl_array(1,j)+aos_array(i)*aos_4th_cross_array(1,j) + 2.d0*(aos_3rd_cross_array(1,i)*aos_grad_array(2,j) + aos_3rd_cross_array(4,i)*aos_grad_array(1,j)+ aos_grad_array(1,i)*aos_3rd_cross_array(4,j) + aos_grad_array(2,i)*aos_3rd_cross_array(1,j)+ 2d0* aos_2nd_cross_array(1,i)*aos_2nd_cross_array(1,j)) 

   d2yd2z = aos_4th_cross_array(2,i)*aos_array(j) + aos_lapl_array(2,i)*aos_lapl_array(3,j)+ aos_lapl_array(3,i)*aos_lapl_array(2,j)+aos_array(i)*aos_4th_cross_array(2,j)+ 2d0*(aos_3rd_cross_array(2,i)*aos_grad_array(3,j) + aos_3rd_cross_array(5,i)*aos_grad_array(2,j)+ aos_grad_array(2,i)*aos_3rd_cross_array(5,j)+ aos_grad_array(3,i)*aos_3rd_cross_array(2,j)+ 2d0*aos_2nd_cross_array(2,i)*aos_2nd_cross_array(2,j))

   d2xd2z = aos_4th_cross_array(3,i)*aos_array(j) + aos_lapl_array(1,i)*aos_lapl_array(3,j)+ aos_lapl_array(3,i)*aos_lapl_array(1,j)+aos_array(i)*aos_4th_cross_array(3,j)+ 2d0*(aos_3rd_cross_array(6,i)*aos_grad_array(3,j) + aos_3rd_cross_array(3,i)*aos_grad_array(1,j)+ aos_grad_array(1,i)*aos_3rd_cross_array(3,j)+ aos_grad_array(3,i)*aos_3rd_cross_array(6,j)+ 2d0* aos_2nd_cross_array(3,i)*aos_2nd_cross_array(3,j))
 

   d4_x = aos_4th_direct_array(1,i)*aos_array(j)+aos_array(i)*aos_4th_direct_array(1,j)+ 2d0*(2d0*aos_3rd_direct_array(1,i)*aos_grad_array(1,j)+2d0*aos_grad_array(1,i)*aos_3rd_direct_array(1,j)+ 3d0*aos_lapl_array(1,i)*aos_lapl_array(1,j) )
   
   d4_y = aos_4th_direct_array(2,i)*aos_array(j)+aos_array(i)*aos_4th_direct_array(2,j)+ 2d0*(2d0*aos_3rd_direct_array(2,i)*aos_grad_array(2,j)+2d0*aos_grad_array(2,i)*aos_3rd_direct_array(2,j)+ 3d0*aos_lapl_array(2,i)*aos_lapl_array(2,j) ) 

   d4_z = aos_4th_direct_array(3,i)*aos_array(j)+aos_array(i)*aos_4th_direct_array(3,j)+ 2d0*(2d0*aos_3rd_direct_array(3,i)*aos_grad_array(3,j)+2d0*aos_grad_array(3,i)*aos_3rd_direct_array(3,j)+ 3d0*aos_lapl_array(3,i)*aos_lapl_array(3,j) )

   nabla_4_at_r_contrib(1,i,j) = d2xd2y 
   nabla_4_at_r_contrib(2,i,j) = d2yd2z
   nabla_4_at_r_contrib(3,i,j) = d2xd2z
   nabla_4_at_r_contrib(4,i,j) = d4_x
   nabla_4_at_r_contrib(5,i,j) = d4_y 
   nabla_4_at_r_contrib(6,i,j) = d4_z 


  enddo
 enddo
 end

 BEGIN_PROVIDER[double precision, aos_nabla_4_in_r_array, (ao_num,ao_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, aos_nabla_4_in_r_array_transp, (n_points_final_grid,ao_num,ao_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k
 double precision :: r(3)
 double precision :: nabla_4_at_r(ao_num,ao_num)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_nabla_4_at_r(r,nabla_4_at_r) 
  do j = 1, ao_num
   do k=1,ao_num
    aos_nabla_4_in_r_array(k,j,i) = nabla_4_at_r(k,j)
   !aos_nabla_4_in_r_array_transp(i,k,j) = nabla_4_at_r(k,j) 
   enddo
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, aos_nabla_2_in_r_array, (ao_num,ao_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, aos_nabla_2_in_r_array_transp, (n_points_final_grid,ao_num,ao_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k
 double precision :: r(3)
 double precision :: nabla_2_at_r(3,ao_num,ao_num)
 double precision :: nabla_2_tot_at_r(ao_num,ao_num)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_nabla_2_at_r(r,nabla_2_at_r,nabla_2_tot_at_r)
  do j = 1, ao_num
   do k=1,ao_num
    aos_nabla_2_in_r_array(k,j,i) = nabla_2_tot_at_r(k,j)
   !aos_nabla_2_in_r_array_transp(i,k,j) = nabla_2_tot_at_r(k,j) 
   enddo
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_nabla_4_in_r_array, (mo_num,mo_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, mos_nabla_4_in_r_array_transp, (n_points_final_grid,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,m,n
 double precision :: r(3)
 double precision :: nabla_4_at_r(ao_num,ao_num)
 mos_nabla_4_in_r_array =0.d0
 mos_nabla_4_in_r_array_transp =0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num 
   do k = 1, mo_num
    do m = 1, ao_num
     do n = 1, ao_num 
      mos_nabla_4_in_r_array(k,j,i) += mo_coef(m,j)*mo_coef(n,k)* aos_nabla_4_in_r_array(n,m,i)  
!     mos_nabla_4_in_r_array_transp(i,k,j) += mo_coef(m,j)*mo_coef(n,k)*aos_nabla_4_in_r_array(n,m,i) 
     enddo
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER

 subroutine give_nabla_2_at_r_mo(r,nabla_2_at_r_mo)
 implicit none
 BEGIN_DOC
! blablablabla
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: nabla_2_at_r_mo(mo_num,mo_num)

 double precision :: nabla_2_at_r(3,ao_num,ao_num)
 double precision :: nabla_2_tot_at_r(ao_num,ao_num)
 integer :: j,k,m,n
 call give_nabla_2_at_r(r,nabla_2_at_r,nabla_2_tot_at_r)

 nabla_2_at_r_mo = 0.d0

 do j = 1, mo_num 
  do k = 1, mo_num
   do m = 1, ao_num
    do n = 1, ao_num 
     nabla_2_at_r_mo(k,j) += mo_coef(m,j)*mo_coef(n,k)* nabla_2_tot_at_r(n,m) 
    enddo
   enddo
  enddo
 enddo

 end

 subroutine give_nabla_2_at_r_act_mo(r,nabla_2_at_r_mo)
 implicit none
 BEGIN_DOC
! blablablabla
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: nabla_2_at_r_mo(n_act_orb,n_act_orb)

 double precision :: nabla_2_at_r(3,ao_num,ao_num)
 double precision :: nabla_2_tot_at_r(ao_num,ao_num)
 integer :: j,k,m,n
 call give_nabla_2_at_r(r,nabla_2_at_r,nabla_2_tot_at_r)

 nabla_2_at_r_mo = 0.d0

 do j = 1, n_act_orb
  do k = 1, n_act_orb
   do m = 1, ao_num
    do n = 1, ao_num 
     nabla_2_at_r_mo(k,j) += mo_coef(m,list_act(j))*mo_coef(n,list_act(k))* nabla_2_tot_at_r(n,m) 
    enddo
   enddo
  enddo
 enddo

 end


 subroutine give_nabla_4_at_r_mo(r,nabla_4_at_r_mo)
 implicit none
 BEGIN_DOC
! blablablabla
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: nabla_4_at_r_mo(mo_num,mo_num)

 double precision :: nabla_4_at_r(ao_num,ao_num)
 integer :: j,k,m,n
 call give_nabla_4_at_r(r,nabla_4_at_r)

 nabla_4_at_r_mo = 0.d0

 do j = 1, mo_num 
  do k = 1, mo_num
   do m = 1, ao_num
    do n = 1, ao_num 
     nabla_4_at_r_mo(k,j) += mo_coef(m,j)*mo_coef(n,k)* nabla_4_at_r(n,m) 
    enddo
   enddo
  enddo
 enddo

 end


 BEGIN_PROVIDER[double precision, mos_nabla_4_in_r_array_2, (mo_num,mo_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, mos_nabla_4_in_r_array_transp_2, (n_points_final_grid,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,m,n
 double precision :: tempo_matrix(mo_num,ao_num)
 mos_nabla_4_in_r_array_transp_2= 0.d0
 do i = 1, n_points_final_grid
  call dgemm('N','N',mo_num,ao_num,ao_num,1.d0,mo_coef_transp,ao_num,aos_nabla_4_in_r_array(1,1,i),ao_num,0.d0,tempo_matrix,mo_num)
  call dgemm('N','N',mo_num,mo_num,ao_num,1.d0,tempo_matrix,mo_num,mo_coef,ao_num,0.d0,mos_nabla_4_in_r_array_2(1,1,i),mo_num)
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_nabla_2_in_r_array_2, (mo_num,mo_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, mos_nabla_2_in_r_array_transp_2, (n_points_final_grid,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,m,n
 double precision :: tempo_matrix(mo_num,ao_num)
 mos_nabla_2_in_r_array_transp_2= 0.d0
 do i = 1, n_points_final_grid
  call dgemm('N','N',mo_num,ao_num,ao_num,1.d0,mo_coef_transp,ao_num,aos_nabla_2_in_r_array(1,1,i),ao_num,0.d0,tempo_matrix,mo_num)
  call dgemm('N','N',mo_num,mo_num,ao_num,1.d0,tempo_matrix,mo_num,mo_coef,ao_num,0.d0,mos_nabla_2_in_r_array_2(1,1,i),mo_num)
 enddo
 END_PROVIDER

