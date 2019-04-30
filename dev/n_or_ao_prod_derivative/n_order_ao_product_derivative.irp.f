program n_order_ao_product_derivative
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
! call test_ao
! call test_grad_ao
! call test_grad_lapl_ao
! call test_fourth_ao
! call test_fourth_ao_cross
! call test_lapl_ao_product
  call test_nabla_4_ao_product
end

subroutine test_ao
 implicit none
 integer :: i,j,k,l,m,n
 integer :: order_der(3)
 integer :: num_ao,power_ao(3)
 double precision :: r(3),center_ao(3)
 double precision :: accu,dx,dy,dz,r2
 double precision :: aos_array(ao_num)
 double precision :: grad_aos_array(3,ao_num)
 double precision :: ao_part_derivative,pol_part_derivative,exp_part_derivative
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test ao'
 accu = 0d0
 order_der=0
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
   call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)
   do j = 1, ao_num
    
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)

    accu += dabs(aos_array(j) - ao_part_derivative) * final_weight_at_r_vector(i)
  enddo

 enddo
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'accu',accu
 print*,'\\\\\\\\\\\\\\\\\'
end

subroutine test_grad_ao
 implicit none
 integer :: i,j,k,l,m,n
 integer :: order_der(3)
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu(3),accu_2(3),accu_3(3)
 double precision :: grad_aos_array_bis(3,ao_num)
 double precision :: aos_array_plus(ao_num),aos_array_minus(ao_num)
 double precision :: lapl_aos_array_bis(3,ao_num),aos_lapl(3,ao_num)
 double precision :: grad_aos_array_plus(3,ao_num),grad_aos_array_minus(3,ao_num)
 double precision :: dr,ao_part_derivative
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test MO'
 print*,'dr,error grad, error lapl'
 do n = 1, 16
!do n = 4, 8
  dr = 10d0**(-n)
  r = 0d0
  accu = 0d0
  accu_2= 0d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   do m = 1,3
    rdx_plus = r
    rdx_plus(m) = r(m) + dr
    rdx_minus = r
    rdx_minus(m) = r(m) - dr
    call give_all_aos_and_grad_at_r(rdx_plus,aos_array_plus,grad_aos_array_plus)
    call give_all_aos_and_grad_at_r(rdx_minus,aos_array_minus,grad_aos_array_minus)
    do j = 1, ao_num
     grad_aos_array_bis(m,j) = (aos_array_plus(j) - aos_array_minus(j))/(2.d0 * dr)

     if (m .eq. 1)then
      order_der(1)=1
      order_der(2)=0
      order_der(3)=0
    !print*,'order1,2,3=  ',order_der(1),order_der(2),order_der(3) 
     ELSE IF (m .eq. 2) then 
      order_der(1)=0
      order_der(2)=1
      order_der(3)=0
     ELSE IF (m .eq. 3) then
      order_der(1)=0
      order_der(2)=0
      order_der(3)=1
     endif
     
     call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)

     accu(m) += dabs(grad_aos_array_bis(m,j) - ao_part_derivative) * final_weight_at_r_vector(i)
    !if (dabs(grad_aos_array_bis(m,j) - ao_part_derivative) .ge. 1.0d-5)then
    ! print*,'dabs diff int num ana  ', dabs(grad_aos_array_bis(m,j) - ao_part_derivative)
    ! print*,'grad ao num  ', grad_aos_array_bis(m,j)
    ! print*,'grad ao ana  ',ao_part_derivative
    ! print*,'r=',r(1),r(2),r(3)
    ! print*,'l value of the ao  ', ao_l(j)
    ! print*,'power ao', ao_power(j,1:3)
    ! print*,'center ao', nucl_coord(j,1:3)
    ! pause
    !endif
    !lapl_mos_array_bis(m,j) = (grad_mos_array_plus(m,j) -grad_mos_array_minus(m,j))/(2.d0 * dr)
    !accu_3(m) += dabs(lapl_mos_array_bis(m,j) - mos_lapl_in_r_array(j,i,m)) * final_weight_at_r_vector(i)
    enddo
   enddo

  enddo
  write(*,'(5(F16.10,x))')dr,accu(1),accu(2),accu(3)
 !write(33,'(100(F16.10,X))'),dr,accu(1),accu(2),accu(3),accu_2(1),accu_2(2),accu_2(3)
 enddo
end

subroutine test_grad_lapl_ao
 implicit none
 integer :: i,j,k,l,m,n
 integer :: order_der(3)
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu(3),accu_2(3),accu_3(3)
 double precision :: grad_aos_array_bis(3,ao_num)
 double precision :: aos_array_plus(ao_num),aos_array_minus(ao_num),aos_array(ao_num)
 double precision :: lapl_aos_array_bis(3,ao_num),aos_lapl(3,ao_num)
 double precision :: grad_aos_array_plus(3,ao_num),grad_aos_array_minus(3,ao_num),grad_aos_array(3,ao_num)
 double precision :: dr,ao_part_derivative
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test MO'
 print*,'dr,error grad, error lapl'
 print*,'yooooolllloooooo',0.d0**0
 do n = 1, 16
  dr = 10d0**(-n)
  r = 0d0
  accu = 0d0
  accu_2= 0d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   do m = 1,3
    rdx_plus = r
    rdx_plus(m) = r(m) + dr
    rdx_minus = r
    rdx_minus(m) = r(m) - dr
    call give_all_aos_and_grad_at_r(rdx_plus,aos_array_plus,grad_aos_array_plus)
    call give_all_aos_and_grad_at_r(rdx_minus,aos_array_minus,grad_aos_array_minus)

    call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,grad_aos_array,aos_lapl)
    do j = 1, ao_num
     grad_aos_array_bis(m,j) = (aos_array_plus(j) - aos_array_minus(j))/(2.d0 * dr)

     if (m .eq. 1)then
      order_der(1)=1
      order_der(2)=0
      order_der(3)=0
     ELSE IF (m .eq. 2) then 
      order_der(1)=0
      order_der(2)=1
      order_der(3)=0
     ELSE IF (m .eq. 3) then
      order_der(1)=0
      order_der(2)=0
      order_der(3)=1
     endif
     
     call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)

     accu(m) += dabs(grad_aos_array_bis(m,j) - ao_part_derivative) * final_weight_at_r_vector(i)
    !if (dabs(grad_aos_array_bis(m,j) - ao_part_derivative) .ge. 1.0d-5)then
    ! print*,'MERDE  ', dabs(grad_aos_array_bis(m,j) - ao_part_derivative)
    !endif

     if (m .eq. 1)then
      order_der(1)=2
      order_der(2)=0
      order_der(3)=0
     ELSE IF (m .eq. 2) then 
      order_der(1)=0
      order_der(2)=2
      order_der(3)=0
     ELSE IF (m .eq. 3) then
      order_der(1)=0
      order_der(2)=0
      order_der(3)=2
     endif
    !print*,'yooollllooooo '
     call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)

     lapl_aos_array_bis(m,j) = (grad_aos_array_plus(m,j) -grad_aos_array_minus(m,j))/(2.d0 * dr)
     accu_2(m) += dabs(lapl_aos_array_bis(m,j) - ao_part_derivative) * final_weight_at_r_vector(i)
       
    !if (dabs(lapl_aos_array_bis(m,j) - ao_part_derivative) * final_weight_at_r_vector(i) .ge. 1.0d-5)then
    ! print*,'dabs diff int num ana  ', dabs(lapl_aos_array_bis(m,j) - ao_part_derivative)!* final_weight_at_r_vector(i) 
    ! print*,'grad ao num  ', lapl_aos_array_bis(m,j)
    ! print*,'grad ao ana  ',ao_part_derivative
    ! print*,'grad ao ana juste ',aos_lapl(m,j)
    ! print*,'r=',r(1),r(2),r(3)
    ! print*,'l value of the ao  ', ao_l(j)
    ! print*,'power ao', ao_power(j,1:3)
    ! print*,'center ao', nucl_coord(j,1:3)
    ! print*,'weight',final_weight_at_r_vector(i) 
    ! pause
    !endif
      
    enddo
   enddo

  enddo
  
  write(*,'(5(F16.10,x))')dr,accu(1),accu_2(1)
  write(33,'(100(F16.10,X))'),dr,accu(1),accu(2),accu(3),accu_2(1),accu_2(2),accu_2(3)
 enddo
end



subroutine test_fourth_ao
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu(3),accu_2(3),accu_3(3),accu_4(3)
 double precision :: grad_aos_array_bis(3,ao_num)
 double precision :: aos_array_plus(ao_num),aos_array_minus(ao_num)
 double precision :: lapl_aos_array_bis(3,ao_num),aos_3rd_alpha_bis(3,ao_num),aos_4th_alpha_bis(3,ao_num)
 double precision :: grad_aos_array_plus(3,ao_num),grad_aos_array_minus(3,ao_num),aos_3rd_alpha_array_minus(3,ao_num),aos_4th_alpha_array_minus(3,ao_num)
 double precision :: aos_lapl_array_plus(3,ao_num),aos_lapl_array_minus(3,ao_num),aos_3rd_alpha_array_plus(3,ao_num),aos_4th_alpha_array_plus(3,ao_num)
 double precision :: dr
 double precision :: aos_array(3,ao_num),aos_grad_array(3,ao_num),aos_lapl_array(3,ao_num),aos_3rd_alpha_array(3,ao_num),aos_4th_alpha_array(3,ao_num)
 integer :: order_der(3)
 double precision :: ao_part_derivative
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test AO AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
 print*,'dr,error grad, error lapl'
 do n = 1, 16
  dr = 10d0**(-n)
  r = 0d0
  accu = 0d0
  accu_2= 0d0
  accu_3= 0d0
  accu_4= 0d0
  do i = 1, n_points_final_grid
!  print*,'point num =',i
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   do m = 1,3
    rdx_plus = r
    rdx_plus(m) = r(m) + dr
    rdx_minus = r
    rdx_minus(m) = r(m) - dr
    call give_all_aos_and_grad_and_lapl_at_r(rdx_plus,aos_array_plus,grad_aos_array_plus,aos_lapl_array_plus)
    call give_all_aos_and_grad_and_lapl_at_r(rdx_minus,aos_array_minus,grad_aos_array_minus,aos_lapl_array_minus)

    call give_all_aos_and_fourth_at_r(rdx_plus,aos_array_plus,grad_aos_array_plus,aos_lapl_array_plus,aos_3rd_alpha_array_plus,aos_4th_alpha_array_plus)
    call give_all_aos_and_fourth_at_r(rdx_minus,aos_array_minus,grad_aos_array_minus,aos_lapl_array_minus,aos_3rd_alpha_array_minus,aos_4th_alpha_array_minus)

    call  give_all_aos_and_fourth_at_r(r,aos_array,aos_grad_array,aos_lapl_array,aos_3rd_alpha_array,aos_4th_alpha_array)

    do j = 1, ao_num
!    print*,'aonum =',j
    !grad_aos_array_bis(m,j) = (aos_array_plus(j) - aos_array_minus(j))/(2.d0 *dr)
    !accu(m) += dabs(grad_aos_array_bis(m,j) - aos_grad_array(m,j)) * final_weight_at_r_vector(i)

    !lapl_aos_array_bis(m,j) = (grad_aos_array_plus(m,j) - grad_aos_array_minus(m,j))/(2.d0 * dr)
    !accu_2(m) += dabs(lapl_aos_array_bis(m,j) - aos_lapl_array(m,j)) * final_weight_at_r_vector(i)

     if (m .eq. 1)then
      order_der(1)=3
      order_der(2)=0
      order_der(3)=0
     ELSE IF (m .eq. 2) then
      order_der(1)=0
      order_der(2)=3
      order_der(3)=0
     ELSE IF (m .eq. 3) then
      order_der(1)=0
      order_der(2)=0
      order_der(3)=3
     endif
     
!    print*,'order_der',order_der
!    print*,'yoooollllooo1'
     call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
!    print*,'yoooollllooo2'

     aos_3rd_alpha_bis(m,j) = (aos_lapl_array_plus(m,j) - aos_lapl_array_minus(m,j))/(2.d0 * dr)
     accu_3(m) += dabs(aos_3rd_alpha_bis(m,j) - ao_part_derivative) * final_weight_at_r_vector(i)
     if (m .eq. 1)then
      order_der(1)=4
      order_der(2)=0
      order_der(3)=0
     ELSE IF (m .eq. 2) then
      order_der(1)=0
      order_der(2)=4
      order_der(3)=0
     ELSE IF (m .eq. 3) then
      order_der(1)=0
      order_der(2)=0
      order_der(3)=4
     endif

     call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)

     aos_4th_alpha_bis(m,j) = (aos_3rd_alpha_array_plus(m,j) - aos_3rd_alpha_array_minus(m,j))/(2.d0 * dr)
     accu_4(m) += dabs(aos_4th_alpha_bis(m,j) - ao_part_derivative) * final_weight_at_r_vector(i)
   
    enddo
   enddo
  enddo
  print*,dr,accu_3(1),accu_4(1)
  
  write(22,'(100(F16.10,X))'),dr,accu(1),accu(2),accu(3)
  write(33,'(100(F16.10,X))'),dr,accu_2(1),accu_2(2),accu_2(3)
  write(44,'(100(F16.10,X))'),dr,accu_3(1),accu_3(2),accu_3(3)
  write(55,'(100(F16.10,X))'),dr,accu_4(1),accu_4(2),accu_4(3)
 enddo
end


subroutine test_fourth_ao_cross
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rdx_plus(3),rdx_minus(3),rdy_plus(3),rdy_minus(3),rdz_plus(3),rdz_minus(3),accu(3),accu_2(3),accu_3(6),accu_4(3)
 double precision :: grad_aos_array_bis(3,ao_num)
 double precision :: aos_array_plus(ao_num),aos_array_minus(ao_num)
 double precision :: lapl_aos_array_bis(3,ao_num),aos_3rd_cross_bis(6,ao_num),aos_4th_alpha_bis(3,ao_num)
 double precision :: grad_aos_array_plus_y(3,ao_num),grad_aos_array_minus_y(3,ao_num),aos_3rd_alpha_array_minus_y(6,ao_num),aos_4th_alpha_array_minus_y(3,ao_num)
 double precision :: aos_lapl_array_plus_y(3,ao_num),aos_lapl_array_minus_y(3,ao_num),aos_3rd_alpha_array_plus_y(6,ao_num),aos_4th_alpha_array_plus_y(3,ao_num)


 double precision :: aos_array_plus_x(3,ao_num),aos_array_plus_y(3,ao_num),aos_array_plus_z(3,ao_num)
 double precision :: aos_array_minus_x(3,ao_num),aos_array_minus_y(3,ao_num),aos_array_minus_z(3,ao_num)
 double precision :: grad_aos_array_plus_x(3,ao_num),grad_aos_array_minus_x(3,ao_num),aos_3rd_alpha_array_minus_x(6,ao_num),aos_4th_alpha_array_minus_x(3,ao_num)
 double precision :: aos_lapl_array_plus_x(3,ao_num),aos_lapl_array_minus_x(3,ao_num),aos_3rd_alpha_array_plus_x(6,ao_num),aos_4th_alpha_array_plus_x(3,ao_num)
 double precision :: grad_aos_array_plus_z_2(3,ao_num),grad_aos_array_minus_z_2(3,ao_num)
 double precision :: grad_aos_array_plus_y_2(3,ao_num),grad_aos_array_minus_y_2(3,ao_num)
 double precision :: grad_aos_array_plus_x_2(3,ao_num),grad_aos_array_minus_x_2(3,ao_num)
 double precision :: aos_lapl_array_plus_z_2(3,ao_num),aos_lapl_array_minus_z_2(3,ao_num)
 double precision :: aos_lapl_array_plus_y_2(3,ao_num),aos_lapl_array_minus_y_2(3,ao_num)
 double precision :: aos_lapl_array_plus_x_2(3,ao_num),aos_lapl_array_minus_x_2(3,ao_num)


 double precision :: grad_aos_array_plus_z(3,ao_num),grad_aos_array_minus_z(3,ao_num),aos_3rd_alpha_array_minus_z(6,ao_num),aos_4th_alpha_array_minus_z(3,ao_num)
 double precision :: aos_lapl_array_plus_z(3,ao_num),aos_lapl_array_minus_z(3,ao_num),aos_3rd_alpha_array_plus_z(6,ao_num),aos_4th_alpha_array_plus_z(3,ao_num)

 double precision :: dr
 double precision :: aos_array(3,ao_num),aos_grad_array(3,ao_num),aos_lapl_array(3,ao_num),aos_3rd_cross_array(6,ao_num),aos_4th_cross_array(3,ao_num)

 integer :: order_der(3)
 double precision :: ao_part_derivative


 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test AO AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
 print*,'dr,error grad, error lapl'
 do n = 1, 16
  dr = 10d0**(-n)
  r = 0d0
  accu = 0d0
  accu_2= 0d0
  accu_3= 0d0
  accu_4= 0d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   rdx_plus = r
   rdx_plus(1) = r(1) + dr
   rdx_minus = r
   rdx_minus(1) = r(1) - dr

   rdy_plus = r
   rdy_plus(2) = r(2) + dr
   rdy_minus = r
   rdy_minus(2) = r(2) - dr


   rdz_plus = r
   rdz_plus(3) = r(3) + dr
   rdz_minus = r
   rdz_minus(3) = r(3) - dr

   call give_all_aos_and_grad_and_lapl_at_r(rdx_plus,aos_array_plus_x,grad_aos_array_plus_x,aos_lapl_array_plus_x)
   call give_all_aos_and_grad_and_lapl_at_r(rdx_minus,aos_array_minus_x,grad_aos_array_minus_x,aos_lapl_array_minus_x)

   call give_all_aos_and_grad_and_lapl_at_r(rdy_plus,aos_array_plus_y,grad_aos_array_plus_y,aos_lapl_array_plus_y)
   call give_all_aos_and_grad_and_lapl_at_r(rdy_minus,aos_array_minus_y,grad_aos_array_minus_y,aos_lapl_array_minus_y)

   call give_all_aos_and_grad_and_lapl_at_r(rdz_plus,aos_array_plus_z,grad_aos_array_plus_z,aos_lapl_array_plus_z)
   call give_all_aos_and_grad_and_lapl_at_r(rdz_minus,aos_array_minus_z,grad_aos_array_minus_z,aos_lapl_array_minus_z)

   call give_all_aos_and_fourth_order_cross_terms_at_r(rdx_plus,aos_array_plus_x,grad_aos_array_plus_x_2,aos_lapl_array_plus_x_2,aos_3rd_alpha_array_plus_x,aos_4th_alpha_array_plus_x)
   call give_all_aos_and_fourth_order_cross_terms_at_r(rdx_minus,aos_array_minus_x,grad_aos_array_minus_x_2,aos_lapl_array_minus_x_2,aos_3rd_alpha_array_minus_x,aos_4th_alpha_array_minus_x)

   call give_all_aos_and_fourth_order_cross_terms_at_r(rdy_plus,aos_array_plus_y,grad_aos_array_plus_y_2,aos_lapl_array_plus_y_2,aos_3rd_alpha_array_plus_y,aos_4th_alpha_array_plus_y)
   call give_all_aos_and_fourth_order_cross_terms_at_r(rdy_minus,aos_array_minus_y,grad_aos_array_minus_y_2,aos_lapl_array_minus_y_2,aos_3rd_alpha_array_minus_y,aos_4th_alpha_array_minus_y)

   call give_all_aos_and_fourth_order_cross_terms_at_r(rdz_plus,aos_array_plus_z,grad_aos_array_plus_z_2,aos_lapl_array_plus_z_2,aos_3rd_alpha_array_plus_z,aos_4th_alpha_array_plus_z)
   call give_all_aos_and_fourth_order_cross_terms_at_r(rdz_minus,aos_array_minus_z,grad_aos_array_minus_z_2,aos_lapl_array_minus_z_2,aos_3rd_alpha_array_minus_z,aos_4th_alpha_array_minus_z)




   call give_all_aos_and_fourth_order_cross_terms_at_r(r,aos_array,aos_grad_array,aos_lapl_array,aos_3rd_cross_array,aos_4th_cross_array) 

   do j = 1, ao_num

    order_der(1)=1
    order_der(2)=1
    order_der(3)=0
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    lapl_aos_array_bis(1,j) = (grad_aos_array_plus_x(2,j) - grad_aos_array_minus_x(2,j))/(2.d0 * dr)
    accu_2(1) += dabs(lapl_aos_array_bis(1,j) - ao_part_derivative) * final_weight_at_r_vector(i)

   !order_der(1)=0
   !order_der(2)=1
   !order_der(3)=1
   !call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
   !lapl_aos_array_bis(2,j) = (grad_aos_array_plus_y(3,j) -grad_aos_array_minus_y(3,j))/(2.d0 * dr)
   !accu_2(2) += dabs(lapl_aos_array_bis(2,j) - ao_part_derivative) *final_weight_at_r_vector(i)

   !order_der(1)=1
   !order_der(2)=0
   !order_der(3)=1
   !call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
   !lapl_aos_array_bis(3,j) = (grad_aos_array_plus_x(3,j) -grad_aos_array_minus_x(3,j))/(2.d0 * dr)
   !accu_2(3) += dabs(lapl_aos_array_bis(3,j) - ao_part_derivative) *final_weight_at_r_vector(i)

    order_der(1)=2
    order_der(2)=1
    order_der(3)=0
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_3rd_cross_bis(1,j) = (aos_lapl_array_plus_y(1,j) - aos_lapl_array_minus_y(1,j))/(2.d0 * dr)
    accu_3(1) += dabs(aos_3rd_cross_bis(1,j) - ao_part_derivative) * final_weight_at_r_vector(i)

    !if (dabs(aos_3rd_cross_bis(1,j) - ao_part_derivative) * final_weight_at_r_vector(i) .ge. 1.0d-5)then
    ! print*,'dabs diff int num ana  ', dabs(aos_3rd_cross_bis(1,j) - ao_part_derivative)!* final_weight_at_r_vector(i) 
    ! print*,'grad ao num  ',aos_3rd_cross_bis(1,j) 
    ! print*,'grad ao ana  ',ao_part_derivative
    ! print*,'grad ao ana juste ',aos_3rd_cross_array(1,j)
    ! print*,'r=',r(1),r(2),r(3)
    ! print*,'l value of the ao  ', ao_l(j)
    ! print*,'power ao', ao_power(j,1:3)
    ! print*,'center ao', nucl_coord(j,1:3)
    ! print*,'weight',final_weight_at_r_vector(i) 
    ! pause
    !endif

    order_der(1)=0
    order_der(2)=2
    order_der(3)=1
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_3rd_cross_bis(2,j) = (aos_lapl_array_plus_z(2,j) - aos_lapl_array_minus_z(2,j))/(2.d0 * dr)
    accu_3(2) += dabs(aos_3rd_cross_bis(2,j) - ao_part_derivative) *final_weight_at_r_vector(i)

    order_der(1)=1
    order_der(2)=0
    order_der(3)=2
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_3rd_cross_bis(3,j) = (aos_lapl_array_plus_x(3,j) - aos_lapl_array_minus_x(3,j))/(2.d0 * dr)
    accu_3(3) += dabs(aos_3rd_cross_bis(3,j) - ao_part_derivative) * final_weight_at_r_vector(i)


    order_der(1)=1
    order_der(2)=2
    order_der(3)=0
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_3rd_cross_bis(4,j) = (aos_lapl_array_plus_x(2,j) - aos_lapl_array_minus_x(2,j))/(2.d0 * dr) 
    accu_3(4) += dabs(aos_3rd_cross_bis(4,j) - ao_part_derivative) * final_weight_at_r_vector(i)
   
    order_der(1)=0
    order_der(2)=1
    order_der(3)=2
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_3rd_cross_bis(5,j) = (aos_lapl_array_plus_y(3,j) - aos_lapl_array_minus_y(3,j))/(2.d0 * dr)
    accu_3(5) += dabs(aos_3rd_cross_bis(5,j) - ao_part_derivative) * final_weight_at_r_vector(i)

    order_der(1)=2
    order_der(2)=0
    order_der(3)=1
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_3rd_cross_bis(6,j) = (aos_lapl_array_plus_z(1,j) - aos_lapl_array_minus_z(1,j))/(2.d0 * dr) 
    accu_3(6) += dabs(aos_3rd_cross_bis(6,j) - ao_part_derivative) * final_weight_at_r_vector(i)

    order_der(1)=2
    order_der(2)=2
    order_der(3)=0
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_4th_alpha_bis(1,j) = (aos_3rd_alpha_array_plus_y(1,j) - aos_3rd_alpha_array_minus_y(1,j))/(2.d0 * dr)
    accu_4(1) += dabs(aos_4th_alpha_bis(1,j) - ao_part_derivative) * final_weight_at_r_vector(i)
   
    order_der(1)=0
    order_der(2)=2
    order_der(3)=2
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_4th_alpha_bis(2,j) = (aos_3rd_alpha_array_plus_z(2,j) - aos_3rd_alpha_array_minus_z(2,j))/(2.d0 * dr)
    accu_4(2) += dabs(aos_4th_alpha_bis(2,j) - ao_part_derivative) * final_weight_at_r_vector(i)

    order_der(1)=2
    order_der(2)=0
    order_der(3)=2
    call n_order_partial_derivative_ao(r,j,order_der,ao_part_derivative)
    aos_4th_alpha_bis(3,j) = (aos_3rd_alpha_array_plus_x(3,j) - aos_3rd_alpha_array_minus_x(3,j))/(2.d0 * dr)
    accu_4(3) += dabs(aos_4th_alpha_bis(3,j) - ao_part_derivative) * final_weight_at_r_vector(i)

   enddo
  enddo
  print*,dr,accu_2(1),accu_3(1),accu_4(1)
  
  !write(22,'(100(F16.10,X))'),dr,accu(1),accu(2),accu(3)
  write(33,'(100(F16.10,X))'),dr,accu_2(1),accu_2(2),accu_2(3)
  write(44,'(100(F16.10,X))'),dr,accu_3(1),accu_3(2),accu_3(3),accu_3(4),accu_3(5),accu_3(6)
  write(55,'(100(F16.10,X))'),dr,accu_4(1),accu_4(2),accu_4(3)
 enddo
end

 subroutine test_lapl_ao_product
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rdx_plus(3),rdx_minus(3),accu(3)
 double precision :: accu_2,lapl_num_tot
 double precision :: grad_aos_array_bis(3,ao_num)
 double precision :: aos_array(ao_num),aos_array_plus(ao_num),aos_array_minus(ao_num)
 double precision :: lapl_aos_product_array_bis(3,ao_num,ao_num)
 double precision :: grad_aos_array_plus(3,ao_num),grad_aos_array_minus(3,ao_num),grad_aos_array(3,ao_num)
 double precision :: dr
 double precision :: nabla_2_at_r(3,ao_num,ao_num)
 double precision :: nabla_2_tot_at_r(ao_num,ao_num)
 integer :: order_der(3)
 double precision :: ao_part_derivative,n_order_delta
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test AO'
 print*,'dr,error grad, error lapl'
 do n = 1, 16
  dr = 10d0**(-n)
  r = 0d0
  accu = 0d0
  accu_2 = 0d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   do m = 1,3
    rdx_plus = r
    rdx_plus(m) = r(m) + dr
    rdx_minus = r
    rdx_minus(m) = r(m) - dr


    !!!!! For Nabla 2

    call give_all_aos_and_grad_at_r(r,aos_array,grad_aos_array)

    call give_all_aos_and_grad_at_r(rdx_plus,aos_array_plus,grad_aos_array_plus)
    call give_all_aos_and_grad_at_r(rdx_minus,aos_array_minus,grad_aos_array_minus)

    call give_nabla_2_at_r(r,nabla_2_at_r,nabla_2_tot_at_r)

    if (m .eq. 1)then
     order_der(1)=2
     order_der(2)=0
     order_der(3)=0
    ELSE IF (m .eq. 2) then
     order_der(1)=0
     order_der(2)=2
     order_der(3)=0
    ELSE IF (m .eq. 3) then
     order_der(1)=0
     order_der(2)=0
     order_der(3)=2
    endif

    do k = 1, ao_num 
     do j = 1, ao_num
      call n_order_partial_derivative_ao_product(r,j,k,order_der,ao_part_derivative)
      lapl_aos_product_array_bis(m,k,j) = (aos_array_plus(k)*aos_array_plus(j)+aos_array_minus(k)*aos_array_minus(j)-2d0*aos_array(k)*aos_array(j) )/(dr**2)
      accu(m) += dabs(lapl_aos_product_array_bis(m,k,j) - ao_part_derivative ) * final_weight_at_r_vector(i)
     enddo
    enddo
   enddo
  
   

   do k = 1, ao_num
    do j = 1, ao_num
     call n_order_delta_ao_product(r,j,k,1,n_order_delta)
     lapl_num_tot = lapl_aos_product_array_bis(1,k,j)+ lapl_aos_product_array_bis(2,k,j) + lapl_aos_product_array_bis(3,k,j) 
     accu_2 += dabs(lapl_num_tot - n_order_delta) * final_weight_at_r_vector(i)
    enddo
   enddo


  enddo
  print*,dr,accu(1),accu_2
  write(33,'(100(F16.10,X))'),dr,accu(1),accu(2),accu(3)
 enddo
end



 subroutine test_nabla_4_ao_product
 implicit none
 integer :: i,j,k
 double precision :: r(3)
 double precision :: accu
 double precision :: nabla_2_tot_at_r(ao_num,ao_num)
 double precision :: n_order_delta
 double precision :: nabla_4_at_r(ao_num,ao_num)

 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '
 print*,'Test AO'
 print*,'dr,error grad, error lapl'
 accu = 0d0
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)

  call give_nabla_4_at_r(r,nabla_4_at_r)   

   do k = 1, ao_num 
    do j = 1, ao_num
     call n_order_delta_ao_product(r,j,k,2,n_order_delta)
     accu += dabs(nabla_4_at_r(j,k) - n_order_delta ) * final_weight_at_r_vector(i)
    enddo
   enddo
 
 enddo
 print*,'accu =',accu
end
