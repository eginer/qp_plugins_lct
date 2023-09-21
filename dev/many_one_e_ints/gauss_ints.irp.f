double precision function overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 double precision  :: overlap_x,overlap_y,overlap_z,overlap
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1
 dim1=100
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 accu = 0.d0
 do lx = 0, iorder_a_new(1)
  coefx = A_new(lx,1)
  if(dabs(coefx).lt.thr)cycle
  iorder_tmp(1) = lx
  do ly = 0, iorder_a_new(2)
   coefy = A_new(ly,2)
   coefxy = coefx * coefy 
   if(dabs(coefxy).lt.thr)cycle
   iorder_tmp(2) = ly
   do lz = 0, iorder_a_new(3)
    coefz = A_new(lz,3)
    coefxyz = coefxy * coefz 
    if(dabs(coefxyz).lt.thr)cycle
    iorder_tmp(3) = lz
    call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
    accu += coefxyz * overlap
   enddo
  enddo
 enddo
 overlap_gauss_r12 = fact_a_new * accu 
end


subroutine overlap_gauss_xyz_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,gauss_ints)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   gauss_ints(m) = \int dr exp(-delta (r - D)^2 ) * x/y/z (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  ! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)
 double precision, intent(out)   :: gauss_ints(3)

 double precision  :: overlap_x,overlap_y,overlap_z,overlap
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 integer           :: power_B_new(3)
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,m
 dim1=100
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 gauss_ints = 0.d0
 do lx = 0, iorder_a_new(1)
  coefx = A_new(lx,1)
  if(dabs(coefx).lt.thr)cycle
  iorder_tmp(1) = lx
  do ly = 0, iorder_a_new(2)
   coefy = A_new(ly,2)
   coefxy = coefx * coefy 
   if(dabs(coefxy).lt.thr)cycle
   iorder_tmp(2) = ly
   do lz = 0, iorder_a_new(3)
    coefz = A_new(lz,3)
    coefxyz = coefxy * coefz 
    if(dabs(coefxyz).lt.thr)cycle
    iorder_tmp(3) = lz
    do m = 1, 3
     ! change (x-Bx)^bx --> (x-Bx)^(bx+1) + Bx(x-Bx)^bx
     power_B_new = power_B
     power_B_new(m) += 1 ! (x-Bx)^(bx+1)
     call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
     gauss_ints(m) += coefxyz * overlap

     power_B_new = power_B
     call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
     gauss_ints(m) += coefxyz * overlap * B_center(m) ! Bx (x-Bx)^(bx)
    enddo
   enddo
  enddo
 enddo
 gauss_ints *= fact_a_new
end

subroutine overlap_gauss_xyz_r12_ao(D_center,delta,i,j,gauss_ints)
 implicit none
 BEGIN_DOC
! gauss_ints(m) = \int dr AO_i(r) AO_j(r) x/y/z e^{-delta |r-D_center|^2}
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in)  :: D_center(3), delta
 double precision, intent(out) :: gauss_ints(3)

 integer :: num_a,num_b,power_A(3), power_B(3),l,k,m
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,gauss_ints_tmp(3)
 gauss_ints = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)     
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)     
   call overlap_gauss_xyz_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,gauss_ints_tmp)
   do m = 1, 3
    gauss_ints(m) += gauss_ints_tmp(m) *  ao_coef_normalized_ordered_transp(l,i)             &
                                       *  ao_coef_normalized_ordered_transp(k,j) 
   enddo
  enddo 
 enddo
end



double precision function overlap_gauss_xyz_r12_specific(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,mx)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !    \int dr exp(-delta (r - D)^2 ) * x/y/z (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  ! with mx == 1 ==> x, mx == 2 ==> y, mx == 3 ==> z
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3),mx

 double precision  :: overlap_x,overlap_y,overlap_z,overlap
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 integer           :: power_B_new(3)
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,m
 dim1=100
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 overlap_gauss_xyz_r12_specific = 0.d0
 do lx = 0, iorder_a_new(1)
  coefx = A_new(lx,1)
  if(dabs(coefx).lt.thr)cycle
  iorder_tmp(1) = lx
  do ly = 0, iorder_a_new(2)
   coefy = A_new(ly,2)
   coefxy = coefx * coefy 
   if(dabs(coefxy).lt.thr)cycle
   iorder_tmp(2) = ly
   do lz = 0, iorder_a_new(3)
    coefz = A_new(lz,3)
    coefxyz = coefxy * coefz 
    if(dabs(coefxyz).lt.thr)cycle
    iorder_tmp(3) = lz
    m = mx
    ! change (x-Bx)^bx --> (x-Bx)^(bx+1) + Bx(x-Bx)^bx
    power_B_new = power_B
    power_B_new(m) += 1 ! (x-Bx)^(bx+1)
    call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
    overlap_gauss_xyz_r12_specific += coefxyz * overlap

    power_B_new = power_B
    call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
    overlap_gauss_xyz_r12_specific += coefxyz * overlap * B_center(m) ! Bx (x-Bx)^(bx)
   enddo
  enddo
 enddo
 overlap_gauss_xyz_r12_specific *= fact_a_new
end

double precision function overlap_gauss_xyz_r12_ao_specific(D_center,delta,i,j,mx)
 implicit none
 BEGIN_DOC
! \int dr AO_i(r) AO_j(r) x/y/z e^{-delta |r-D_center|^2}
!
! with mx == 1 ==> x, mx == 2 ==> y, mx == 3 ==> z
 END_DOC
 integer, intent(in) :: i,j,mx
 double precision, intent(in)  :: D_center(3), delta

 integer :: num_a,num_b,power_A(3), power_B(3),l,k
 double precision :: gauss_int
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta
 double precision :: overlap_gauss_xyz_r12_specific
 overlap_gauss_xyz_r12_ao_specific = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)     
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)     
   gauss_int = overlap_gauss_xyz_r12_specific(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,mx)
   overlap_gauss_xyz_r12_ao_specific = gauss_int *  ao_coef_normalized_ordered_transp(l,i)             &
                                                 *  ao_coef_normalized_ordered_transp(k,j) 
  enddo 
 enddo
end


subroutine overlap_gauss_r12_all_ao(D_center,delta,aos_ints)
 implicit none
 double precision, intent(in) :: D_center(3), delta
 double precision, intent(out):: aos_ints(ao_num,ao_num)

 integer :: num_a,num_b,power_A(3), power_B(3),l,k,i,j
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,analytical_j
 aos_ints = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   if(ao_overlap_abs(j,i).lt.1.d-12)cycle
   num_A = ao_nucl(i)
   power_A(1:3)= ao_power(i,1:3)
   A_center(1:3) = nucl_coord(num_A,1:3)
   num_B = ao_nucl(j)
   power_B(1:3)= ao_power(j,1:3)
   B_center(1:3) = nucl_coord(num_B,1:3)
   do l=1,ao_prim_num(i)
    alpha = ao_expo_ordered_transp(l,i)     
    do k=1,ao_prim_num(j)
     beta = ao_expo_ordered_transp(k,j)     
     analytical_j = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
     aos_ints(j,i) += analytical_j *  ao_coef_normalized_ordered_transp(l,i)             &
                                   *  ao_coef_normalized_ordered_transp(k,j) 
    enddo 
   enddo
  enddo
 enddo
end

double precision function overlap_gauss_r12_ao(D_center,delta,i,j)
 implicit none
 BEGIN_DOC
! \int dr AO_i(r) AO_j(r) e^{-delta |r-D_center|^2}
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: D_center(3), delta

 integer :: num_a,num_b,power_A(3), power_B(3),l,k
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,analytical_j
 overlap_gauss_r12_ao = 0.d0
 if(ao_overlap_abs(j,i).lt.1.d-12)then
  return
 endif
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)     
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)     
   analytical_j = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
   overlap_gauss_r12_ao += analytical_j *  ao_coef_normalized_ordered_transp(l,i)             &
                                        *  ao_coef_normalized_ordered_transp(k,j) 
  enddo 
 enddo
end


