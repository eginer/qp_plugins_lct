
BEGIN_PROVIDER [ double precision, three_body_3_index, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix element of the -L  three-body operator 
!
! three_body_3_index(k,l,n) = < phi_k phi_l phi_n | phi_k phi_l phi_n >
!
! notice the -1 sign: in this way three_body_3_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 name_file = 'six_index_tensor'
 provide x_W_ij_erf_rk
 three_body_3_index = 0.d0
 print*,'Providing the three_body_3_index ...'
 call wall_time(wall0)
 call wall_time(wall1)
 print*,'wall time for three_body_3_index',wall1 - wall0
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_body_3_index)
 !$OMP DO SCHEDULE (guided) COLLAPSE(3)
  do m = 1, mo_num ! 3
   do j = 1, mo_num ! 2 
    do i = 1, mo_num ! 1 
     integral = 0.d0
     !                          1 2 3 1 2 3
     call give_integrals_3_body(i,j,m,i,j,m,integral)

     three_body_3_index(i,j,m) = -1.d0 * integral 
  
    enddo
   enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_12, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_3_index_exch_12(k,l,n) = < phi_k phi_l phi_n | phi_l phi_k phi_n >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 name_file = 'six_index_tensor'
 provide x_W_ij_erf_rk
 three_body_3_index_exch_12 = 0.d0
 print*,'Providing the three_body_3_index_exch_12 ...'
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_body_3_index_exch_12)
 !$OMP DO SCHEDULE (guided) COLLAPSE(3)
  do m = 1, mo_num ! 3
   do j = 1, mo_num ! 2 
    do i = 1, mo_num ! 1 
     integral = 0.d0
     !                          1 2 3 1 2 3
     call give_integrals_3_body(i,j,m,j,i,m,integral)

     three_body_3_index_exch_12(i,j,m) = -1.d0 * integral 
  
    enddo
   enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'wall time for three_body_3_index_exch_12',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_23, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_3_index_exch_12(k,l,n) = < phi_k phi_l phi_n | phi_k phi_n phi_l >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 name_file = 'six_index_tensor'
 provide x_W_ij_erf_rk
 three_body_3_index_exch_23 = 0.d0
 print*,'Providing the three_body_3_index_exch_23 ...'
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_body_3_index_exch_23)
 !$OMP DO SCHEDULE (guided) COLLAPSE(3)
  do m = 1, mo_num ! 3
   do j = 1, mo_num ! 2 
    do i = 1, mo_num ! 1 
     integral = 0.d0
     !                          1 2 3 1 2 3
     call give_integrals_3_body(i,j,m,i,m,j,integral)

     three_body_3_index_exch_23(i,j,m) = -1.d0 * integral 
  
    enddo
   enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'wall time for three_body_3_index_exch_23',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_13, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_3_index_exch_12(k,l,n) = < phi_k phi_l phi_n | phi_k phi_n phi_l >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 name_file = 'six_index_tensor'
 provide x_W_ij_erf_rk
 three_body_3_index_exch_13 = 0.d0
 print*,'Providing the three_body_3_index_exch_13 ...'
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_body_3_index_exch_13)
 !$OMP DO SCHEDULE (guided)
  do m = 1, mo_num ! 3
   do j = 1, mo_num ! 2 
    do i = 1, mo_num ! 1 
     integral = 0.d0
     !                          1 2 3 1 2 3
     call give_integrals_3_body(i,j,m,m,j,i,integral)

     three_body_3_index_exch_13(i,j,m) = -1.d0 * integral 
  
    enddo
   enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'wall time for three_body_3_index_exch_13',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [ double precision, three_body_3_index_exch_231, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix element of the -L  three-body operator 
!
! three_body_3_index_exch_231(k,l,n) = < phi_k phi_l phi_n | phi_l phi_n phi_k >
!
! notice the -1 sign: in this way three_body_3_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 name_file = 'six_index_tensor'
 provide x_W_ij_erf_rk
 three_body_3_index_exch_231 = 0.d0
 print*,'Providing the three_body_3_index_231 ...'
 call wall_time(wall0)
 call wall_time(wall1)
 print*,'wall time for three_body_3_index',wall1 - wall0
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_body_3_index_exch_231)
 !$OMP DO SCHEDULE (guided) COLLAPSE(3)
  do m = 1, mo_num ! 3
   do j = 1, mo_num ! 2 
    do i = 1, mo_num ! 1 
     integral = 0.d0
     !                          1 2 3 1 2 3
     call give_integrals_3_body(i,j,m,j,m,i,integral)

     three_body_3_index_exch_231(i,j,m) = -1.d0 * integral 
  
    enddo
   enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_312, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix element of the -L  three-body operator 
!
! three_body_3_index(k,l,n) = < phi_k phi_l phi_n | phi_l phi_n phi_k >
!
! notice the -1 sign: in this way three_body_3_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 name_file = 'six_index_tensor'
 provide x_W_ij_erf_rk
 three_body_3_index_exch_312 = 0.d0
 print*,'Providing the three_body_3_index_312 ...'
 call wall_time(wall0)
 call wall_time(wall1)
 print*,'wall time for three_body_3_index',wall1 - wall0
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_body_3_index_exch_312)
 !$OMP DO SCHEDULE (guided) COLLAPSE(3)
  do m = 1, mo_num ! 3
   do j = 1, mo_num ! 2 
    do i = 1, mo_num ! 1 
     integral = 0.d0
     !                          1 2 3 1 2 3
     call give_integrals_3_body(i,j,m,m,i,j,integral)

     three_body_3_index_exch_312(i,j,m) = -1.d0 * integral 
  
    enddo
   enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER 
