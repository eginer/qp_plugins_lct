program debug_test_ao_tc

  implicit none
  !call main()
  call print_1eInteg()
  !call print_nuc()
  !call j1b_grd()

end


subroutine main()

  implicit none 
  integer :: j, l

  do j = 1, nucl_num
    print *, j, ao_nucl(j), j1b_gauss_pen(j)
  enddo

  do j = 1, ao_num
    print *, ao_nucl(j)
    print *, ao_power(j,1:3)
    print *, nucl_coord(ao_nucl(j),1:3)
    print *, ao_prim_num(j)
    do l = 1, ao_prim_num(j)
      print *, ao_expo_ordered_transp(l,j), ao_coef_normalized_ordered_transp(l,j)
    enddo
  enddo

  return
end subroutine main



subroutine print_1eInteg()
 
  implicit none
  integer :: i, j

  do j = 1, ao_num
    do i = 1, ao_num
      print *, j1b_gauss_hermI(i,j), j1b_gauss_hermII(i,j), j1b_gauss_nonherm(i,j)
    enddo
  enddo

  return
end subroutine print_1eInteg



subroutine print_nuc()

  implicit none 
  integer :: j, l

  do j = 1, nucl_num
    print *, j, j1b_gauss_pen(j), nucl_coord(j,1:3)
  enddo

  do j = 1, ao_num
    print *, ao_nucl(j), nucl_coord(ao_nucl(j),1:3)
  enddo

  return
end subroutine print_nuc



subroutine j1b_grd()

  implicit none

  integer          :: num_A, num_B
  integer          :: power_A(3), power_B(3)
  integer          :: i, j, k1, k2, l, m
  double precision :: alpha, beta, gama1, gama2
  double precision :: A_center(3), B_center(3), C_center1(3), C_center2(3)
  double precision :: c1, c
  double precision :: integ_grd

  integer          :: dim1
  double precision :: overlap_y, d_a_2, overlap_z, overlap

  double precision :: int_gauss_4G_deb

  PROVIDE j1b_gauss_pen

  ! --------------------------------------------------------------------------------
  ! -- Dummy call to provide everything
  dim1        = 100
  A_center(:) = 0.d0
  B_center(:) = 1.d0
  alpha       = 1.d0
  beta        = 0.1d0
  power_A(:)  = 1
  power_B(:)  = 0
  call overlap_gaussian_xyz( A_center, B_center, alpha, beta, power_A, power_B &
                           , overlap_y, d_a_2, overlap_z, overlap, dim1 )
  ! --------------------------------------------------------------------------------

  j = 1
  num_A         = ao_nucl(j)
  power_A(1:3)  = ao_power(j,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)

  i = 42
  num_B         = ao_nucl(i)
  power_B(1:3)  = ao_power(i,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)

  integ_grd = 0.d0
  do l = 1, ao_prim_num(j)
    alpha = ao_expo_ordered_transp(l,j)

    do m = 1, ao_prim_num(i)
      beta = ao_expo_ordered_transp(m,i)

      c = 0.d0
      do k1 = 1, nucl_num
        gama1          = j1b_gauss_pen(k1)
        C_center1(1:3) = nucl_coord(k1,1:3)

        do k2 = 1, nucl_num
          gama2          = j1b_gauss_pen(k2)
          C_center2(1:3) = nucl_coord(k2,1:3)

          c1 = int_gauss_4G_deb( A_center, B_center, C_center1, C_center2     &
                               , power_A, power_B, alpha, beta, gama1, gama2  )
          c = c - 2.d0 * gama1 * gama2 * c1
        enddo
      enddo

      integ_grd += ao_coef_normalized_ordered_transp(l,j) &
                 * ao_coef_normalized_ordered_transp(m,i) * c
    enddo
  enddo

  print *, i, j, integ_grd

  return
end subroutine j1b_grd






!_____________________________________________________________________________________________________________
!
!
double precision function int_gauss_4G_deb( A_center, B_center, C_center1, C_center2, power_A, power_B &
                                          , alpha, beta, gama1, gama2 )

  include 'constants.include.F'

  implicit none

  integer         , intent(in) :: power_A(3), power_B(3)
  double precision, intent(in) :: A_center(3), B_center(3), C_center1(3), C_center2(3)
  double precision, intent(in) :: alpha, beta, gama1, gama2

  integer                      :: i, j, dim1, power_C
  integer                      :: iorder_AB(3), iorder_PQ(3), power_P(3), power_Q(3)
  double precision             :: AB_expo, fact_AB, AB_center(3), P_AB(0:max_dim,3)
  double precision             :: PQ_expo, fact_PQ, PQ_center(3), P_PQ(0:max_dim,3)
  double precision             :: gama, fact_C, C_center(3)
  double precision             :: cx0, cy0, cz0, c_tmp1, c_tmp2, cx, cy, cz
  double precision             :: int_tmp

  double precision             :: overlap_gaussian_x

  dim1 = 100

  call give_explicit_poly_and_gaussian( P_AB, AB_center, AB_expo, fact_AB &
                                      , iorder_AB, alpha, beta, power_A, power_B, A_center, B_center, dim1)

  int_tmp = 0.d0

  power_P = (/ 1, 0, 0 /)
  power_Q = (/ 1, 0, 0 /) 
  call give_explicit_poly_and_gaussian( P_PQ, PQ_center, PQ_expo, fact_PQ &
                                      , iorder_PQ, gama1, gama2, power_P, power_Q, C_center1, C_center2, dim1)
  cx = 0.d0
  do i = 0, iorder_AB(1)
    do j = 0, iorder_PQ(1)
      cx += P_AB(i,1) * P_PQ(j,1) * overlap_gaussian_x( AB_center(1), PQ_center(1), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  cy = 0.d0
  do i = 0, iorder_AB(2)
    do j = 0, iorder_PQ(2)
      cy += P_AB(i,2) * P_PQ(j,2) * overlap_gaussian_x( AB_center(2), PQ_center(2), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  cz = 0.d0
  do i = 0, iorder_AB(3)
    do j = 0, iorder_PQ(3)
      cz += P_AB(i,3) * P_PQ(j,3) * overlap_gaussian_x( AB_center(3), PQ_center(3), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  int_tmp += fact_PQ * cx * cy * cz  

  power_P = (/ 0, 1, 0 /)
  power_Q = (/ 0, 1, 0 /) 
  call give_explicit_poly_and_gaussian( P_PQ, PQ_center, PQ_expo, fact_PQ &
                                      , iorder_PQ, gama1, gama2, power_P, power_Q, C_center1, C_center2, dim1)
  cx = 0.d0
  do i = 0, iorder_AB(1)
    do j = 0, iorder_PQ(1)
      cx += P_AB(i,1) * P_PQ(j,1) * overlap_gaussian_x( AB_center(1), PQ_center(1), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  cy = 0.d0
  do i = 0, iorder_AB(2)
    do j = 0, iorder_PQ(2)
      cy += P_AB(i,2) * P_PQ(j,2) * overlap_gaussian_x( AB_center(2), PQ_center(2), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  cz = 0.d0
  do i = 0, iorder_AB(3)
    do j = 0, iorder_PQ(3)
      cz += P_AB(i,3) * P_PQ(j,3) * overlap_gaussian_x( AB_center(3), PQ_center(3), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  int_tmp += fact_PQ * cx * cy * cz  

  power_P = (/ 0, 0, 1 /)
  power_Q = (/ 0, 0, 1 /) 
  call give_explicit_poly_and_gaussian( P_PQ, PQ_center, PQ_expo, fact_PQ &
                                      , iorder_PQ, gama1, gama2, power_P, power_Q, C_center1, C_center2, dim1)
  cx = 0.d0
  do i = 0, iorder_AB(1)
    do j = 0, iorder_PQ(1)
      cx += P_AB(i,1) * P_PQ(j,1) * overlap_gaussian_x( AB_center(1), PQ_center(1), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  cy = 0.d0
  do i = 0, iorder_AB(2)
    do j = 0, iorder_PQ(2)
      cy += P_AB(i,2) * P_PQ(j,2) * overlap_gaussian_x( AB_center(2), PQ_center(2), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  cz = 0.d0
  do i = 0, iorder_AB(3)
    do j = 0, iorder_PQ(3)
      cz += P_AB(i,3) * P_PQ(j,3) * overlap_gaussian_x( AB_center(3), PQ_center(3), AB_expo, PQ_expo, i, j, dim1)
    enddo
  enddo
  int_tmp += fact_PQ * cx * cy * cz  

  int_gauss_4G_deb = fact_AB * int_tmp 

  return
end function int_gauss_4G_deb
!_____________________________________________________________________________________________________________
!_____________________________________________________________________________________________________________
