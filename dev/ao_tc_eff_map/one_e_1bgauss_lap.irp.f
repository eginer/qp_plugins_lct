BEGIN_PROVIDER [ double precision, j1b_gauss_hermI, (ao_num,ao_num)]

  BEGIN_DOC
  !
  ! Hermitian part of 1-body Jastrow factow in the |AO| basis set.
  !
  !  :math:`\langle \chi_A | -0.5 \Delta \tau_{1b} | \chi_B \rangle` 
  !
  END_DOC

  implicit none

  integer          :: num_A, num_B
  integer          :: power_A(3), power_B(3)
  integer          :: i, j, k, l, m
  double precision :: alpha, beta, gama
  double precision :: A_center(3), B_center(3), C_center(3)
  double precision :: c1, c2, c

  integer          :: dim1
  double precision :: overlap_y, d_a_2, overlap_z, overlap

  double precision :: int_gauss_r0, int_gauss_r2

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
  
  j1b_gauss_hermI(1:ao_num,1:ao_num) = 0.d0

 !$OMP PARALLEL                                                 &
 !$OMP DEFAULT (NONE)                                           &
 !$OMP PRIVATE (i, j, k, l, m, alpha, beta, gama,               &
 !$OMP          A_center, B_center, C_center, power_A, power_B, &
 !$OMP          num_A, num_B, c1, c2, c)                        &
 !$OMP SHARED (ao_num, ao_prim_num, ao_expo_ordered_transp,     & 
 !$OMP         ao_power, ao_nucl, nucl_coord,                   &
 !$OMP         ao_coef_normalized_ordered_transp,               &
 !$OMP         nucl_num, j1b_gauss_pen, j1b_gauss_hermI)

 !$OMP DO SCHEDULE (dynamic)

  do j = 1, ao_num

    num_A         = ao_nucl(j)
    power_A(1:3)  = ao_power(j,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)

    do i = 1, ao_num

      num_B         = ao_nucl(i)
      power_B(1:3)  = ao_power(i,1:3)
      B_center(1:3) = nucl_coord(num_B,1:3)

      do l = 1, ao_prim_num(j)
        alpha = ao_expo_ordered_transp(l,j)

        do m = 1, ao_prim_num(i)
          beta = ao_expo_ordered_transp(m,i)

          c = 0.d0
          do k = 1, nucl_num

            gama          = j1b_gauss_pen(k)
            C_center(1:3) = nucl_coord(k,1:3)

            ! < XA | exp[-gama r_C^2] | XB >
            c1 = int_gauss_r0( A_center, B_center, C_center        &
                             , power_A, power_B, alpha, beta, gama )

            ! < XA | r_A^2 exp[-gama r_C^2] | XB >
            c2 = int_gauss_r2( A_center, B_center, C_center        &
                             , power_A, power_B, alpha, beta, gama )

            c = c + 3.d0 * gama * c1 - 2.d0 * gama * gama * c2
          enddo

          j1b_gauss_hermI(i,j) = j1b_gauss_hermI(i,j)      & 
                  + ao_coef_normalized_ordered_transp(l,j) &
                  * ao_coef_normalized_ordered_transp(m,i) * c

        enddo
      enddo
    enddo
  enddo

 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER


!_____________________________________________________________________________________________________________
!
!                             < XA | exp[-gama r_C^2] | XB >
!
double precision function int_gauss_r0(A_center, B_center, C_center, power_A, power_B, alpha, beta, gama)

  ! for max_dim
  include 'constants.include.F'

  implicit none

  integer         , intent(in) :: power_A(3), power_B(3)
  double precision, intent(in) :: A_center(3), B_center(3), C_center(3)
  double precision, intent(in) :: alpha, beta, gama 

  integer                      :: i, power_C, dim1
  integer                      :: iorder(3)
  integer                      :: nmax
  double precision             :: AB_expo, fact_AB, AB_center(3), P_AB(0:max_dim,3)
  double precision             :: cx, cy, cz 
  double precision             :: fact_p, P_center(3), p, F_integral_tab(0:max_dim)

  double precision             :: overlap_gaussian_x 
  double precision             :: F_integral

  dim1 = 100

  ! P_AB(0:max_dim,3) polynomial
  ! AB_center(3)      new center
  ! AB_expo           new exponent
  ! fact_AB           constant factor
  ! iorder(3)         i_order(i) = order of the polynomials
  call give_explicit_poly_and_gaussian( P_AB, AB_center, AB_expo, fact_AB &
                                      , iorder, alpha, beta, power_A, power_B, A_center, B_center, dim1)

  if( fact_AB .lt. 1d-20 ) then
    int_gauss_r0 = 0.d0
    return
  endif

   ! method 1
!
!  power_C = 0
!
!  cx = 0.d0
!  do i = 0, iorder(1)
!    cx = cx + P_AB(i,1) * overlap_gaussian_x(AB_center(1), C_center(1), AB_expo, gama, i, power_C, dim1)
!  enddo
!  cy = 0.d0
!  do i = 0, iorder(2)
!    cy = cy + P_AB(i,2) * overlap_gaussian_x(AB_center(2), C_center(2), AB_expo, gama, i, power_C, dim1)
!  enddo
!  cz = 0.d0
!  do i = 0, iorder(3)
!    cz = cz + P_AB(i,3) * overlap_gaussian_x(AB_center(3), C_center(3), AB_expo, gama, i, power_C, dim1)
!  enddo
  ! !!!


  ! method 2
  call gaussian_product(AB_expo, AB_center, gama, C_center, fact_p, p, P_center)
  fact_AB = fact_AB * fact_p

  nmax = maxval(iorder)
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i,p)
  enddo

  cx = P_AB(0,1) * F_integral_tab(0)
  cy = P_AB(0,2) * F_integral_tab(0)
  cz = P_AB(0,3) * F_integral_tab(0)

  do i = 1, iorder(1)
    cx = cx + P_AB(i,1) * F_integral_tab(i)
  enddo
  do i = 1, iorder(2)
    cy = cy + P_AB(i,2) * F_integral_tab(i)
  enddo
  do i = 1, iorder(3)
    cz = cz + P_AB(i,3) * F_integral_tab(i)
  enddo

  ! !!!

  int_gauss_r0 = fact_AB * cx * cy * cz

  return
end function int_gauss_r0 
!_____________________________________________________________________________________________________________
!_____________________________________________________________________________________________________________



!_____________________________________________________________________________________________________________
!
!                             < XA | r_C^2 exp[-gama r_C^2] | XB >
!
double precision function int_gauss_r2(A_center, B_center, C_center, power_A, power_B, alpha, beta, gama)

  ! for max_dim
  include 'constants.include.F'

  implicit none

  integer,          intent(in) :: power_A(3), power_B(3)
  double precision, intent(in) :: A_center(3), B_center(3), C_center(3)
  double precision, intent(in) :: alpha, beta, gama 

  integer                      :: i, power_C, dim1
  integer                      :: iorder(3)
  double precision             :: AB_expo, fact_AB, AB_center(3), P_AB(0:max_dim,3)
  double precision             :: cx0, cy0, cz0, cx, cy, cz

  double precision             :: overlap_gaussian_x

  dim1 = 100

  ! P_AB(0:max_dim,3) polynomial
  ! AB_center(3)      new center
  ! AB_expo           new exponent
  ! fact_AB           constant factor
  ! iorder(3)         i_order(i) = order of the polynomials
  call give_explicit_poly_and_gaussian( P_AB, AB_center, AB_expo, fact_AB &
                                      , iorder, alpha, beta, power_A, power_B, A_center, B_center, dim1)

  ! <<<
  ! to avoid multi-evaluation
  power_C = 0

  cx0 = 0.d0
  do i = 0, iorder(1)
    cx0 = cx0 + P_AB(i,1) * overlap_gaussian_x(AB_center(1), C_center(1), AB_expo, gama, i, power_C, dim1)
  enddo
  cy0 = 0.d0
  do i = 0, iorder(2)
    cy0 = cy0 + P_AB(i,2) * overlap_gaussian_x(AB_center(2), C_center(2), AB_expo, gama, i, power_C, dim1)
  enddo
  cz0 = 0.d0
  do i = 0, iorder(3)
    cz0 = cz0 + P_AB(i,3) * overlap_gaussian_x(AB_center(3), C_center(3), AB_expo, gama, i, power_C, dim1)
  enddo
  ! >>>

  int_gauss_r2 = 0.d0

  power_C = 2

  ! ( x - XC)^2
  cx = 0.d0
  do i = 0, iorder(1)
    cx = cx + P_AB(i,1) * overlap_gaussian_x(AB_center(1), C_center(1), AB_expo, gama, i, power_C, dim1)
  enddo
  int_gauss_r2 = int_gauss_r2 + fact_AB * cx * cy0 * cz0

  ! ( y - YC)^2
  cy = 0.d0
  do i = 0, iorder(2)
    cy = cy + P_AB(i,2) * overlap_gaussian_x(AB_center(2), C_center(2), AB_expo, gama, i, power_C, dim1)
  enddo
  int_gauss_r2 = int_gauss_r2 + fact_AB * cx0 * cy * cz0

  ! ( z - ZC)^2
  cz = 0.d0
  do i = 0, iorder(3)
    cz = cz + P_AB(i,3) * overlap_gaussian_x(AB_center(3), C_center(3), AB_expo, gama, i, power_C, dim1)
  enddo
  int_gauss_r2 = int_gauss_r2 + fact_AB * cx0 * cy0 * cz

  return
end function int_gauss_r2
!_____________________________________________________________________________________________________________
!_____________________________________________________________________________________________________________
