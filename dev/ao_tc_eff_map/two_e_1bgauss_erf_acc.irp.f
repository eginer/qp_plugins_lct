double precision function j1b_gauss_erf_acc(i, j, k, l)

  BEGIN_DOC
  ! 
  !  integral in the AO basis:
  !     i(r1) j(r1) f(r12) k(r2) l(r2)
  !
  !  with:
  !     f(r12) = - [ -0.5 erf(mu r12) / r12 ] (r1-r2) \cdot \sum_A (-2 a_A) [ r1A exp(-aA r1A^2) - r2A exp(-aA r2A^2) ]
  !            = - [ erf(mu r12) / r12 ] \sum_A a_A [ (r1-RA)^2 exp(-aA r1A^2)
  !                                                 + (r2-RA)^2 exp(-aA r2A^2) 
  !                                                 - (r1-RA) \cdot (r2-RA) exp(-aA r1A^2)
  !                                                 - (r1-RA) \cdot (r2-RA) exp(-aA r2A^2) ]
  !
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer, intent(in) :: i, j, k, l

  integer             :: p, q, r, s, ii
  integer             :: num_i, num_j, num_k, num_l, num_ii 
  integer             :: I_power(3), J_power(3), K_power(3), L_power(3)
  integer             :: iorder_p(3), iorder_q(3)
  integer             :: shift_P(3), shift_Q(3)
  integer             :: dim1

  double precision    :: coef1, coef2, coef3, coef4
  double precision    :: expo1, expo2, expo3, expo4
  double precision    :: p1_inv, q1_inv, p2_inv, q2_inv
  double precision    :: P1_new(0:max_dim,3), P1_center(3), fact_p1, pp1
  double precision    :: P2_new(0:max_dim,3), P2_center(3), fact_p2, pp2
  double precision    :: Q1_new(0:max_dim,3), Q1_center(3), fact_q1, qq1
  double precision    :: Q2_new(0:max_dim,3), Q2_center(3), fact_q2, qq2
  double precision    :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision    :: expoii, factii, Centerii(3)
  double precision    :: ff, gg, cx, cy, cz

  double precision    :: general_primitive_integral_erf_shifted
  !double precision    :: j1b_gauss_erf_schwartz_accel
  
  PROVIDE j1b_gauss_pen

  dim1 = n_pt_max_integrals

  ! TODO
  !if( ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024 ) then
  !  j1b_gauss_erf_schwartz_accel = j1b_gauss_erf_schwartz_accel(i, j, k, l)
  !  return
  !endif

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)

  do p = 1, 3
    I_power(p)  = ao_power(i,p)
    J_power(p)  = ao_power(j,p)
    K_power(p)  = ao_power(k,p)
    L_power(p)  = ao_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = nucl_coord(num_j,p)
    K_center(p) = nucl_coord(num_k,p)
    L_center(p) = nucl_coord(num_l,p)
  enddo

  j1b_gauss_erf_acc = 0.d0

  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p, i)
    expo1 = ao_expo_ordered_transp(p, i)

    do q = 1, ao_prim_num(j)
      coef2 = coef1 * ao_coef_normalized_ordered_transp(q, j)
      expo2 = ao_expo_ordered_transp(q, j)

      call give_explicit_poly_and_gaussian( P1_new, P1_center, pp1, fact_p1, iorder_p, expo1, expo2 &
                                          , I_power, J_power, I_center, J_center, dim1 )
      p1_inv = 1.d0 / pp1

      do r = 1, ao_prim_num(k)
        coef3 = coef2 * ao_coef_normalized_ordered_transp(r, k)
        expo3 = ao_expo_ordered_transp(r, k)

        do s = 1, ao_prim_num(l)
          coef4 = coef3 * ao_coef_normalized_ordered_transp(s, l)
          expo4 = ao_expo_ordered_transp(s, l)
 
          call give_explicit_poly_and_gaussian( Q1_new, Q1_center, qq1, fact_q1, iorder_q, expo3, expo4 &
                                              , K_power, L_power, K_center, L_center, dim1 )
          q1_inv = 1.d0 / qq1

          cx = 0.d0
          cy = 0.d0
          cz = 0.d0
          do ii = 1, nucl_num
            expoii        = j1b_gauss_pen(ii)
            Centerii(1:3) = nucl_coord(ii, 1:3)

            call gaussian_product(pp1, P1_center, expoii, Centerii, factii, pp2, P2_center)
            fact_p2 = fact_p1 * factii
            p2_inv  = 1.d0 / pp2
            call pol_modif_center( P1_center, P2_center, iorder_p, P1_new, P2_new)

            call gaussian_product(qq1, Q1_center, expoii, Centerii, factii, qq2, Q2_center)
            fact_q2 = fact_q1 * factii
            q2_inv  = 1.d0 / qq2
            call pol_modif_center( Q1_center, Q2_center, iorder_q, Q1_new, Q2_new)


            ! ----------------------------------------------------------------------------------------------------
            !                         [ erf(mu r12) / r12 ] \sum_A a_A [ (r1-RA)^2 exp(-aA r1A^2)
            ! ----------------------------------------------------------------------------------------------------

            shift_Q = (/ 0, 0, 0 /)

            ! x term:
            ff = P2_center(1) - Centerii(1) 

            shift_P = (/ 2, 0, 0 /)
            cx = cx + expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_P = (/ 1, 0, 0 /)
            cx = cx + expoii * 2.d0 * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P       &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q       )

            shift_P = (/ 0, 0, 0 /)
            cx = cx + expoii * ff * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P     &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q     )

            ! y term:
            ff = P2_center(2) - Centerii(2) 

            shift_P = (/ 0, 2, 0 /)
            cy = cy + expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_P = (/ 0, 1, 0 /)
            cy = cy + expoii * 2.d0 * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P       &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q       )

            shift_P = (/ 0, 0, 0 /)
            cy = cy + expoii * ff * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P     &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q     )

            ! z term:
            ff = P2_center(3) - Centerii(3) 

            shift_P = (/ 0, 0, 2 /)
            cz = cz + expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_P = (/ 0, 0, 1 /)
            cz = cz + expoii * 2.d0 * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P       &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q       )

            shift_P = (/ 0, 0, 0 /)
            cz = cz + expoii * ff * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P     &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q     )

            ! ----------------------------------------------------------------------------------------------------



            ! ----------------------------------------------------------------------------------------------------
            !                         [ erf(mu r12) / r12 ] \sum_A a_A [ (r2-RA)^2 exp(-aA r2A^2)
            ! ----------------------------------------------------------------------------------------------------

            shift_P = (/ 0, 0, 0 /)

            ! x term:
            ff = Q2_center(1) - Centerii(1) 

            shift_Q = (/ 2, 0, 0 /)
            cx = cx + expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_Q = (/ 1, 0, 0 /)
            cx = cx + expoii * 2.d0 * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P       &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q       )

            shift_Q = (/ 0, 0, 0 /)
            cx = cx + expoii * ff * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P     &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q     )

            ! y term:
            ff = Q2_center(2) - Centerii(2) 

            shift_Q = (/ 0, 2, 0 /)
            cy = cy + expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_Q = (/ 0, 1, 0 /)
            cy = cy + expoii * 2.d0 * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P       &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q       )

            shift_Q = (/ 0, 0, 0 /)
            cy = cy + expoii * ff * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P     &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q     )

            ! z term:
            ff = Q2_center(3) - Centerii(3) 

            shift_Q = (/ 0, 0, 2 /)
            cz = cz + expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_Q = (/ 0, 0, 1 /)
            cz = cz + expoii * 2.d0 * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P       &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q       )

            shift_Q = (/ 0, 0, 0 /)
            cz = cz + expoii * ff * ff * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P     &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q     )

            ! ----------------------------------------------------------------------------------------------------



            ! ----------------------------------------------------------------------------------------------------
            !                 - [ erf(mu r12) / r12 ] \sum_A a_A [ (r1-RA) \cdot (r2-RA) exp(-aA r1A^2) ]
            ! ----------------------------------------------------------------------------------------------------

            ! x term:
            ff = P2_center(1) - Centerii(1) 
            gg = Q1_center(1) - Centerii(1) 

            shift_p = (/ 1, 0, 0 /)
            shift_Q = (/ 1, 0, 0 /)
            cx = cx - expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 1, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cx = cx - expoii * gg * general_primitive_integral_erf_shifted( dim1  &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 1, 0, 0 /)
            cx = cx - expoii * ff * general_primitive_integral_erf_shifted( dim1  &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cx = cx - expoii * ff * gg * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P     &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q     )

            ! y term:
            ff = P2_center(2) - Centerii(2) 
            gg = Q1_center(2) - Centerii(2) 

            shift_p = (/ 0, 1, 0 /)
            shift_Q = (/ 0, 1, 0 /)
            cy = cy - expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 1, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cy = cy - expoii * gg * general_primitive_integral_erf_shifted( dim1  &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 1, 0 /)
            cy = cy - expoii * ff * general_primitive_integral_erf_shifted( dim1  &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cy = cy - expoii * ff * gg * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P     &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q     )

            ! z term:
            ff = P2_center(3) - Centerii(3) 
            gg = Q1_center(3) - Centerii(3) 

            shift_p = (/ 0, 0, 1 /)
            shift_Q = (/ 0, 0, 1 /)
            cz = cz - expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 1 /)
            shift_Q = (/ 0, 0, 0 /)
            cz = cz - expoii * gg * general_primitive_integral_erf_shifted( dim1  &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 1 /)
            cz = cz - expoii * ff * general_primitive_integral_erf_shifted( dim1  &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cz = cz - expoii * ff * gg * general_primitive_integral_erf_shifted( dim1 &
                     , P2_new, P2_center, fact_p2, pp2, p2_inv, iorder_p, shift_P     &
                     , Q1_new, Q1_center, fact_q1, qq1, q1_inv, iorder_q, shift_Q     )

            ! ----------------------------------------------------------------------------------------------------



            ! ----------------------------------------------------------------------------------------------------
            !                 - [ erf(mu r12) / r12 ] \sum_A a_A [ (r1-RA) \cdot (r2-RA) exp(-aA r2A^2) ]
            ! ----------------------------------------------------------------------------------------------------

            ! x term:
            ff = P1_center(1) - Centerii(1) 
            gg = Q2_center(1) - Centerii(1) 

            shift_p = (/ 1, 0, 0 /)
            shift_Q = (/ 1, 0, 0 /)
            cx = cx - expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 1, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cx = cx - expoii * gg * general_primitive_integral_erf_shifted( dim1  &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 1, 0, 0 /)
            cx = cx - expoii * ff * general_primitive_integral_erf_shifted( dim1  &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cx = cx - expoii * ff * gg * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P     &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q     )

            ! y term:
            ff = P1_center(2) - Centerii(2) 
            gg = Q2_center(2) - Centerii(2) 

            shift_p = (/ 0, 1, 0 /)
            shift_Q = (/ 0, 1, 0 /)
            cy = cy - expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 1, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cy = cy - expoii * gg * general_primitive_integral_erf_shifted( dim1  &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 1, 0 /)
            cy = cy - expoii * ff * general_primitive_integral_erf_shifted( dim1  &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cy = cy - expoii * ff * gg * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P     &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q     )

            ! z term:
            ff = P1_center(3) - Centerii(3) 
            gg = Q2_center(3) - Centerii(3) 

            shift_p = (/ 0, 0, 1 /)
            shift_Q = (/ 0, 0, 1 /)
            cz = cz - expoii * general_primitive_integral_erf_shifted( dim1       &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 1 /)
            shift_Q = (/ 0, 0, 0 /)
            cz = cz - expoii * gg * general_primitive_integral_erf_shifted( dim1  &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 1 /)
            cz = cz - expoii * ff * general_primitive_integral_erf_shifted( dim1  &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q )

            shift_p = (/ 0, 0, 0 /)
            shift_Q = (/ 0, 0, 0 /)
            cz = cz - expoii * ff * gg * general_primitive_integral_erf_shifted( dim1 &
                     , P1_new, P1_center, fact_p1, pp1, p1_inv, iorder_p, shift_P     &
                     , Q2_new, Q2_center, fact_q2, qq2, q2_inv, iorder_q, shift_Q     )

            ! ----------------------------------------------------------------------------------------------------

          enddo

          j1b_gauss_erf_acc = j1b_gauss_erf_acc - coef4 * ( cx + cy + cz )
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

  return
end function j1b_gauss_erf_acc
