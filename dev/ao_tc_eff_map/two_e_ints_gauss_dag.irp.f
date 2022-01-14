double precision function ao_tc_sym_two_e_pot_dag(i,j,k,l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) (tc_pot_dag(r12,mu)) k(r2) l(r2)
  !
  ! where (tc_pot_dag(r12,mu)) is the scalar part of the potential dagger excluding erf(mu r12)/r12
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer, intent(in)            :: i, j, k, l

  integer                        :: p, q, r, s
  integer                        :: num_i, num_j, num_k, num_l, dim1, I_power(3), J_power(3), K_power(3), L_power(3)
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: integral, schwartz_ij, thr
  double precision               :: p_inv, q_inv
  double precision               :: coef1, coef2, coef3, coef4 
  double precision               :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision               :: P_new(0:max_dim,3), P_center(3), fact_p, pp
  double precision               :: Q_new(0:max_dim,3), Q_center(3), fact_q, qq
  double precision, allocatable  :: schwartz_kl(:,:)

  double precision               :: scw_gauss_int, general_primitive_integral_gauss_dag

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)

  thr = ao_integrals_threshold*ao_integrals_threshold

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

  allocate( schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)) )
  schwartz_kl(0,0) = 0.d0

  ao_tc_sym_two_e_pot_dag = 0.d0
  do r = 1, ao_prim_num(k)
    coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
    schwartz_kl(0,r) = 0.d0
    do s = 1, ao_prim_num(l)
      coef2 = coef1 * ao_coef_normalized_ordered_transp(s,l) * ao_coef_normalized_ordered_transp(s,l)

      call give_explicit_poly_and_gaussian(Q_new, Q_center, qq, fact_q, iorder_q, &
            ao_expo_ordered_transp(r,k), ao_expo_ordered_transp(s,l),             &
            K_power, L_power, K_center, L_center, dim1)
      q_inv = 1.d0/qq

      scw_gauss_int = general_primitive_integral_gauss_dag(dim1, &
                Q_new, Q_center, fact_q, qq, q_inv, iorder_q,    &
                Q_new, Q_center, fact_q, qq, q_inv, iorder_q)

      schwartz_kl(s,r) = dabs(scw_gauss_int * coef2)
      schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
    enddo
    schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
  enddo

  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p,i)
    do q = 1, ao_prim_num(j)
      coef2 = coef1 * ao_coef_normalized_ordered_transp(q,j)

      call give_explicit_poly_and_gaussian(P_new, P_center, pp, fact_p, iorder_p, &
            ao_expo_ordered_transp(p,i), ao_expo_ordered_transp(q,j),             &
            I_power, J_power, I_center, J_center, dim1)
      p_inv = 1.d0/pp

      scw_gauss_int = general_primitive_integral_gauss_dag(dim1, &
                P_new, P_center, fact_p, pp, p_inv, iorder_p,    &
                P_new, P_center, fact_p, pp, p_inv, iorder_p)

      schwartz_ij = dabs(scw_gauss_int * coef2*coef2)
      if( schwartz_kl(0,0)*schwartz_ij < thr ) cycle

      do r = 1, ao_prim_num(k)
        if( schwartz_kl(0,r)*schwartz_ij < thr ) cycle
        coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
        do s = 1, ao_prim_num(l)
          if( schwartz_kl(s,r)*schwartz_ij < thr ) cycle
          coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)

          call give_explicit_poly_and_gaussian(Q_new, Q_center, qq, fact_q, iorder_q, &
                ao_expo_ordered_transp(r,k), ao_expo_ordered_transp(s,l),             &
                K_power, L_power, K_center, L_center, dim1)
          q_inv = 1.d0/qq

          integral = general_primitive_integral_gauss_dag(dim1, &
                P_new, P_center, fact_p, pp, p_inv, iorder_p,   &
                Q_new, Q_center, fact_q, qq, q_inv, iorder_q)

          ao_tc_sym_two_e_pot_dag = ao_tc_sym_two_e_pot_dag + coef4 * integral
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

  deallocate( schwartz_kl )

end function ao_tc_sym_two_e_pot_dag



double precision function general_primitive_integral_gauss_dag(dim, &
      P_new, P_center, fact_p, p, p_inv, iorder_p,                  &
      Q_new, Q_center, fact_q, q, q_inv, iorder_q)

  BEGIN_DOC
  ! Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives
  END_DOC

  include 'utils/constants.include.F'

  implicit none

  integer,          intent(in)   :: dim
  integer,          intent(in)   :: iorder_p(3)
  integer,          intent(in)   :: iorder_q(3)
  double precision, intent(in)   :: P_new(0:max_dim,3), P_center(3), fact_p, p, p_inv
  double precision, intent(in)   :: Q_new(0:max_dim,3), Q_center(3), fact_q, q, q_inv

  integer                        :: n_Ix, n_Iy, n_Iz, nx, ny, nz
  integer                        :: ix, iy, iz, jx, jy, jz, i
  integer                        :: n_pt_tmp, n_pt_out, iorder, m
  double precision               :: aa, c_a, t_a, rho_old, w_a, pi_3, prefactor, inv_pq_3_2
  double precision               :: gauss_int
  double precision               :: bla, thr, r_cut, gama_r_cut, rho, dist, rint
  double precision               :: a, b, c, d, e, f, accu, pq, const
  double precision               :: pq_inv, p10_1, p10_2, p01_1, p01_2, pq_inv_2
  double precision               :: Ix_pol(0:max_dim), Iy_pol(0:max_dim), Iz_pol(0:max_dim)
  double precision               :: dx(0:max_dim), dy(0:max_dim), dz(0:max_dim)
  double precision               :: d1(0:max_dim), d_poly(0:max_dim), d1_screened(0:max_dim)

  thr = ao_integrals_threshold

  general_primitive_integral_gauss_dag = 0.d0

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx, Ix_pol, dy, Iy_pol, dz, Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  ! Gaussian Product
  ! ----------------

  pq = p_inv*0.5d0*q_inv
  pq_inv = 0.5d0/(p+q)
  p10_1 = q*pq  ! 1/(2p)
  p01_1 = p*pq  ! 1/(2q)
  pq_inv_2 = pq_inv+pq_inv
  p10_2 = pq_inv_2 * p10_1*q !0.5d0*q/(pq + p*p)
  p01_2 = pq_inv_2 * p01_1*p !0.5d0*p/(q*q + pq)


  accu = 0.d0
  iorder = iorder_p(1)+iorder_q(1)+iorder_p(1)+iorder_q(1)
  do ix = 0, iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)
    if (abs(P_new(ix,1)) < thr) cycle
    a = P_new(ix,1)
    do jx = 0, iorder_q(1)
      d = a*Q_new(jx,1)
      if (abs(d) < thr) cycle
      !DIR$ FORCEINLINE
      call give_polynom_mult_center_x(P_center(1),Q_center(1),ix,jx,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dx,nx)
      !DIR$ FORCEINLINE
      call add_poly_multiply(dx,nx,d,Ix_pol,n_Ix)
    enddo
  enddo
  if (n_Ix == -1) then
    return
  endif
  iorder = iorder_p(2)+iorder_q(2)+iorder_p(2)+iorder_q(2)
  do ix=0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)
    if (abs(P_new(iy,2)) > thr) then
      b = P_new(iy,2)
      do jy = 0, iorder_q(2)
        e = b*Q_new(jy,2)
        if (abs(e) < thr) cycle
        !DIR$ FORCEINLINE
        call give_polynom_mult_center_x(P_center(2),Q_center(2),iy,jy,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dy,ny)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dy,ny,e,Iy_pol,n_Iy)
      enddo
    endif
  enddo
  if (n_Iy == -1) then
    return
  endif

  iorder = iorder_p(3)+iorder_q(3)+iorder_p(3)+iorder_q(3)
  do ix = 0, iorder
    Iz_pol(ix) = 0.d0
  enddo
  n_Iz = 0
  do iz = 0, iorder_p(3)
    if (abs(P_new(iz,3)) > thr) then
      c = P_new(iz,3)
      do jz = 0, iorder_q(3)
        f = c*Q_new(jz,3)
        if (abs(f) < thr) cycle
        !DIR$ FORCEINLINE
        call give_polynom_mult_center_x(P_center(3),Q_center(3),iz,jz,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dz,nz)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dz,nz,f,Iz_pol,n_Iz)
      enddo
    endif
  enddo
  if (n_Iz == -1) then
    return
  endif

  rho  = p * q * pq_inv_2
  dist = ( P_center(1) - Q_center(1) ) * ( P_center(1) - Q_center(1) ) + &
         ( P_center(2) - Q_center(2) ) * ( P_center(2) - Q_center(2) ) + &
         ( P_center(3) - Q_center(3) ) * ( P_center(3) - Q_center(3) )
  const = dist * rho

  n_pt_tmp = n_Ix + n_Iy
  do i = 0, n_pt_tmp
    d_poly(i) = 0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(Ix_pol, n_Ix, Iy_pol, n_Iy, d_poly, n_pt_tmp)
  if (n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp + n_Iz
  do i = 0, n_pt_out
    d1(i) = 0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(d_poly, n_pt_tmp, Iz_pol, n_Iz, d1, n_pt_out)

  pi_3       = pi * pi * pi
  inv_pq_3_2 = (p_inv * q_inv)**(1.5d0)
  rho_old    = (p*q)/(p+q)
  prefactor  = pi_3 * inv_pq_3_2 * fact_p * fact_q 
  gauss_int  = 0.d0
  do i = 1, n_gauss_eff_pot ! browse the gaussians with different expo/coef
 
    ! ---   ---   ---   ---   ---   ---
    aa  = expo_gauss_eff_pot_dag(i) 
    c_a = coef_gauss_eff_pot_dag(i)
    ! ---   ---   ---   ---   ---   ---

    t_a = dsqrt( aa /(rho_old + aa) ) 
    w_a = dexp(-t_a*t_a*rho_old*dist)

    accu = 0.d0
    ! evaluation of the polynom Ix(t_a) * Iy(t_a) * Iz(t_a)
    do m = 0, n_pt_out,2
      accu += d1(m) * (t_a)**(dble(m)) 
    enddo
    ! equation A8 of PRA-70-062505 (2004) of Toul. Col. Sav. 
    gauss_int = gauss_int + c_a * prefactor * (1.d0 - t_a*t_a)**(1.5d0) * w_a * accu
  enddo

  general_primitive_integral_gauss_dag = gauss_int

end function general_primitive_integral_gauss_dag

