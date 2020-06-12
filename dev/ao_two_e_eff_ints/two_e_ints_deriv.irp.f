double precision function ao_two_e_integral_schwartz_accel_dr12_coul(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|1/r12 * r_{12}.d/dr12 jl> or (ij|kl)
  !  i(r1) k(r2) 1/r12 (x_1 - x_2) (d/dx1 - d/dx2) + (y_1 - y_2) (d/dy1 - d/dy2) + (z_1 - z_2) (d/dz1 - d/dz2) j(r1) l(r2)
  !  
  ! WARNING <ik|jl> IS NOT EQUAL TO <kl|ik> because of the first-order differential operator
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision, allocatable  :: schwartz_kl(:,:)
  double precision               :: schwartz_ij
  double precision :: scw_gauss_int,general_primitive_integral_gauss

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_two_e_integral_schwartz_accel_gauss = 0.d0
  double precision               :: thr
  thr = ao_integrals_threshold*ao_integrals_threshold

  allocate(schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)))

      double precision               :: coef3
      double precision               :: coef2
      double precision               :: p_inv,q_inv
      double precision               :: coef1
      double precision               :: coef4

    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(s,l) * ao_coef_normalized_ordered_transp(s,l)
        call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
            ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),                 &
            K_power,L_power,K_center,L_center,dim1)
        q_inv = 1.d0/qq
        scw_gauss_int = general_primitive_integral_gauss(dim1,              &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)

        schwartz_kl(s,r) = scw_gauss_int * coef2
        schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
            ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
            I_power,J_power,I_center,J_center,dim1)
        p_inv = 1.d0/pp
        scw_gauss_int = general_primitive_integral_gauss(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                P_new,P_center,fact_p,pp,p_inv,iorder_p)
        schwartz_ij = scw_gauss_int * coef2*coef2
        if (schwartz_kl(0,0)*schwartz_ij < thr) then
           cycle
        endif
        do r = 1, ao_prim_num(k)
          if (schwartz_kl(0,r)*schwartz_ij < thr) then
             cycle
          endif
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            if (schwartz_kl(s,r)*schwartz_ij < thr) then
               cycle
            endif
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q, &
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),            &
                K_power,L_power,K_center,L_center,dim1)
            q_inv = 1.d0/qq
            integral = general_primitive_integral_gauss(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
            ao_two_e_integral_schwartz_accel_gauss = ao_two_e_integral_schwartz_accel_gauss + coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  deallocate (schwartz_kl)

end

