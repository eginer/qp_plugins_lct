
BEGIN_PROVIDER [ double precision, ao_erf_integrals_n_e, (ao_num,ao_num)]
  BEGIN_DOC
  !  Nucleus-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  !
  !  These integrals also contain the pseudopotential integrals.
  END_DOC
  implicit none
  double precision               :: alpha, beta, gama, delta
  integer                        :: num_A,num_B
  double precision               :: A_center(3),B_center(3),C_center(3)
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt_in,m
  double precision               :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult_erf,mu_in

  mu_in = mu_one_e_j
   ao_erf_integrals_n_e = 0.d0

    !        _
    ! /|  / |_)
    !  | /  | \
    !

    !$OMP PARALLEL                                                   &
        !$OMP DEFAULT (NONE)                                         &
        !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,C_center,power_A,power_B,&
        !$OMP          num_A,num_B,Z,c,n_pt_in)                      &
        !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
        !$OMP         n_pt_max_integrals,ao_erf_integrals_n_e,nucl_num,nucl_charge,mu_in)

    n_pt_in = n_pt_max_integrals

    !$OMP DO SCHEDULE (dynamic)

    do j = 1, ao_num
      num_A = ao_nucl(j)
      power_A(1:3)= ao_power(j,1:3)
      A_center(1:3) = nucl_coord(num_A,1:3)

      do i = 1, ao_num

        num_B = ao_nucl(i)
        power_B(1:3)= ao_power(i,1:3)
        B_center(1:3) = nucl_coord(num_B,1:3)

        do l=1,ao_prim_num(j)
          alpha = ao_expo_ordered_transp(l,j)

          do m=1,ao_prim_num(i)
            beta = ao_expo_ordered_transp(m,i)

            double precision               :: c
            c = 0.d0

            do  k = 1, nucl_num
              double precision               :: Z
              Z = nucl_charge(k)

              C_center(1:3) = nucl_coord(k,1:3)

              c = c - Z * NAI_pol_mult_erf(A_center,B_center,        &
                  power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)

            enddo
            ao_erf_integrals_n_e(i,j) = ao_erf_integrals_n_e(i,j)  &
                + ao_coef_normalized_ordered_transp(l,j)             &
                * ao_coef_normalized_ordered_transp(m,i) * c
          enddo
        enddo
      enddo
    enddo

    !$OMP END DO
    !$OMP END PARALLEL
!    IF (DO_PSEUDO) THEN
!       ao_erf_integrals_n_e += ao_pseudo_integrals
!    ENDIF



END_PROVIDER

 BEGIN_PROVIDER [double precision, mo_erf_integrals_n_e, (mo_num, mo_num)]
 implicit none
    call ao_to_mo(ao_erf_integrals_n_e,size(ao_erf_integrals_n_e,1),mo_erf_integrals_n_e,size(mo_erf_integrals_n_e,1))
 END_PROVIDER 
