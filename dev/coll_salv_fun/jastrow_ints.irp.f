subroutine give_jastrow2_ovlp_ints_ao(mu,r1,n_taylor,ao_ints)
 implicit none
 BEGIN_DOC
! ao_ints(i,j) = \int dr2 \phi_i(r2) \phi_j(r2) ( J(r1,r2;\mu (r1)) )^2
!
! where \phi_i and \phi_j are AOs
 END_DOC
 double precision, intent(in) :: mu,r1(3) ! mu(r1) and r1
 integer, intent(in) :: n_taylor ! Order of the Taylor expansion of exp(x)
 double precision, intent(out):: ao_ints(ao_num,ao_num) ! integrals
 ao_ints = 0.d0
 integer :: i,j,k,l,num_A,num_B,power_A(3), power_B(3)
 double precision :: alpha, beta, A_center(3), B_center(3),c
 double precision :: ovlp_exp_f_phi_ij
 !$OMP PARALLEL                                                   &
    !$OMP DEFAULT (NONE)                                         &
    !$OMP PRIVATE (i,j,k,l,alpha,beta,A_center,B_center,power_A,power_B,&
    !$OMP          num_A,num_B,c)                      &
    !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
    !$OMP         ao_ints,n_taylor,r1,mu)
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
          do k=1,ao_prim_num(i)
            beta = ao_expo_ordered_transp(k,i)
            c = ovlp_exp_f_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)

            ao_ints(i,j) =  ao_ints(i,j) + ao_coef_normalized_ordered_transp(l,j)             &
                                         * ao_coef_normalized_ordered_transp(k,i) * c
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
end

subroutine give_jastrow2_erf_ints_ao(mu,r1,n_taylor,ao_ints)
 implicit none
 BEGIN_DOC
! ao_ints(i,j) = \int dr2 \phi_i(r2) \phi_j(r2) ( J(r1,r2;\mu (r1)) )^2 erf(mu(r1) r12)/r12
!
! where \phi_i and \phi_j are AOs
 END_DOC
 double precision, intent(in) :: mu,r1(3) ! mu(r1) and r1
 integer, intent(in) :: n_taylor ! Order of the Taylor expansion of exp(x)
 double precision, intent(out):: ao_ints(ao_num,ao_num) ! integrals
 ao_ints = 0.d0
 integer :: i,j,k,l,num_A,num_B,power_A(3), power_B(3)
 double precision :: alpha, beta, A_center(3), B_center(3),c
 double precision :: erf_exp_f_phi_ij
 !$OMP PARALLEL                                                   &
    !$OMP DEFAULT (NONE)                                         &
    !$OMP PRIVATE (i,j,k,l,alpha,beta,A_center,B_center,power_A,power_B,&
    !$OMP          num_A,num_B,c)                      &
    !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
    !$OMP         ao_ints,n_taylor,r1,mu)
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
          do k=1,ao_prim_num(i)
            beta = ao_expo_ordered_transp(k,i)
            c = erf_exp_f_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)

            ao_ints(i,j) =  ao_ints(i,j) + ao_coef_normalized_ordered_transp(l,j)             &
                                         * ao_coef_normalized_ordered_transp(k,i) * c
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
end

subroutine give_jastrow2_ovlp_ints_mo(mu,r1,n_taylor,mo_ints)
 implicit none
 BEGIN_DOC
! ao_ints(i,j) = \int dr2 \phi_i(r2) \phi_j(r2) ( J(r1,r2;\mu (r1)) )^2
!
! where \phi_i and \phi_j are MOs
 END_DOC
 double precision, intent(in) :: mu,r1(3) ! mu(r1) and r1
 integer, intent(in) :: n_taylor ! Order of the Taylor expansion of exp(x)
 double precision, intent(out):: mo_ints(mo_num,mo_num) ! integrals
 double precision :: ao_ints(ao_num,ao_num)
 mo_ints = 0.d0
 call give_jastrow2_ovlp_ints_ao(mu,r1,n_taylor,ao_ints)
 call ao_to_mo(ao_ints,ao_num,mo_ints,mo_num)

end

subroutine give_jastrow2_erf_ints_mo(mu,r1,n_taylor,mo_ints)
 implicit none
 BEGIN_DOC
! ao_ints(i,j) = \int dr2 \phi_i(r2) \phi_j(r2) ( J(r1,r2;\mu (r1)) )^2 erf(mu(r1) r12)/r12
!
! where \phi_i and \phi_j are MOs
 END_DOC
 double precision, intent(in) :: mu,r1(3) ! mu(r1) and r1
 integer, intent(in) :: n_taylor ! Order of the Taylor expansion of exp(x)
 double precision, intent(out):: mo_ints(mo_num,mo_num) ! integrals
 double precision :: ao_ints(ao_num,ao_num)
 mo_ints = 0.d0
 call give_jastrow2_erf_ints_ao(mu,r1,n_taylor,ao_ints)
 call ao_to_mo(ao_ints,ao_num,mo_ints,mo_num)

end


BEGIN_PROVIDER [double precision, ovlp_jastrow2_ij_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),mu,mo_ints(mo_num,mo_num)
 n_taylor = 4
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  mu = mu_of_r_prov(ipoint,1)
  call give_jastrow2_ovlp_ints_mo(mu,r1,n_taylor,mo_ints)
  do i = 1, mo_num
   do j = 1, mo_num
    ovlp_jastrow2_ij_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, erf_jastrow2_ij_mo, (mo_num,mo_num,n_points_final_grid)]
 implicit none
 integer :: i,j,ipoint,n_taylor
 double precision :: r1(3),mu,mo_ints(mo_num,mo_num)
 n_taylor = 4
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  mu = mu_of_r_prov(ipoint,1)
  call give_jastrow2_erf_ints_mo(mu,r1,n_taylor,mo_ints)
  do i = 1, mo_num
   do j = 1, mo_num
    erf_jastrow2_ij_mo(j,i,ipoint) = mo_ints(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 

