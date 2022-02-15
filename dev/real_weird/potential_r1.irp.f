 BEGIN_PROVIDER [ double precision, inv_dx_deriv]
&BEGIN_PROVIDER [ double precision, inv_dx_deriv_2]
 implicit none
 inv_dx_deriv = 1.d0/dx_deriv
 inv_dx_deriv_2 = inv_dx_deriv * inv_dx_deriv
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, potential_r1, (n_points_final_grid) ]
&BEGIN_PROVIDER [ double precision, kin_r1_plus,  (3,n_points_final_grid) ]
&BEGIN_PROVIDER [ double precision, kin_r1_minus, (3,n_points_final_grid) ]
 implicit none
 integer :: i,j,k
 double precision :: r(3),e0,r_tmp(3),u_dot_v,coef_0,sign_0
 double precision, allocatable :: v0(:),v0_plus(:), v0_minus(:)
 allocate(v0(mo_num), v0_plus(mo_num), v0_minus(mo_num))
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call diag_eff_h_at_r(r,e0,v0)
  potential_r1(i) = e0
  do k = 1, 3
   r_tmp = r
   r_tmp(k) += dx_deriv
   call diag_eff_h_at_r(r_tmp,e0,v0_plus)
   kin_r1_plus(k,i) = dabs(u_dot_v(v0,v0_plus,mo_num) )

   r_tmp = r
   r_tmp(k) -= dx_deriv
   call diag_eff_h_at_r(r_tmp,e0,v0_minus)
   kin_r1_minus(k,i) = dabs(u_dot_v(v0,v0_minus,mo_num) )
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, deriv_scal_prod, (3,n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, lapl_scal_prod, (3,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! deriv_scal_prod(m,i) = d/dx' <Psi_R|Psi_R'> |R=R'
!
! lapl_scal_prod(m,i) = d^2/dx'^2 <Psi_R|Psi_R'> |R=R'
 END_DOC
 integer :: i,j,k
 do i = 1, n_points_final_grid
  do k = 1, 3
   deriv_scal_prod(k,i) = 0.5d0 * ( kin_r1_plus(k,i) - kin_r1_minus(k,i)) * inv_dx_deriv
   lapl_scal_prod(k,i)  = (kin_r1_plus(k,i) + kin_r1_minus(k,i) - 2.d0) * inv_dx_deriv_2
  enddo
 enddo
END_PROVIDER 

subroutine diag_eff_h_at_r(r,e0,v0)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out):: e0,v0(mo_num)

 double precision mu 
 double precision, allocatable ::  integrals_ao(:,:), integrals_mo(:,:), h_r1(:,:)
 double precision, allocatable :: eigvalues(:),eigvectors(:,:)

 allocate( integrals_ao(ao_num,ao_num),integrals_mo(mo_num,mo_num),h_r1(mo_num, mo_num))
 allocate(eigvalues(mo_num),eigvectors(mo_num,mo_num))
 mu = 1.d+9
 call give_all_erf_kl_ao(integrals_ao,mu,r)
 call ao_to_mo(integrals_ao,size(integrals_ao,1),&
               integrals_mo,size(integrals_mo,1))
 h_r1 = mo_one_e_integrals
 h_r1 += integrals_mo 
 call lapack_diagd(eigvalues,eigvectors,h_r1,mo_num,mo_num)
 e0 = eigvalues(1)
 v0(:) = eigvectors(:,1)
end

 BEGIN_PROVIDER [ double precision, effective_H, (mo_num, mo_num)]
 implicit none
 effective_H = mo_integrals_n_e
 effective_H += kinetic_potential 
! effective_H += mo_kinetic_integrals
 effective_H += potential_from_r2 
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, e0_effective_H]
&BEGIN_PROVIDER [ double precision, effective_orbs, (mo_num,mo_num)]
 implicit none
 double precision, allocatable :: eigvalues(:),eigvectors(:,:)
 allocate(eigvalues(mo_num),eigvectors(mo_num,mo_num))
 print*,'HF approx = ',effective_H(1,1)
 call lapack_diagd(eigvalues,eigvectors,effective_H,mo_num,mo_num)
 e0_effective_H = eigvalues(1)
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, potential_from_r2,(mo_num, mo_num) ]
 implicit none
 integer :: i,j,ipoint
 double precision :: phi
 potential_from_r2 = 0.d0
 do ipoint = 1, n_points_final_grid
  do i = 1, mo_num
   do j = 1, mo_num
    potential_from_r2(j,i) += mos_in_r_array(i,ipoint) * mos_in_r_array(j,ipoint) & 
                            * potential_r1(ipoint) * final_weight_at_r_vector(ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, kinetic_potential, (mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l,ipoint
 kinetic_potential = 0.d0
 do ipoint = 1, n_points_final_grid
  do i = 1, mo_num
   do j = 1, mo_num
    do k = 1, 3
     kinetic_potential(j,i) += mos_in_r_array(i,ipoint) * final_weight_at_r_vector(ipoint) * (      & 
                             mo_num_lapl(k,j,ipoint)        * 1.d0                                  & 
                           + mos_in_r_array(j,ipoint)       * lapl_scal_prod(k,ipoint)              ) 
!                           + 2.d0 * deriv_scal_prod(k,ipoint) * mos_grad_in_r_array_tranp(k,j,ipoint) ) 
    enddo
   enddo
  enddo
 enddo
 kinetic_potential *= -0.5d0

 END_PROVIDER 
