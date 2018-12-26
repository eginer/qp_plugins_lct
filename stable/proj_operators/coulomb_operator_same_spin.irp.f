double precision function f_HF_aa(r1,r2)
 implicit none
 BEGIN_DOC
! f_HF_aa(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha alpha spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > = 0.5 * \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 f_HF_aa = 0.d0
 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_aa += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r2(j) & 
                * (mos_array_r1(k) * mos_array_r2(l) - mos_array_r1(l) * mos_array_r2(k)) ! direct and exchange term
    enddo
   enddo
  enddo
 enddo
end


subroutine f_HF_aa_spherical_averaged(r1,r12,laplacian,value_f_HF_aa)
 implicit none
 BEGIN_DOC
 ! second order derivative of f_HF_aa(r1,r2) around r1 and spherically averaged 
 END_DOC
 double precision, intent(in) :: r1(3),r12
 double precision, intent(out):: laplacian,value_f_HF_aa
 integer :: i,j,k,l
 double precision, allocatable :: mos_grad_array_r1(:,:),mos_array_r1(:)
 double precision, allocatable :: mos_lapl_array_r1(:,:)
 double precision :: lapl_j,grad_j,lapl_k,grad_k,lapl_l,grad_l,grad_jl,grad_jk,direct_term,exchange_term
 allocate(mos_array_r1(mo_tot_num),mos_grad_array_r1(mo_tot_num,3),mos_lapl_array_r1(mo_tot_num,3))
 call give_all_mos_and_grad_and_lapl_at_r(r1,mos_array_r1,mos_grad_array_r1,mos_lapl_array_r1)
 laplacian = 0.d0
 value_f_HF_aa = 0.d0
 direct_term = 0.d0
 exchange_term = 0.d0
 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     lapl_j  = mos_lapl_array_r1(j,1)    + mos_lapl_array_r1(j,2)    + mos_lapl_array_r1(j,3)
     grad_j  = mos_grad_array_r1(j,1)**2 + mos_grad_array_r1(j,2)**2 + mos_grad_array_r1(j,3)**2
     lapl_k  = mos_lapl_array_r1(k,1)    + mos_lapl_array_r1(k,2)    + mos_lapl_array_r1(k,3)
     grad_k  = mos_grad_array_r1(k,1)**2 + mos_grad_array_r1(k,2)**2 + mos_grad_array_r1(k,3)**2
     lapl_l  = mos_lapl_array_r1(l,1)    + mos_lapl_array_r1(l,2)    + mos_lapl_array_r1(l,3)
     grad_l  = mos_grad_array_r1(l,1)**2 + mos_grad_array_r1(l,2)**2 + mos_grad_array_r1(l,3)**2
     grad_jl = mos_grad_array_r1(j,1) * mos_grad_array_r1(l,1) + mos_grad_array_r1(j,2) * mos_grad_array_r1(l,2) + mos_grad_array_r1(j,3) * mos_grad_array_r1(l,3) 
     grad_jk = mos_grad_array_r1(j,1) * mos_grad_array_r1(k,1) + mos_grad_array_r1(j,2) * mos_grad_array_r1(k,2) + mos_grad_array_r1(j,3) * mos_grad_array_r1(k,3) 
     direct_term   += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(k)   &
                   * ( 2.d0 * grad_jl + lapl_j * mos_array_r1(l) + lapl_l * mos_array_r1(j) )  ! direct term
                                                                                                                             
     exchange_term += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(l)   & 
                   * ( 2.d0 * grad_jk + lapl_j * mos_array_r1(k) + lapl_k * mos_array_r1(j) )  ! exchange term
                
    enddo
   enddo
  enddo
 enddo
 laplacian     =  0.33333333333333D0 * ( direct_term - exchange_term )
!laplacian     =  0.33333333333333D0 * ( direct_term )
 value_f_HF_aa =  0.5d0 * laplacian * r12*r12
end

double precision function f_HF_aa_integrated(r1)
 implicit none
 BEGIN_DOC
! f_HF_aa_integrated(X_1) = function int(X_2) f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha alpha spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > = 0.5 * \int (X1) f_HF_aa_integrated(X1)
 END_DOC
 double precision, intent(in) :: r1(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 f_HF_aa_integrated = 0.d0
 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = j, j ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_aa_integrated += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(k)  
    enddo
   enddo
  enddo
 enddo

 do j = 1, elec_alpha_num ! electron 2
  do i = 1, elec_alpha_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = j, j ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_aa_integrated -= integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(l)  
    enddo
   enddo
  enddo
 enddo

end


double precision function f_HF_bb(r1,r2)
 implicit none
 BEGIN_DOC
! f_HF_bb(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for beta beta spins and an HF wave function
! < HF | wee_{\beta\beta} | HF > = 0.5 * \int (X1,X2) f_HF_bb(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(mo_tot_num)
 double precision :: mos_array_r2(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 f_HF_bb = 0.d0
 do j = 1, elec_beta_num ! electron 2
  do i = 1, elec_beta_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_bb += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r2(j) & 
                * (mos_array_r1(k) * mos_array_r2(l) - mos_array_r1(l) * mos_array_r2(k)) ! direct and exchange term
    enddo
   enddo
  enddo
 enddo
end

double precision function f_HF_bb_integrated(r1)
 implicit none
 BEGIN_DOC
! f_HF_bb_integrated(X_1) = function int(X_2) f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for beta beta spins and an HF wave function
! < HF | wee_{\beta\beta} | HF > = 0.5 * \int (X1) f_HF_bb_integrated(X1)
 END_DOC
 double precision, intent(in) :: r1(3)
 integer :: i,j,k,l
 double precision :: mos_array_r1(mo_tot_num)
 call give_all_mos_at_r(r1,mos_array_r1) 
 f_HF_bb_integrated = 0.d0
 do j = 1, elec_beta_num ! electron 2
  do i = 1, elec_beta_num ! electron 1 
   do l = j, j ! electron 2 
    do k = 1, mo_tot_num ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_bb_integrated += integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(k)  
    enddo
   enddo
  enddo
 enddo

 do j = 1, elec_beta_num ! electron 2
  do i = 1, elec_beta_num ! electron 1 
   do l = 1, mo_tot_num ! electron 2 
    do k = j, j ! electron 1 
     !                                     1 2 1 2 : <ij|kl> 
     f_HF_bb_integrated -= integrals_for_hf_potential(k,l,i,j) *  mos_array_r1(i) * mos_array_r1(l)  
    enddo
   enddo
  enddo
 enddo

end


