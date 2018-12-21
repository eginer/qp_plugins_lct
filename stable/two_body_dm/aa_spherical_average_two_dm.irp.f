double precision function HF_two_body_dm_aa(r1,r2)
 implicit none
 BEGIN_DOC
! HF_two_body_dm_aa(r1,r2) = two body density for alpha/alpha spins in real space of the HF wave function (no 1/2 factor)
! < HF | 1/r12 | HF > = 1/2 \int(r1,r2) HF_two_body_dm_aa(r1,r2) 
 END_DOC
 double precision, intent(in) :: r1(2), r2(3)
 integer :: i,j
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 allocate(mos_array_r2(mo_tot_num), mos_array_r1(mo_tot_num))
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 HF_two_body_dm_aa = 0.d0
 do i =1, elec_alpha_num
  do j = 1, elec_alpha_num
   HF_two_body_dm_aa += mos_array_r1(j)**2 * mos_array_r2(i)**2 & ! direct term  
                      - mos_array_r1(i) * mos_array_r1(j) * mos_array_r2(i) * mos_array_r2(j)  
  enddo
 enddo
end

subroutine HF_two_body_dm_aa_spherical_laplacian(r1,r12,HF_two_bod,laplacian)
 implicit none
 double precision, intent(in)  :: r1(3),r12
 double precision, intent(out) :: HF_two_bod,laplacian
 BEGIN_DOC
! laplacian = spherical average of the second order derivative around r1 and evaluated at r12=0 of the two body density for alpha/alpha spins of HF wave function
! HF_two_bod = 1/2 * laplacian * r12**2
 END_DOC
 integer :: i,j
 double precision, allocatable :: mos_array_r1(:)
 double precision, allocatable :: mos_grad_array_r1(:,:)
 double precision, allocatable :: mos_lapl_array_r1(:,:)
 double precision :: lapl_j,lapl_i,grad_i,grad_j,grad_ij,direct_term,exchange_term
 allocate(mos_array_r1(mo_tot_num),mos_grad_array_r1(mo_tot_num,3),mos_lapl_array_r1(mo_tot_num,3))
 call give_all_mos_and_grad_and_lapl_at_r(r1,mos_array_r1,mos_grad_array_r1,mos_lapl_array_r1)
 direct_term  = 0.d0
 exchange_term = 0.d0
 ! direct term 
 do i = 1, elec_alpha_num
  do j = 1, elec_alpha_num
   lapl_j  = mos_lapl_array_r1(j,1)    + mos_lapl_array_r1(j,2)    + mos_lapl_array_r1(j,3)
   grad_j  = mos_grad_array_r1(j,1)**2 + mos_grad_array_r1(j,2)**2 + mos_grad_array_r1(j,3)**2
   lapl_i  = mos_lapl_array_r1(i,1)    + mos_lapl_array_r1(i,2)    + mos_lapl_array_r1(i,3)
   grad_i  = mos_grad_array_r1(i,1)**2 + mos_grad_array_r1(i,2)**2 + mos_grad_array_r1(i,3)**2
   grad_ij = mos_grad_array_r1(j,1) * mos_grad_array_r1(i,1) + mos_grad_array_r1(j,2) * mos_grad_array_r1(i,2) + mos_grad_array_r1(j,3) * mos_grad_array_r1(i,3) 
   direct_term   += 2.d0 * mos_array_r1(i) * (mos_array_r1(j) * lapl_j + 2.d0 * grad_j)
   exchange_term += mos_array_r1(i) * mos_array_r1(j) * ( 2.d0 * grad_ij + lapl_j * mos_array_r1(i) + lapl_i * mos_array_r1(j) )
  enddo
 enddo
 laplacian  =  0.33333333333333D0 * ( direct_term - exchange_term )
 HF_two_bod =  0.5d0 * laplacian * r12*r12
end
