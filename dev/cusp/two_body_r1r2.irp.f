
 subroutine two_body_r1r2(r1,r2,istate,two_rdm)
 implicit none
 BEGIN_DOC
 ! gives the purely active on-top pair density for a given state 
 END_DOC
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out) :: two_rdm
 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)
 integer :: i,j,k,l
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)

 two_rdm = 0.d0
  do l = 1, mo_num ! 2
   do k = 1, mo_num ! 1
    do j = 1, mo_num ! 2
     do i = 1, mo_num ! 1
      !                               1 2 1 2
      two_rdm += full_occ_2_rdm_ab_mo(i,j,k,l,istate) * mos_array_r1(i) * mos_array_r2(j) * mos_array_r1(k) * mos_array_r2(l)
     enddo
    enddo
   enddo
  enddo
 end

 subroutine two_body_dm_hf_r1r2(r1,r2,two_rdm)
 implicit none
 BEGIN_DOC
 ! gives the purely active on-top pair density for a given state 
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out) :: two_rdm
 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)
 integer :: i,j
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)

 two_rdm = 0.d0
  do i = 1, elec_alpha_num
   do j = 1, elec_beta_num
    two_rdm +=  mos_array_r1(i)**2.d0 * mos_array_r2(j)**2.d0
   enddo
  enddo
 end
