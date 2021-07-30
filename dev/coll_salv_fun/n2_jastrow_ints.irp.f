double precision function int_ovlp_n2_jaswtrow2(r1,muj,istate,n_taylor,exponent_exp)
 double precision, intent(in) :: r1(3),muj,exponent_exp
 integer, intent(in) :: istate,n_taylor
 integer :: i,j,k,l
 double precision :: phi_ik
 double precision, allocatable :: mo_ints(:,:),mos_array(:)
 allocate(mo_ints(mo_num,mo_num),mos_array(mo_num))
 provide full_occ_2_rdm_ab_chemist_mo 
 print*,'n_taylor = ',n_taylor
 print*,'exponent_exp = ',exponent_exp
 call give_jastrow2_ovlp_ints_mo(muj,r1,n_taylor,mo_ints,exponent_exp)
 call give_all_mos_at_r(r1,mos_array)
 int_ovlp_n2_jaswtrow2 = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   phi_ik = mos_array(i) * mos_array(k)
   do j = 1, mo_num
    do l = 1, mo_num
     int_ovlp_n2_jaswtrow2 += mo_ints(l,k) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate) * phi_ik
    enddo
   enddo
  enddo
 enddo

end

double precision function int_erf_n2_jaswtrow2(r1,muj,muc,istate,n_taylor)
 double precision, intent(in) :: r1(3),muj,muc
 integer, intent(in) :: istate,n_taylor
 integer :: i,j,k,l
 double precision :: phi_ik
 double precision, allocatable :: mo_ints(:,:),mos_array(:)
 allocate(mo_ints(mo_num,mo_num),mos_array(mo_num))
 provide full_occ_2_rdm_ab_chemist_mo 
 call give_jastrow2_erf_ints_mo(muj,muc,r1,n_taylor,mo_ints)
 call give_all_mos_at_r(r1,mos_array)
 int_erf_n2_jaswtrow2 = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   phi_ik = mos_array(i) * mos_array(k)
   do j = 1, mo_num
    do l = 1, mo_num
     int_erf_n2_jaswtrow2 += mo_ints(l,k) * full_occ_2_rdm_ab_chemist_mo(l,j,k,i,istate) * phi_ik
    enddo
   enddo
  enddo
 enddo

end
