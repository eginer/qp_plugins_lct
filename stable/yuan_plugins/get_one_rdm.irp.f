BEGIN_PROVIDER [ double precision, one_rdm_from_two_rdm, (mo_num, mo_num)]
 implicit none
 integer :: i,j,k
 one_rdm_from_two_rdm = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    one_rdm_from_two_rdm(j,i) += 2.d0 * full_occ_2_rdm_spin_trace_mo(j,k,i,k,1) / (dble(elec_num)-1.d0)
   enddo
  enddo
 enddo
END_PROVIDER 

program print_one_rdm
 implicit none
 read_wf = .True.
 touch read_wf
 call pouet 
end

subroutine pouet
 implicit none
 integer :: i
 print*,'one_rdm_from_two_rdm'
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')one_rdm_from_two_rdm(i,:)
 enddo
 print*,'one_e_dm_mo'
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')one_e_dm_mo(i,:)
 enddo
end
