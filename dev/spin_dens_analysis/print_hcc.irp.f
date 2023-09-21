program print_hcc
 implicit none
 read_wf = .True.
 touch read_wf
 print*,'Isotropic Hyperfine coupling constants in Gauss'
 call test
end
subroutine test
 implicit none
 double precision :: accu
 integer :: i,j
 print*,'Z               AU           GAUSS              MHZ             cm^-1'
 do i = 1, nucl_num
  write(*,'(I2,X,F3.1,X,4(F16.6,X))')i,nucl_charge(i),spin_density_at_nucleus(i),iso_hcc_gauss(i),iso_hcc_mhz(i),iso_hcc_cm_1(i)
 enddo
! print*,'Check MO/AO calculation'
! do i = 1, nucl_num
!  write(*,'(1000(F16.6,X))')spin_density_at_nucleus(i), spin_density_at_nucleus_from_mo(i),spin_density_at_nucleus_contrib_mo_test(i)
! enddo
! accu = 0.d0
! print*,'contribution per MOs '
! do i = 1, mo_tot_num
!  accu += spin_density_at_nucleus_contrib_per_mo(atom_number_hcc,i)
!  write(*,'(I3,X,1000(F16.6,X))')i,spin_density_at_nucleus_contrib_per_mo(atom_number_hcc,i),accu
! do j = 1, mo_tot_num
!  if(dabs(spin_density_at_nucleus_contrib_mo(atom_number_hcc,i,j)).gt.1.d-3)then
!   accu += spin_density_at_nucleus_contrib_mo(atom_number_hcc,i,j)
!   if(j.ge.i)then
!    print*,'i,j = ',i,j
!    print*,spin_density_at_nucleus_contrib_mo(atom_number_hcc,i,j)
!   endif
!  endif
! enddo
! enddo
! print*,'accu = ',accu

end


