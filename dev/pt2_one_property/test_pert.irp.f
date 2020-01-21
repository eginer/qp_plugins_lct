program test
 implicit none
  read_wf = .True.
  touch read_wf
  call routine
end

subroutine routine
 implicit none
 integer :: i,j,degree,nsingles,idx_hf
 integer, allocatable :: idx(:)
 double precision, allocatable :: coef(:)
 double precision :: hij,hjj,eref,oij
 double precision :: e2
 allocate(idx(N_det),coef(N_det))
 nsingles = 0
 eref = ref_bitmask_energy
 e2 = 0.d0
 do i = 1, N_det
  call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree == 0)then
   idx_hf = i
!   coef_hf = psi_coef(i,1)
   coef_hf = 1.d0
  else 
   nsingles += 1
   idx(nsingles) = i
    call i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hjj)
    call i_H_j(ref_bitmask,psi_det(1,1,i),N_int,hij)
    coef(nsingles) = hij/(eref - hjj)
    e2 += coef(nsingles) * hij 
!    coef(nsingles) = psi_coef(i,1)
  endif
 enddo
 double precision :: accu0j,accu00,accujj_diag,accujj_off_diag,normpsi1,z_dip,coef_hf,normpsi0
 double precision :: accu_diag,accu_off_diag

 ! part <HF|O|HF>
 call i_H_j_eff_pot_diag(ref_bitmask,mo_dipole_z,mo_dipole_z,mo_num,N_int,accu00)
 accu00 = - accu00
 normpsi0 = coef_hf * coef_hf

 accu_diag = accu00 * normpsi0

 ! part <HF|O|PSI_1>
 normpsi1 = 0.d0
 do i = 1, nsingles
  call i_H_j_eff_pot(ref_bitmask,psi_det(1,1,idx(i)),mo_dipole_z,mo_dipole_z,mo_num,N_int,oij)
  normpsi1 += coef(i) * coef(i)
  accu0j -= 2.d0 * oij * coef(i) * coef_hf
  accu_off_diag -= 2.d0 * oij * coef(i) * coef_hf
 enddo

 ! part <PSI_1|O|PSI_1>
 accujj_diag = 0.d0
 accujj_off_diag = 0.d0
 do i = 1, nsingles
  j = i
  call i_H_j_eff_pot(psi_det(1,1,idx(j)),psi_det(1,1,idx(i)),mo_dipole_z,mo_dipole_z,mo_num,N_int,oij)
  accujj_diag -= oij * coef(i) * coef(j)
  accu_diag -= oij * coef(i) * coef(j)
!  print*,'coef,ojj',coef(i),oij
  do j = i+1, nsingles
   call i_H_j_eff_pot(psi_det(1,1,idx(j)),psi_det(1,1,idx(i)),mo_dipole_z,mo_dipole_z,mo_num,N_int,oij)
!   accujj_off_diag -= 2.d0 * oij * coef(i) * coef(j)
!   accu_off_diag -= 2.d0 * oij * coef(i) * coef(j)
  enddo
 enddo
 z_dip = (accu00 * normpsi0 + accu0j + accujj_diag + accujj_off_diag)/(1.d0 + normpsi1)

 double precision :: nucl_contrib

 nucl_contrib = 0.d0
 do i = 1,nucl_num
  nucl_contrib += nucl_charge(i) * nucl_coord(i,3)
 enddo
 print*,'**************************'
 print*,'**************************'
 print*,'**************************'
 print*,'z_dip           = ',z_dip + nucl_contrib
! print*,'z_dipole_moment = ',z_dipole_moment
 print*,'**************************'
 print*,'electron contrib= ',z_dip
 print*,'nucl_contrib    = ',nucl_contrib
 print*,''
! print*,'z_dip diag      = ',z_dipole_moment_diag(1)
! print*,'z_dip off diag  = ',z_dipole_moment_off_diag(1)
 print*,'accu_diag       = ',accu_diag
 print*,'accu_off_diag   = ',accu_off_diag
 print*,''
 print*,'accu00          = ',accu00 * normpsi0
 print*,'accu0j          = ',accu0j
 print*,'accujj_diag     = ',accujj_diag
 print*,'accujj_off_diag = ',accujj_off_diag
 print*,''
 print*,'e2              =',e2
 print*,'normpsi1        =',normpsi1
 print*,''

end
