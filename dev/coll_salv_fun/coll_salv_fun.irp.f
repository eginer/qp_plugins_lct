program coll_salv_fun
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 read_wf = .True.
!touch read_wf
!  print *, 'Hello world'
!  call test_int_f_tilde
! call test_shank
  call test_int_f_special
! call routine_test_n2_j
! call routine_test_n2_j_full
! call print_energy
end

subroutine test_shank
 implicit none
 integer :: i,n,k
 double precision, allocatable :: array(:)
 double precision :: accu,shank_general,test,pi
 pi = dacos(-1.d0)
 n = 6
 accu = 0.d0
 allocate(array(0:n))
 do i = 0, n
  accu += 4.d0 * (-1.d0)**dble(i)/(2.d0*dble(i)+1.d0)
  array(i) = accu
 enddo
 test = shank_general(array,n,n)
 print*,'shank, error = ',test,dabs(test -pi )
 print*,'accu , error = ',accu,dabs(accu -pi )
 print*,'pi           = ',pi


end

subroutine print_energy
 implicit none
 print*,'************************************************'
 print*,'psi_energy             = ',psi_energy
 print*,'psi_energy_two_e       = ',psi_energy_two_e
 print*,''
 print*,'************************************************'
! print*,'psi_wee_mu_of_r        = ',psi_wee_mu_of_r
 print*,'---------> '
! print*,'psi_wee_mu_of_r_sr     = ',psi_wee_mu_of_r_sr
 print*,''
 print*,'coulomb_n2_jastrow     = ',coulomb_n2_jastrow
 print*,'Delta                  = ',coulomb_n2_jastrow  - psi_energy_two_e
 print*,'0.5 * Delta            = ',0.5d0*(coulomb_n2_jastrow  - psi_energy_two_e)
 print*,'E_FCI + correction     = ',psi_energy + 0.5d0*(coulomb_n2_jastrow  - psi_energy_two_e)
 print*,''
 print*,'coulomb_n2_jastrow_reno= ',coulomb_n2_jastrow_renorm
 print*,'Delta                  = ',coulomb_n2_jastrow_renorm  - psi_energy_two_e
 print*,'0.5 * Delta            = ',0.5d0*(coulomb_n2_jastrow_renorm  - psi_energy_two_e)
 print*,'E_FCI + correction     = ',psi_energy + 0.5d0*(coulomb_n2_jastrow_renorm  - psi_energy_two_e)
 print*,''
 print*,'norm_n2_jastrow        = ',norm_n2_jastrow
 print*,'coulomb_n2_jastrow/norm= ',coulomb_n2_jastrow/norm_n2_jastrow
 print*,'Delta                  = ',coulomb_n2_jastrow/norm_n2_jastrow  - psi_energy_two_e
 print*,'0.5 * Delta            = ',0.5d0*(coulomb_n2_jastrow/norm_n2_jastrow  - psi_energy_two_e)
 print*,'E_FCI + correction     = ',psi_energy + 0.5d0*(coulomb_n2_jastrow/norm_n2_jastrow  - psi_energy_two_e)
 print*,''
 print*,''
 print*,'E_FCI + PBEoTSu        = ',psi_energy + ecmd_pbe_on_top_su_mu_of_r
! print*,'coulomb_n2_jastrow_reno= ',coulomb_n2_jastrow_renorm
! print*,'wee_mu_of_r_n2_jastrow = ',wee_mu_of_r_n2_jastrow
! print*,'wee_mu_of_r_n2_jastrow_r ',wee_mu_of_r_n2_jastrow_renorm
! print*,'---------> '
! print*,'wee_mu_of_r_sr_n2_jastrow',wee_mu_of_r_sr_n2_jastrow
! print*,'wee_mu_of_r_sr_n2Rjastrow',wee_mu_of_r_sr_n2_jastrow_renorm
! print*,''
! print*,''
! print*,'Delta sr               = ',wee_mu_of_r_sr_n2_jastrow - psi_wee_mu_of_r_sr
! print*,'0.5 Delta sr           = ',0.5d0*(wee_mu_of_r_sr_n2_jastrow - psi_wee_mu_of_r_sr)
! print*,''
! print*,'Delta sr renorm        = ',wee_mu_of_r_sr_n2_jastrow_renorm - psi_wee_mu_of_r_sr
! print*,'0.5 Delta sr renorm    = ',0.5d0*(wee_mu_of_r_sr_n2_jastrow_renorm - psi_wee_mu_of_r_sr)


end
