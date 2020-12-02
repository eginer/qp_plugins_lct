program test_pot
 implicit none
 read_wf = .True.
 touch read_wf
 print*,'e_c_md_basis_pbe_ueg = ',e_c_md_basis_pbe_ueg
 double precision :: de_barth
 call compute_func_der(delta_gamma_i_j_alpha, delta_gamma_i_j_beta, potential_c_alpha_mo_md_sr_pbe, potential_c_beta_mo_md_sr_pbe, de_barth)
 print*,'ec_pbe_at_n              = ',ec_pbe_at_n
 print*,'ec_pbe_at_n_plus_delta_n = ',ec_pbe_at_n_plus_delta_n
 print*,'De                       = ',ec_pbe_at_n_plus_delta_n - ec_pbe_at_n
 print*,'-----'
 print*,'int_vc_pbe_at_n          = ',int_vc_pbe_at_n
 print*,'-----'
 print*,'de_barth                 = ',de_barth
end

subroutine routine_mo
 implicit none
 integer :: i,j
 double precision :: err, accu, thr
 thr = 1.d-6
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   err = dabs(pot_basis_alpha_mo_sp_pbe_ueg(j,i,1) - potential_c_alpha_mo_md_sr_pbe(j,i))
   if(err.gt.thr)then
    print*,'j,i',j,i
    print*,pot_basis_alpha_mo_sp_pbe_ueg(j,i,1),potential_c_alpha_mo_md_sr_pbe(j,i),err
   endif
   accu += err
  enddo
 enddo
 print*,'err = ',err

end

subroutine routine_ao
 implicit none
 integer :: i,j
 double precision :: err, accu, thr
 thr = 1.d-6
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   err = dabs(pot_basis_alpha_ao_sp_pbe_ueg(j,i,1) - potential_c_alpha_ao_md_sr_pbe(j,i,1))
   if(err.gt.thr)then
    print*,'j,i',j,i
    print*,pot_basis_alpha_ao_sp_pbe_ueg(j,i,1),potential_c_alpha_ao_md_sr_pbe(j,i,1),err
    print*,pot_basis_alpha_ao_sp_pbe_ueg(j,i,1)/potential_c_alpha_ao_md_sr_pbe(j,i,1)
   endif
   accu += err
  enddo
 enddo
 print*,'err = ',err

end
