program ecmd_ontop_test
 implicit none
 BEGIN_DOC
 !
 END_DOC

 integer           :: istate, i
 double precision  :: mu
 double precision  :: r(3), r_norm, weight
 double precision  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, rho2, two_dm_in_r_exact
 double precision  :: grad_rho_ax, grad_rho_ay, grad_rho_az
 double precision  :: grad_rho_bx, grad_rho_by, grad_rho_bz
 double precision  :: ec_srmuPBE,decdrho_a,decdrho_b,decdrho,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2, decdrho2
 double precision  :: ec_md_sr_pbe_n2(N_states)
!--------A supprimer après les test-------
 double precision  :: rho, grad_rho_2, mu_array(20), test_r, r_norm_prec(n_points_final_grid), mu_tab(n_points_final_grid)
 integer :: p, k
!-----------------------------------------
! mu = 0.5d0
 r_norm_prec = 1.d-12
 mu_tab = 1.d-12
 mu_array = (/ 0.d0, 0.125d0, 0.25d0, 0.375d0, 0.5d0, 0.625d0, 0.75d0, 0.875d0, 1.d0, 1.5d0, 2.d0, 2.5d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0 /)
print*, '#r_norm, rho, rho2, grad_rho_2, mu, ec_srmuPBE, decdrho, decdgrad_rho_2, decdrho2'

do p = 1, 20  ! loop over mu_array  !do1 
! print*, mu_array(p)
 ec_md_sr_pbe_n2 = 1.d-12
 mu = mu_array(p)
 r_norm_prec = 1.d-12
 do i=1, n_points_final_grid !do2
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  r_norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
  weight = final_weight_at_r_vector(i)

  do istate=1, N_states !do3
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   rho = rho_a + rho_b
   !rho2 = two_dm_in_r_exact(r,r,istate)
   call give_on_top_in_r_one_state(r,istate,rho2)
   grad_rho_ax = one_e_dm_and_grad_alpha_in_r(1,i,istate)
   grad_rho_ay = one_e_dm_and_grad_alpha_in_r(2,i,istate)
   grad_rho_az = one_e_dm_and_grad_alpha_in_r(3,i,istate)
   grad_rho_bx = one_e_dm_and_grad_alpha_in_r(1,i,istate)
   grad_rho_by = one_e_dm_and_grad_alpha_in_r(2,i,istate)
   grad_rho_bz = one_e_dm_and_grad_alpha_in_r(3,i,istate)  
   
   grad_rho_a_2 = grad_rho_ax * grad_rho_ax + grad_rho_ay * grad_rho_ay + grad_rho_az * grad_rho_az
   grad_rho_b_2 = grad_rho_bx * grad_rho_bx + grad_rho_by * grad_rho_by + grad_rho_bz * grad_rho_bz
   grad_rho_a_b = grad_rho_ax * grad_rho_bx + grad_rho_ay * grad_rho_by + grad_rho_az * grad_rho_bz
   grad_rho_2 = grad_rho_a_2 + grad_rho_b_2 + 2.d0*grad_rho_a_b
  
  call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

  ec_md_sr_pbe_n2(istate) += ec_srmuPBE * weight
     test_r = 0
     do k=1, i 
        if ((r_norm - r_norm_prec(k)) < 1.d-10) then
            test_r = 1
           exit
        endif
     enddo
     r_norm_prec(i) = r_norm
     mu_tab(i) = mu
     if(test_r==0)then
     print*, r_norm, rho, rho2, grad_rho_2, mu, ec_srmuPBE, decdrho, decdgrad_rho_2, decdrho2
 !     print*, r_norm, rho2
      endif
enddo !enddo3
enddo !enddo 2
!print*, 'here 3'
! print*, mu, ec_md_sr_pbe_n2(1), energy_c_md_sr_pbe_n2(1)
enddo ! enddo 1

end program
