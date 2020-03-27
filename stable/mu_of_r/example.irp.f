
subroutine test_f_HF_valence_ab
 implicit none
 BEGIN_DOC
! routine to test the function f_HF(r1,r2) 
!
! the integral over r1,r2 should be equal to the alpha/beta interaction of HF determinant
 END_DOC
 integer :: ipoint,i,j,i_i,j_j,jpoint
 double precision :: accu_val,accu_ful, weight1,weight2, r1(3),integral_psi_val,integral_psi,r2(3),two_bod
 accu_2 = 0.d0
 ! You compute the coulomb repulsion between alpha-beta electrons for HF
 do i = 1, n_occ_val_orb_for_hf(1)
  i_i = list_valence_orb_for_hf(i,1)
  do j = 1, n_occ_val_orb_for_hf(2)
   j_j = list_valence_orb_for_hf(j,2)
   accu_2 += mo_two_e_integrals_jj(j_j,i_i)
  enddo
 enddo
 print*,''
 print*,''
 print*,''
 print*,'**************************'
 print*,'**************************'
 print*,'Routine to test the f_HF(r1,r2) function'
 print*,'**************************'
 print*,''
 print*,''
 print*,''
 print*,'**************************'
 print*,'<HF| We_ee^{ab}|HF>     = ',accu_2
 print*,'**************************'

 print*,'semi analytical form '
 accu_val = 0.d0
 ! You integrate on r2 the analytical integral over r1 of f_HF(r1,r2) 
 do ipoint  = 1, n_points_final_grid
  weight1 =final_weight_at_r_vector(ipoint)
  r2(1)   = final_grid_points(1,ipoint)
  r2(2)   = final_grid_points(2,ipoint)
  r2(3)   = final_grid_points(3,ipoint)
  call integral_f_HF_valence_ab(r2,integral_psi_val)
  accu_val += integral_psi_val * weight1
 enddo
 print*,'**************************'
 ! Should give you the alpha-beta repulsion of HF, excluding core contributions,  
 print*,'int dr1 dr2 f_HF(r1,r2) = ',accu_val
 double precision :: accu_2


 print*,'pure numerical form (might take quite some time ...)'
 ! You integrate brut force on r1 and r2 
 accu_val = 0.d0
 do jpoint = 1, n_points_final_grid
  weight1 =final_weight_at_r_vector(jpoint)
  r1(1)   = final_grid_points(1,jpoint)
  r1(2)   = final_grid_points(2,jpoint)
  r1(3)   = final_grid_points(3,jpoint)
  do ipoint  = 1, n_points_final_grid
   weight2 =final_weight_at_r_vector(ipoint)
   r2(1)   = final_grid_points(1,ipoint)
   r2(2)   = final_grid_points(2,ipoint)
   r2(3)   = final_grid_points(3,ipoint)
   call f_HF_valence_ab(r1,r2,integral_psi_val,two_bod)
   accu_val += integral_psi_val * weight1 * weight2
  enddo
 enddo
 print*,'int dr1 dr2 f_HF(r1,r2) = ',accu_val


 print*,'**************************'
 print*,'**************************'
 print*,'**************************'
 accu_val = 0.d0
 r1 = 0.d0
 r1(1) = 0.5d0
 print*,'r1 = ',r1
 ! You compute the integral over r2 of f_HF(r1,r2) 
 call integral_f_HF_valence_ab(r1,integral_psi)
 do ipoint  = 1, n_points_final_grid
  weight1 =final_weight_at_r_vector(ipoint)
  r2(1)   = final_grid_points(1,ipoint)
  r2(2)   = final_grid_points(2,ipoint)
  r2(3)   = final_grid_points(3,ipoint)
  call f_HF_valence_ab(r1,r2,integral_psi_val,two_bod)
  accu_val += integral_psi_val * weight1
 enddo
 print*,'int dr2 f_HF(r1,r2)     = ',integral_psi
 print*,'analytical form         = ',accu_val
 print*,'**************************'
end


