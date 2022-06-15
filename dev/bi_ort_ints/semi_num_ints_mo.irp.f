
BEGIN_PROVIDER [ double precision, mo_v_ki_bi_ortho_erf_rk_cst_mu_naive, ( mo_num, mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_k(r) phi_i(r) (erf(mu(R) |r - R|) - 1 )/(2|r - R|) on the MO basis 
! 
! WHERE phi_k(r) is a LEFT MOs and phi_i(r) is a RIGHT MO 
!
! WARNING :: computed in a VERY NAIVE WAY ::> ONLY FOR DEBUGING 
 END_DOC
 integer :: i,k,p,q,ipoint
 do ipoint = 1, n_points_final_grid
  mo_v_ki_bi_ortho_erf_rk_cst_mu_naive(:,:,ipoint) = 0.d0
  do i = 1, mo_num
   do k = 1, mo_num
    do p = 1, ao_num
     do q = 1, ao_num
      mo_v_ki_bi_ortho_erf_rk_cst_mu_naive(k,i,ipoint) += mo_l_coef(p,k) * 0.5d0 * v_ij_erf_rk_cst_mu(q,p,ipoint) * mo_r_coef(q,i)
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_v_ki_bi_ortho_erf_rk_cst_mu, ( mo_num, mo_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_k(r) phi_i(r) (erf(mu(R) |r - R|) - 1 )/(2|r - R|) on the MO basis 
! 
! WHERE phi_k(r) is a LEFT MOs and phi_i(r) is a RIGHT MO
 END_DOC
 integer :: ipoint
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (ipoint) & 
 !$OMP SHARED (n_points_final_grid,v_ij_erf_rk_cst_mu,mo_v_ki_bi_ortho_erf_rk_cst_mu)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   call ao_to_mo_bi_ortho(v_ij_erf_rk_cst_mu(1,1,ipoint),size(v_ij_erf_rk_cst_mu,1),mo_v_ki_bi_ortho_erf_rk_cst_mu(1,1,ipoint),size(mo_v_ki_bi_ortho_erf_rk_cst_mu,1))
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 mo_v_ki_bi_ortho_erf_rk_cst_mu = mo_v_ki_bi_ortho_erf_rk_cst_mu * 0.5d0
END_PROVIDER 


BEGIN_PROVIDER [ double precision, mo_v_ki_bi_ortho_erf_rk_cst_mu_transp, ( n_points_final_grid,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/(2|r - R|) on the MO basis
 END_DOC
 integer :: ipoint,i,j
 do i = 1, mo_num
  do j = 1, mo_num
   do ipoint = 1, n_points_final_grid
    mo_v_ki_bi_ortho_erf_rk_cst_mu_transp(ipoint,j,i) = mo_v_ki_bi_ortho_erf_rk_cst_mu(j,i,ipoint)
   enddo
  enddo
 enddo
! FREE mo_v_ki_bi_ortho_erf_rk_cst_mu
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_x_v_ki_bi_ortho_erf_rk_cst_mu_naive, ( mo_num, mo_num,3,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr x * phi_k(r) phi_i(r) (erf(mu(R) |r - R|) - 1 )/(2|r - R|) on the MO basis 
! 
! WHERE phi_k(r) is a LEFT MOs and phi_i(r) is a RIGHT MO 
!
! WARNING :: computed in a VERY NAIVE WAY ::> ONLY FOR DEBUGING 
 END_DOC
 integer :: i,k,p,q,ipoint,m
 do ipoint = 1, n_points_final_grid
  mo_x_v_ki_bi_ortho_erf_rk_cst_mu_naive(:,:,:,ipoint) = 0.d0
  do i = 1, mo_num
   do k = 1, mo_num
    do m = 1, 3
     do p = 1, ao_num
      do q = 1, ao_num
       mo_x_v_ki_bi_ortho_erf_rk_cst_mu_naive(k,i,m,ipoint) += mo_l_coef(p,k) * 0.5d0 * x_v_ij_erf_rk_cst_mu_transp(q,p,m,ipoint) * mo_r_coef(q,i)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_x_v_ki_bi_ortho_erf_rk_cst_mu, ( mo_num, mo_num,3,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr x * phi_i(r) phi_j(r) (erf(mu(R) |r - R|) - 1)/2|r - R| on the MO basis
 END_DOC
 integer :: ipoint,m
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (ipoint,m) & 
 !$OMP SHARED (n_points_final_grid,x_v_ij_erf_rk_cst_mu_transp,mo_x_v_ki_bi_ortho_erf_rk_cst_mu)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
  do m = 1, 3
   call ao_to_mo_bi_ortho(x_v_ij_erf_rk_cst_mu_transp(1,1,m,ipoint),size(x_v_ij_erf_rk_cst_mu_transp,1),mo_x_v_ki_bi_ortho_erf_rk_cst_mu(1,1,m,ipoint),size(mo_x_v_ki_bi_ortho_erf_rk_cst_mu,1))
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 mo_x_v_ki_bi_ortho_erf_rk_cst_mu = 0.5d0 * mo_x_v_ki_bi_ortho_erf_rk_cst_mu

END_PROVIDER 


BEGIN_PROVIDER [ double precision, x_W_ki_bi_ortho_erf_rk, ( n_points_final_grid,3,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! W_ki_bi_ortho^X(R) = \int dr phi_k(r) (1 - erf(mu |r-R|)) (x-X) phi_i(r) 
 END_DOC
 include 'constants.include.F'
 integer :: ipoint,m,i,k
 double precision :: xyz,cst
 double precision :: wall0, wall1

 cst = 0.5d0 * inv_sq_pi
 print*,'providing x_W_ki_bi_ortho_erf_rk ...'
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (ipoint,m,i,k,xyz) & 
 !$OMP SHARED (x_W_ki_bi_ortho_erf_rk,n_points_final_grid,mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp,mo_v_ki_bi_ortho_erf_rk_cst_mu_transp,mo_num,final_grid_points) 
 !$OMP DO SCHEDULE (dynamic)
 do i = 1, mo_num
  do k = 1, mo_num
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
     xyz = final_grid_points(m,ipoint)
     x_W_ki_bi_ortho_erf_rk(ipoint,m,k,i)  =  mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp(ipoint,m,k,i) - xyz * mo_v_ki_bi_ortho_erf_rk_cst_mu_transp(ipoint,k,i)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
! FREE mo_v_ki_bi_ortho_erf_rk_cst_mu_transp 
! FREE mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp
 call wall_time(wall1)
 print*,'time to provide x_W_ki_bi_ortho_erf_rk = ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp, (n_points_final_grid,3, mo_num, mo_num)]
 implicit none
 integer :: i,j,m,ipoint
 do i = 1, mo_num
  do j = 1, mo_num
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
     mo_x_v_ki_bi_ortho_erf_rk_cst_mu_transp(ipoint,m,j,i) = mo_x_v_ki_bi_ortho_erf_rk_cst_mu(j,i,m,ipoint)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER 
