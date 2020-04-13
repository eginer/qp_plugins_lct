!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part) with exact on top pair density
  END_DOC
 
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, rho2
  double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2, decdrho, decdrho2!, decdrho2_a, decdrho2_b
  double precision              :: ecPBE,decPBEdrho_a,decPBEdrho_b,decPBEdgrad_rho_2,decPBEdrho, decPBEdgrad_rho_a_2,decPBEdgrad_rho_b_2,decPBEdgrad_rho_a_b
  double precision              :: rho_c, rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_o_c,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o
  double precision              :: beta, dbetadrho, dbetadgrad_rho_2, denom, ddenomdrho, ddenomdgrad_rho_2, ddenomdrho2
  double precision              :: pi, c, thr
  double precision              :: rho, m  
 
  if(abs(rho_a-rho_b) > 1.d-12)then
  stop "routine implemented only for closed-shell systems"
  endif 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b
  thr = 1.d-12
  
! correlation PBE standard and on-top pair distribution 
  call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_c_2,grad_rho_o_c)

  call ec_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_o_c,ecPBE,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o)

  call v_rho_oc_to_v_rho_ab(decPBEdrho_o, decPBEdrho_c, decPBEdrho_a, decPBEdrho_b)
  call v_grad_rho_oc_to_v_grad_rho_ab(decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_2, decPBEdgrad_rho_c_o, decPBEdgrad_rho_a_2, decPBEdgrad_rho_b_2, decPBEdgrad_rho_a_b)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
   
  beta = ecPBE/(c*rho2)
  if(dabs(beta).lt.thr)then
   beta = 1.d-12
  endif

  denom = 1.d0 + beta*mu**3
  ec_srmuPBE=ecPBE/denom

! calculation of derivatives 
  !dec/dn
  decPBEdrho = 0.5d0 *(decPBEdrho_a + decPBEdrho_b)

  dbetadrho = decPBEdrho/(c*rho2) ! - (ecPBE/(c*rho2**2))*dn2_UEGdrho
  ddenomdrho = dbetadrho*mu**3

  decdrho = decPBEdrho/denom - ecPBE*ddenomdrho/(denom**2)
  decdrho_a = decdrho
  decdrho_b = decdrho

  !dec/((dgradn)^2)
  decPBEdgrad_rho_2 = 0.25d0 *(decPBEdgrad_rho_a_2 + decPBEdgrad_rho_b_2 + decPBEdgrad_rho_a_b) 
 
  dbetadgrad_rho_2 = decPBEdgrad_rho_2/(c*rho2)
  ddenomdgrad_rho_2 = dbetadgrad_rho_2*mu**3
  
  decdgrad_rho_2 = decPBEdgrad_rho_2/denom - ecPBE*ddenomdgrad_rho_2/(denom**2)
  decdgrad_rho_a_2 = decdgrad_rho_2 ! + decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_b_2 = decdgrad_rho_2 ! - decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_a_b = decdgrad_rho_2 ! - decdgrad_m_2 

  !dec/dn2
  
  ddenomdrho2 = - (mu**3)* ecPBE/(c*rho2**2)

  decdrho2 = - ecPBE*ddenomdrho2/(denom**2)
  ! decdrho2_a = decdrho2
  ! decdrho2_b = decdrho2

  end subroutine ecmdsrPBE

!-----------------------------------------------------------------Integrales------------------------------------------------------------------

BEGIN_PROVIDER[double precision, energy_c_md_sr_pbe_n2, (N_states) ]
 implicit none
 BEGIN_DOC
 ! exchange / correlation energies  with the short-range version Perdew-Burke-Ernzerhof GGA functional 
 !
 ! defined in Chem. Phys.329, 276 (2006)
 END_DOC 
 BEGIN_DOC
! exchange/correlation energy with the short range pbe functional
 END_DOC
 integer :: istate,ipoint,j,m
 double precision :: two_dm_in_r_exact
 double precision :: weight, r(3)
 double precision :: ec_srmuPBE, mu
 double precision :: rho2, rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: decdrho_a, decdrho_b, decdrho, decdrho2
 double precision :: decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 energy_c_md_sr_pbe = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)

   weight = final_weight_at_r_vector(ipoint)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   !rho2 = two_dm_in_r_exacti(r,r,istate)
   call give_on_top_in_r_one_state(r,istate,rho2)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

   rho2 = rho2*2.d0 ! normalization
   mu = mu_of_r_prov(ipoint,istate)
   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
   
   decdrho2 = 2.d0*decdrho2 ! normalization

   energy_c_md_sr_pbe_n2(istate) += ec_srmuPBE * weight
  enddo
 enddo

END_PROVIDER

!1
 BEGIN_PROVIDER [double precision, potential_c_alpha_ao_md_sr_pbe_n2,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_md_sr_pbe_n2,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional 
 !
 ! defined in Chem. Phys.329, 276 (2006)
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_c_alpha_ao_md_sr_pbe_n2(j,i,istate) = pot_scal_c_alpha_ao_md_sr_pbe_n2(j,i,istate) + pot_grad_c_alpha_ao_md_sr_pbe_n2(j,i,istate) + pot_grad_c_alpha_ao_md_sr_pbe_n2(i,j,istate)
      potential_c_beta_ao_md_sr_pbe_n2(j,i,istate) = pot_scal_c_beta_ao_md_sr_pbe_n2(j,i,istate) + pot_grad_c_beta_ao_md_sr_pbe_n2(j,i,istate) + pot_grad_c_beta_ao_md_sr_pbe_n2(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 

!3
 BEGIN_PROVIDER[double precision, aos_vc_alpha_md_sr_pbe_w_n2  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_md_sr_pbe_w_n2   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_alpha_md_sr_pbe_w_n2  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_beta_md_sr_pbe_w_n2   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! intermediates to compute the sr_pbe potentials 
 END_DOC
 integer :: istate,ipoint,j,m
 double precision :: weight, r(3)
 double precision :: ec_srmuPBE,mu
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2
 double precision :: contrib_grad_ca(3),contrib_grad_cb(3)
 double precision :: decdrho_a, decdrho_b, decdrho, decdrho2
 double precision :: decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 aos_d_vc_alpha_md_sr_pbe_w_n2 = 0.d0
 aos_d_vc_beta_md_sr_pbe_w_n2 = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   call give_on_top_in_r_one_state(r,istate,rho2)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo
   
   rho2 = 2.d0*rho2
   
   ! mu_erf_dft -> mu_b
   mu = mu_of_r_prov(ipoint,istate)
   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
   
   decdrho_a *= weight
   decdrho_b *= weight

   do m= 1,3
    contrib_grad_ca(m) = weight * (2.d0 * decdgrad_rho_a_2 *  grad_rho_a(m) + decdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_cb(m) = weight * (2.d0 * decdgrad_rho_b_2 *  grad_rho_b(m) + decdgrad_rho_a_b  * grad_rho_a(m) )
   enddo

   do j = 1, ao_num
    aos_vc_alpha_md_sr_pbe_w_n2(j,ipoint,istate) = decdrho_a * aos_in_r_array(j,ipoint)
    aos_vc_beta_md_sr_pbe_w_n2 (j,ipoint,istate) = decdrho_b * aos_in_r_array(j,ipoint)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_d_vc_alpha_md_sr_pbe_w_n2(j,ipoint,istate) += contrib_grad_ca(m) * aos_grad_in_r_array_transp(m,j,ipoint)
     aos_d_vc_beta_md_sr_pbe_w_n2 (j,ipoint,istate) += contrib_grad_cb(m) * aos_grad_in_r_array_transp(m,j,ipoint)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER

!4
 BEGIN_PROVIDER [double precision, pot_scal_c_alpha_ao_md_sr_pbe_n2, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_beta_ao_md_sr_pbe_n2, (ao_num,ao_num,N_states)]
 implicit none
! intermediates to compute the sr_pbe potentials 
! 
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_scal_c_alpha_ao_md_sr_pbe_n2 = 0.d0
   pot_scal_c_beta_ao_md_sr_pbe_n2 = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_alpha_md_sr_pbe_w_n2(1,1,istate),size(aos_vc_alpha_md_sr_pbe_w_n2,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_alpha_ao_md_sr_pbe_n2(1,1,istate),size(pot_scal_c_alpha_ao_md_sr_pbe_n2,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_beta_md_sr_pbe_w_n2(1,1,istate),size(aos_vc_beta_md_sr_pbe_w_n2,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_beta_ao_md_sr_pbe_n2(1,1,istate),size(pot_scal_c_beta_ao_md_sr_pbe_n2,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 

!5
 BEGIN_PROVIDER [double precision, pot_grad_c_alpha_ao_md_sr_pbe_n2,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_c_beta_ao_md_sr_pbe_n2,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_c_alpha_ao_md_sr_pbe_n2 = 0.d0
   pot_grad_c_beta_ao_md_sr_pbe_n2 = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_alpha_md_sr_pbe_w_n2(1,1,istate),size(aos_d_vc_alpha_md_sr_pbe_w_n2,1),  &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_c_alpha_ao_md_sr_pbe_n2(1,1,istate),size(pot_grad_c_alpha_ao_md_sr_pbe_n2,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_beta_md_sr_pbe_w_n2(1,1,istate),size(aos_d_vc_beta_md_sr_pbe_w_n2,1),    &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_c_beta_ao_md_sr_pbe_n2(1,1,istate),size(pot_grad_c_beta_ao_md_sr_pbe_n2,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER
