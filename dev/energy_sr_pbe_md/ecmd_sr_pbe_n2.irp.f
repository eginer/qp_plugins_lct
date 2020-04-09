
 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_alpha_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! qp_plugins_lct/garbage/no_omp_2rdm/all_states_prov.irp.f
! all_states_act_two_rdm_alpha_beta_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta electrons 
! 
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! !!!!! WARNING !!!!! For efficiency reasons, electron 1 is alpha, electron 2 is beta
!
!  all_states_act_two_rdm_alpha_beta_mo(i,j,k,l,istate) = i:alpha, j:beta, j:alpha, l:beta
!                      
!                      Therefore you don't necessayr have symmetry between electron 1 and 2 
 END_DOC 
 integer :: ispin
 double precision :: wall_1, wall_2
 ! condition for alpha/beta spin
 call wall_time(wall_1)
 print*,'providing all_states_act_two_rdm_alpha_beta_mo ...'
 ispin = 3 
 print*,'ispin = ',ispin
 all_states_act_two_rdm_alpha_beta_mo = 0.d0
 call wall_time(wall_1)
 call orb_range_all_states_two_rdm(all_states_act_two_rdm_alpha_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 call wall_time(wall_2)
 print*,'time to provide all_states_act_two_rdm_alpha_beta_mo',wall_2 - wall_1
 END_PROVIDER 


double precision function two_dm_in_r_exact(r1,r2,istate)
 implicit none
 BEGIN_DOC
 ! two body density evaluated at two points in real space 
 END_DOC
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 allocate(mos_array_r2(mo_num), mos_array_r1(mo_num))
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 two_dm_in_r_exact = 0.d0
 do l = 1, mo_num
  do k = 1, mo_num
    do j = 1, mo_num
     do i = 1, mo_num
     !                                                   1 2 1 2 
     two_dm_in_r_exact += all_states_act_two_rdm_alpha_beta_mo(i,j,k,l,istate) * mos_array_r1(i) * mos_array_r1(k) * mos_array_r2(j) * mos_array_r2(l)
    enddo
   enddo
  enddo
 enddo
 two_dm_in_r_exact = max(two_dm_in_r_exact,1.d-15)
end
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
  decPBEdgrad_rho_2 = 0.25d0 *(decPBEdgrad_rho_a_2 + decPBEdgrad_rho_b_2 + 2.d0*decPBEdgrad_rho_a_b) !! Vérifier le facteur 2 
 
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
 integer :: istate,i,j,m
 double precision :: two_dm_in_r_exact
 double precision :: weight, r(3)
 double precision :: ec_srmuPBE
 double precision :: rho2, rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: decdrho_a, decdrho_b, decdrho, decdrho2
 double precision :: decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 energy_c_md_sr_pbe = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)

   weight = final_weight_at_r_vector(i)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   rho2 = two_dm_in_r_exact(r,r,istate)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,i,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,i,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

  call ecmdsrPBEn2(mu_erf_dft,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

   energy_c_md_sr_pbe_n2(istate) += ec_srmuPBE * weight
  enddo
 enddo

END_PROVIDER
