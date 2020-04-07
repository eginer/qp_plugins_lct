 ! routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction
 ! From Emmanuel's plugins: dft_utils_one_e/utils.irp.f
 !
 !-----------------------------------------------------------------------------------------------------------------------------------
 ! Parameter             : Definition                                                                   ; Source
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! a, b and c            :                                                                              ; from (2), eq. ??
 ! B1, C1, D1, E1, F1    :                                                                              ; from (2)
 ! beta                  :                                                                              ; from (2), eq. ??
 ! delta                 :                                                                              ; from (2), eq. ??
 ! dbetadrho             : dbeta/dn                                                                     ;
 ! ddeltadrho            : ddelta/dn                                                                    ;
 ! dgammadrho            : dgamma/dn                                                                    ;
 ! dg0drho               : dg0/dn                                                                       ; from (5), eq. 12
 ! dg0drs                : dg0/drs                                                                      ; from (5), eq. 14
 ! decdrho, decdgrad_rho : decsrmuPBE/dn and decsrmuPBE/dgradrho                                        ; from (3) and (4)
 ! decPBEdrho            : decPBE/dn                                                                    ;
 ! decPBEdgrad_rho_2     :                                                                              ; 
 ! dexdrho, dexdgrad_rho : dexsrmuPBE/dn and dexsrmuPBE/dgradrho                                        ; from (3) and (4)
 ! dexPBEdrho            : dexPBE/dn                                                                    ;
 ! dexPBEdgrad_rho_2     : dexPBE/dgrad(n)^2                                                            ;
 ! ecPBE                 : Usual density of correlation energy                                          ; from (1), already done in QP
 ! ec_srmuPBE            :
 ! exPBE                 : Usual density of exchange energy                                             ; from (1), already done in QP
 ! ex_srmuPBE            :
 ! gamma                 :                                                                              ; from (2), eq. ??
 ! g0                    : On-top pair-distribution function of the spin-unpolarized UEG                ;
 ! g0_UEG_mu_inf         :                                                                              ; rsdft_ecmd/ueg_on_top.irp.f
 ! grad_rho              : gradient of density                                                          ; 
 ! n2_UEG                : On-top pair density of the uniform electron gas                              ; from (2), eq. 51
 ! n2xc_UEG              : On-top exchange-correlation pair density of the uniform electron gas         ; from (2), below eq. 55
 ! rho                   : rho_a + rho_b (both densities of spins alpha and beta)                       ;
 ! rs                    : Seitz radius                                                                 ; rsdft_ecmd/ueg_on_top.irp.f
 !-----------------------------------------------------------------------------------------------------------------------------------
 ! SOURCES
 !-----------------------------------------------------------------------------------------------------------------------------------
 ! (1) : Generalized Gradient Approximation Made Simple - J. P. Perdew, K. Burke, M. Ernzerhof - PRL(77),18 (28/10/96)
 ! (2) : Short-range exchange and correlation density functionnals - J. Toulouse (Odense-Paris collaboration)
 ! (3) : A. Ferté's Thesis
 ! (4) : Developpement of eq. to be done
 ! (5) : Supplementary Materials for 'A density-based basis-set incomleteness correction for GW Methods' 
 !       - P.-F. Loos, B. Pradines, A. Scemama, E. Giner, J. Toulouse - ???????????
 ! 
 ! n = na + nb  | m = na - nb
 ! na = (n+m)/2 | nb = (n-m)/2
 ! gn2 = ga2 + gb2 + 2 gab | ga2 = (gn2 + gm2 + 2gnm)/4
 ! gm2 = ga2 + gb2 - 2 gab | gb2 = (gn2 + gm2 - 2gnm)/4
 ! gnm = ga2 - gb2         | gab = (gn2 - gm2)/4
 ! dedn = dedna * dnadn + dednb * dnbdn = (1/2) * (dedna + dednb)
 ! dedgn2 = dedga2 * dga2/dgn2 + dedgb2 * dgb2dgn2 + dedgab * dgabdgn2 = 
 !----------------------------------------------------------------------------------------------------------------------------------- 

!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
  
  implicit none
  BEGIN_DOC
  ! Give the on-top pair distribution function g0 and its derivative according to rho dg0drho
  END_DOC

  double precision, intent (in) :: rho, rho_a, rho_b
  double precision, intent (out) :: g0, dg0drho
  double precision :: pi
  double precision :: g0_UEG_mu_inf, dg0drs
  double precision :: C1, F1, D1, E1, B1, rs

  pi = dacos(-1.d0)
  C1 = 0.0819306d0
  F1 = 0.752411d0
  D1 = -0.0127713d0
  E1 = 0.00185898d0
  B1 = 0.7317d0 - F1
  rs = (3.d0 / (4.d0*pi*rho))**(1.d0/3.d0) 
   
  g0 = g0_UEG_mu_inf(rho_a, rho_b)
  dg0drs = 0.5d0*((-B1 + 2.d0*C1*rs + 3.d0*D1*rs**2 + 4.d0*E1*rs**3)-F1*(1.d0 - B1*rs + C1*rs**2 + D1*rs**3 + E1*rs**4))*exp(-F1*rs)
  dg0drho = -((6.d0*dsqrt(pi)*rho**2)**(-2.d0/3.d0))*dg0drs
 
  end subroutine g0_dg0 
  
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine exmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex_srmuPBE,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,dexdgrad_rho_2)
  
  implicit none
  BEGIN_DOC
  ! Calculation of exchange energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part)
  END_DOC
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision, intent(out) :: ex_srmuPBE,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2
  double precision              :: exPBE,dexPBEdrho_a,dexPBEdrho_b, dexPBEdrho, dexPBEdgrad_rho_a_2,dexPBEdgrad_rho_b_2,dexPBEdgrad_rho_a_b, dexPBEdgrad_rho_2
  double precision              :: gamma, dgammadrho, dgammadgrad_rho_2, delta, ddeltadrho, ddeltadgrad_rho_2, denom, ddenomdrho, ddenomdgrad_rho_2
  double precision              :: pi, a, b, thr
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho, n2xc_UEG, dn2xc_UEGdrho, g0, dg0drho

  if(abs(rho_a-rho_b) > 1.d-12)then
  stop "routine implemented only for closed-shell systems"
  endif 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b
  thr = 1.d-12

! exchange PBE standard and on-top pair distribution
  call ex_pbe_sr(1.d-12,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,exPBE,dexPBEdrho_a,dexPBEdrho_b,dexPBEdgrad_rho_a_2,dexPBEdgrad_rho_b_2,dexPBEdgrad_rho_a_b)
  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
  
  if(dabs(exPBE).lt.thr)then
   exPBE = 1.d-12
  endif

  if(dabs(rho).lt.thr)then
   rho = 1.d-12
  endif

  if(dabs(g0).lt.thr)then
    g0 = 1.d-12
  endif

! calculation of energy
  a = pi / 2.d0
  b = 2*dsqrt(pi)*(2*dsqrt(2.d0) - 1.d0)/3.d0   
  
  n2_UEG = (rho**2)*g0
  if(dabs(n2_UEG).lt.thr)then
   n2_UEG = 1.d-12
  endif
  
  n2xc_UEG = n2_UEG - rho**2
  if(dabs(n2xc_UEG).lt.thr)then
   n2xc_UEG = 1.d-12
  endif

  gamma = exPBE / (a*n2xc_UEG)
  if(dabs(gamma).lt.thr)then
   gamma = 1.d-12
  endif

  delta = -(b*n2_UEG*gamma**2) / exPBE
  if(dabs(delta).lt.thr)then
   delta = 1.d-12
  endif
  
  denom = 1.d0 + delta*mu + gamma*(mu**2)

  ex_srmuPBE=exPBE/denom

! calculation of derivatives
  !dex/dn
  dexPBEdrho = 0.5d0 *(dexPBEdrho_a + dexPBEdrho_b)
  dn2_UEGdrho = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2xc_UEGdrho = dn2_UEGdrho - 2.d0*rho

  dgammadrho = (1.d0/(a*n2xc_UEG))*dexPBEdrho  - (exPBE/(a*n2xc_UEG**2))*dn2xc_UEGdrho
  ddeltadrho = -((b*gamma**2)/exPBE)*dn2_UEGdrho -(b*n2_UEG/exPBE)*2.d0*gamma*dgammadrho + (b*n2_UEG*gamma**2)/(exPBE**2)*dexPBEdrho

  ddenomdrho = ddeltadrho * mu + dgammadrho * (mu**2)
  dexdrho = dexPBEdrho/denom - exPBE*ddenomdrho/denom**2
  dexdrho_a = dexdrho
  dexdrho_b = dexdrho
 
  !dex/d((gradn)^2)
  dexPBEdgrad_rho_2 = 0.25d0 *(dexPBEdgrad_rho_a_2 + dexPBEdgrad_rho_b_2 + 2.d0*dexPBEdgrad_rho_a_b) !! 2*ab??!
 
  !dexPBEdgrad_rho_2 = 0.25d0 *(dexPBEdgrad_rho_a_2 + dexPBEdgrad_rho_b_2 + dexPBEdgrad_rho_a_b) !! 2*ab??!
 
  dgammadgrad_rho_2 = dexPBEdgrad_rho_2/(a*n2xc_UEG)
  ddeltadgrad_rho_2 = ((b*n2_UEG*gamma**2)/(exPBE**2))*dexPBEdgrad_rho_2 - b*(n2_UEG/exPBE)*2*gamma*dgammadgrad_rho_2

  ddenomdgrad_rho_2 = ddeltadgrad_rho_2*mu + dgammadgrad_rho_2*mu**2
  dexdgrad_rho_2 = dexPBEdgrad_rho_2/denom - exPBE*ddenomdgrad_rho_2/(denom**2)
  dexdgrad_rho_a_2 = dexdgrad_rho_2 ! + decdgrad_n_m + decdgrad_m_2
  dexdgrad_rho_b_2 = dexdgrad_rho_2 ! - decdgrad_n_m + decdgrad_m_2
  dexdgrad_rho_a_b = dexdgrad_rho_2 ! - decdgrad_m_2 
  end subroutine exmdsrPBE
!---------------------------------------------------------------------------------------------------------------------------------------------

  subroutine ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part)
  END_DOC
 
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2, decdrho
  double precision              :: ecPBE,decPBEdrho_a,decPBEdrho_b,decPBEdgrad_rho_2,decPBEdrho, decPBEdgrad_rho_a_2,decPBEdgrad_rho_b_2,decPBEdgrad_rho_a_b
  double precision              :: rho_c, rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_o_c,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o
  double precision              :: beta, dbetadrho, dbetadgrad_rho_2, denom, ddenomdrho, ddenomdgrad_rho_2
  double precision              :: pi, c, thr
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho, n2xc_UEG, dn2xc_UEGdrho, g0, dg0drho
 
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

  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
  
  n2_UEG = (rho**2)*g0
  if(dabs(n2_UEG).lt.thr)then
   n2_UEG = 1.d-12
  endif  
  
  beta = ecPBE/(c*n2_UEG)
  if(dabs(beta).lt.thr)then
   beta = 1.d-12
  endif

  denom = 1.d0 + beta*mu**3
  ec_srmuPBE=ecPBE/denom

! calculation of derivatives 
  !dec/dn
  decPBEdrho = 0.5d0 *(decPBEdrho_a + decPBEdrho_b)

  dn2_UEGdrho = 2.d0*rho*g0 + (rho**2)*dg0drho

  dbetadrho = decPBEdrho/(c*n2_UEG) - (ecPBE/(c*n2_UEG**2))*dn2_UEGdrho
  ddenomdrho = dbetadrho*mu**3

  decdrho = decPBEdrho/denom - ecPBE*ddenomdrho/(denom**2)
  decdrho_a = decdrho
  decdrho_b = decdrho

  !dec/((dgradn)^2)
  decPBEdgrad_rho_2 = 0.25d0 *(decPBEdgrad_rho_a_2 + decPBEdgrad_rho_b_2 + 2.d0*decPBEdgrad_rho_a_b) !! Vérifier le facteur 2
 
  !decPBEdgrad_rho_2 = 0.25d0 *(decPBEdgrad_rho_a_2 + decPBEdgrad_rho_b_2 + decPBEdgrad_rho_a_b) !! Vérifier le facteur 2
 
  dbetadgrad_rho_2 = decPBEdgrad_rho_2/(c*n2_UEG)
  ddenomdgrad_rho_2 = dbetadgrad_rho_2*mu**3
  
  decdgrad_rho_2 = decPBEdgrad_rho_2/denom - ecPBE*ddenomdgrad_rho_2/(denom**2)
  decdgrad_rho_a_2 = decdgrad_rho_2 ! + decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_b_2 = decdgrad_rho_2 ! - decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_a_b = decdgrad_rho_2 ! - decdgrad_m_2 
  end subroutine ecmdsrPBE

!---------------------------------------------------------------------------------------------------------------------------------------------

subroutine exc_dexc_md_sr_PBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdrho,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdrho,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2)
 
 implicit none
 BEGIN_DOC
 ! Give exchange and correlation energies and chemical potentials
 ! Use qp_plugins_dtraore/sr_md_energies/utils_modified_jt.irp.f's plugins : exmdsrPBE and ecmdsrPBE
 END_DOC
 double precision, intent(in)  :: mu
 double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision, intent(out) :: ex_srmuPBE,dexdrho_a,dexdrho_b,dexdrho,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2
 double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdrho,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 call exmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex_srmuPBE,dexdrho_a,dexdrho_b,dexdrho,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,dexdgrad_rho_2)

 call ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b,decdrho,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2)

 end subroutine excmdsrPBE


!-----------------------------------------------------------------Integrales------------------------------------------------------------------
 BEGIN_PROVIDER[double precision, energy_x_md_sr_pbe, (N_states) ]
&BEGIN_PROVIDER[double precision, energy_c_md_sr_pbe, (N_states) ]
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
 double precision :: weight
 double precision :: ex_srmuPBE, ec_srmuPBE
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: decdrho_a, decdrho_b, dexdrho_a, dexdrho_b, dexdrho, decdrho
 double precision :: dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2



 energy_x_md_sr_pbe = 0.d0
 energy_c_md_sr_pbe = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
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

   call exc_dexc_md_sr_PBE(mu_erf_dft,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdrho,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdrho,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2)

   energy_x_md_sr_pbe(istate) += ex_srmuPBE * weight
   energy_c_md_sr_pbe(istate) += ec_srmuPBE * weight
  enddo
 enddo

END_PROVIDER


!1
 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional 
 !
 ! defined in Chem. Phys.329, 276 (2006)
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_x_alpha_ao_md_sr_pbe(j,i,istate) = pot_scal_x_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_x_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_x_alpha_ao_md_sr_pbe(i,j,istate)
      potential_x_beta_ao_md_sr_pbe(j,i,istate) = pot_scal_x_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_x_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_x_beta_ao_md_sr_pbe(i,j,istate)

      potential_c_alpha_ao_md_sr_pbe(j,i,istate) = pot_scal_c_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_c_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_c_alpha_ao_md_sr_pbe(i,j,istate)
      potential_c_beta_ao_md_sr_pbe(j,i,istate) = pot_scal_c_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_c_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_c_beta_ao_md_sr_pbe(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 

!2
 BEGIN_PROVIDER [double precision, potential_xc_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_xc_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_xc_alpha_ao_md_sr_pbe(j,i,istate) = pot_scal_xc_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_xc_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_xc_alpha_ao_md_sr_pbe(i,j,istate)
      potential_xc_beta_ao_md_sr_pbe(j,i,istate)  = pot_scal_xc_beta_ao_md_sr_pbe(j,i,istate)  + pot_grad_xc_beta_ao_md_sr_pbe(j,i,istate)  + pot_grad_xc_beta_ao_md_sr_pbe(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 


!3
 BEGIN_PROVIDER[double precision, aos_vc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_md_sr_pbe_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_beta_md_sr_pbe_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_beta_md_sr_pbe_w   ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vx_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vx_beta_md_sr_pbe_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! intermediates to compute the sr_pbe potentials 
! 
! aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: weight
 double precision :: ex_srmuPBE, ec_srmuPBE
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: contrib_grad_xa(3),contrib_grad_xb(3),contrib_grad_ca(3),contrib_grad_cb(3)
 double precision :: decdrho_a, decdrho_b, dexdrho_a, dexdrho_b, dexdrho, decdrho
 double precision :: dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 aos_d_vc_alpha_md_sr_pbe_w= 0.d0
 aos_d_vc_beta_md_sr_pbe_w = 0.d0
 aos_d_vx_alpha_md_sr_pbe_w= 0.d0
 aos_d_vx_beta_md_sr_pbe_w = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)

   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
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

   call exc_dexc_md_sr_PBE(mu_erf_dft,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdrho,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdrho,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2)

   dexdrho_a *= weight
   decdrho_a *= weight
   dexdrho_b *= weight
   decdrho_b *= weight

   do m= 1,3
    contrib_grad_ca(m) = weight * (2.d0 * decdgrad_rho_a_2 *  grad_rho_a(m) + decdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_xa(m) = weight * (2.d0 * dexdgrad_rho_a_2 *  grad_rho_a(m) + dexdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_cb(m) = weight * (2.d0 * decdgrad_rho_b_2 *  grad_rho_b(m) + decdgrad_rho_a_b  * grad_rho_a(m) )
    contrib_grad_xb(m) = weight * (2.d0 * dexdgrad_rho_b_2 *  grad_rho_b(m) + dexdgrad_rho_a_b  * grad_rho_a(m) )    
   enddo

   do j = 1, ao_num
    aos_vc_alpha_md_sr_pbe_w(j,i,istate) = decdrho_a * aos_in_r_array(j,i)
    aos_vc_beta_md_sr_pbe_w (j,i,istate) = decdrho_b * aos_in_r_array(j,i)
    aos_vx_alpha_md_sr_pbe_w(j,i,istate) = dexdrho_a * aos_in_r_array(j,i)
    aos_vx_beta_md_sr_pbe_w (j,i,istate) = dexdrho_b * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_d_vc_alpha_md_sr_pbe_w(j,i,istate) += contrib_grad_ca(m) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vc_beta_md_sr_pbe_w (j,i,istate) += contrib_grad_cb(m) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vx_alpha_md_sr_pbe_w(j,i,istate) += contrib_grad_xa(m) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vx_beta_md_sr_pbe_w (j,i,istate) += contrib_grad_xb(m) * aos_grad_in_r_array_transp(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER

!4
 BEGIN_PROVIDER [double precision, pot_scal_x_alpha_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_alpha_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_x_beta_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_beta_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
 implicit none
! intermediates to compute the sr_pbe potentials 
! 
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_scal_c_alpha_ao_md_sr_pbe = 0.d0
   pot_scal_x_alpha_ao_md_sr_pbe = 0.d0
   pot_scal_c_beta_ao_md_sr_pbe = 0.d0
   pot_scal_x_beta_ao_md_sr_pbe = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_alpha_md_sr_pbe_w(1,1,istate),size(aos_vc_alpha_md_sr_pbe_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_alpha_ao_md_sr_pbe(1,1,istate),size(pot_scal_c_alpha_ao_md_sr_pbe,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_beta_md_sr_pbe_w(1,1,istate),size(aos_vc_beta_md_sr_pbe_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_beta_ao_md_sr_pbe(1,1,istate),size(pot_scal_c_beta_ao_md_sr_pbe,1))
     ! exchange alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vx_alpha_md_sr_pbe_w(1,1,istate),size(aos_vx_alpha_md_sr_pbe_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_x_alpha_ao_md_sr_pbe(1,1,istate),size(pot_scal_x_alpha_ao_md_sr_pbe,1))
     ! exchange beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vx_beta_md_sr_pbe_w(1,1,istate),size(aos_vx_beta_md_sr_pbe_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_x_beta_ao_md_sr_pbe(1,1,istate), size(pot_scal_x_beta_ao_md_sr_pbe,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 

!5
 BEGIN_PROVIDER [double precision, pot_grad_x_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_x_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_c_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_c_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_c_alpha_ao_md_sr_pbe = 0.d0
   pot_grad_x_alpha_ao_md_sr_pbe = 0.d0
   pot_grad_c_beta_ao_md_sr_pbe = 0.d0
   pot_grad_x_beta_ao_md_sr_pbe = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_alpha_md_sr_pbe_w(1,1,istate),size(aos_d_vc_alpha_md_sr_pbe_w,1),  &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_c_alpha_ao_md_sr_pbe(1,1,istate),size(pot_grad_c_alpha_ao_md_sr_pbe,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_beta_md_sr_pbe_w(1,1,istate),size(aos_d_vc_beta_md_sr_pbe_w,1),    &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_c_beta_ao_md_sr_pbe(1,1,istate),size(pot_grad_c_beta_ao_md_sr_pbe,1))
       ! exchange alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vx_alpha_md_sr_pbe_w(1,1,istate),size(aos_d_vx_alpha_md_sr_pbe_w,1),  &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_x_alpha_ao_md_sr_pbe(1,1,istate),size(pot_grad_x_alpha_ao_md_sr_pbe,1))
       ! exchange beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vx_beta_md_sr_pbe_w(1,1,istate),size(aos_d_vx_beta_md_sr_pbe_w,1),    &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_x_beta_ao_md_sr_pbe(1,1,istate),size(pot_grad_x_beta_ao_md_sr_pbe,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

!6
 BEGIN_PROVIDER[double precision, aos_vxc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]  ! sr_pbe
&BEGIN_PROVIDER[double precision, aos_vxc_beta_md_sr_pbe_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vxc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vxc_beta_md_sr_pbe_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: weight
 double precision :: ex_srmuPBE, ec_srmuPBE
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: contrib_grad_xa(3),contrib_grad_xb(3),contrib_grad_ca(3),contrib_grad_cb(3)
 double precision :: decdrho_a, decdrho_b, dexdrho_a, dexdrho_b, dexdrho, decdrho
 double precision :: dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2
 
 aos_d_vxc_alpha_md_sr_pbe_w = 0.d0
 aos_d_vxc_beta_md_sr_pbe_w = 0.d0

 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
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

   call exc_dexc_md_sr_PBE(mu_erf_dft,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdrho,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdrho,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2)
 
   dexdrho_a *= weight
   decdrho_a *= weight
   dexdrho_b *= weight
   decdrho_b *= weight
   do m= 1,3
    contrib_grad_ca(m) = weight * (2.d0 * decdgrad_rho_a_2 *  grad_rho_a(m) + decdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_xa(m) = weight * (2.d0 * dexdgrad_rho_a_2 *  grad_rho_a(m) + dexdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_cb(m) = weight * (2.d0 * decdgrad_rho_b_2 *  grad_rho_b(m) + decdgrad_rho_a_b  * grad_rho_a(m) )
    contrib_grad_xb(m) = weight * (2.d0 * dexdgrad_rho_b_2 *  grad_rho_b(m) + dexdgrad_rho_a_b  * grad_rho_a(m) )
   enddo
   do j = 1, ao_num
    aos_vxc_alpha_md_sr_pbe_w(j,i,istate) = ( decdrho_a + dexdrho_a ) * aos_in_r_array(j,i)
    aos_vxc_beta_md_sr_pbe_w (j,i,istate) = ( decdrho_b + dexdrho_b ) * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_d_vxc_alpha_md_sr_pbe_w(j,i,istate) += ( contrib_grad_ca(m) + contrib_grad_xa(m) ) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vxc_beta_md_sr_pbe_w (j,i,istate) += ( contrib_grad_cb(m) + contrib_grad_xb(m) ) * aos_grad_in_r_array_transp(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER

!7
 BEGIN_PROVIDER [double precision, pot_scal_xc_alpha_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_xc_beta_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_scal_xc_alpha_ao_md_sr_pbe = 0.d0
   pot_scal_xc_beta_ao_md_sr_pbe = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! exchange - correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                 aos_vxc_alpha_md_sr_pbe_w(1,1,istate),size(aos_vxc_alpha_md_sr_pbe_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                        &
                 pot_scal_xc_alpha_ao_md_sr_pbe(1,1,istate),size(pot_scal_xc_alpha_ao_md_sr_pbe,1))
     ! exchange - correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                 aos_vxc_beta_md_sr_pbe_w(1,1,istate),size(aos_vxc_beta_md_sr_pbe_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                        &
                 pot_scal_xc_beta_ao_md_sr_pbe(1,1,istate),size(pot_scal_xc_beta_ao_md_sr_pbe,1))
   enddo
 call wall_time(wall_2)

END_PROVIDER 

!8
 BEGIN_PROVIDER [double precision, pot_grad_xc_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_xc_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_xc_alpha_ao_md_sr_pbe = 0.d0
   pot_grad_xc_beta_ao_md_sr_pbe = 0.d0
   do istate = 1, N_states
       ! exchange - correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                      &
                  aos_d_vxc_alpha_md_sr_pbe_w(1,1,istate),size(aos_d_vxc_alpha_md_sr_pbe_w,1), &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,            &
                  pot_grad_xc_alpha_ao_md_sr_pbe(1,1,istate),size(pot_grad_xc_alpha_ao_md_sr_pbe,1))
       ! exchange - correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                      &
                  aos_d_vxc_beta_md_sr_pbe_w(1,1,istate),size(aos_d_vxc_beta_md_sr_pbe_w,1),   &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,            &
                  pot_grad_xc_beta_ao_md_sr_pbe(1,1,istate),size(pot_grad_xc_beta_ao_md_sr_pbe,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

