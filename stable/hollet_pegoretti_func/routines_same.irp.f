BEGIN_PROVIDER [double precision, k_holl_peg]
 implicit none
 BEGIN_DOC
! k parameter of J. Chem. Phys. 148, 164111 (2018), needed for the same spin correlation 
 END_DOC
 k_holl_peg = 1.76d0

END_PROVIDER 

double precision function delta_l0(k,l0)
 implicit none
 double precision, intent(in) :: k,l0
 delta_l0 = k * l0**(1.d0/8.d0)
end

double precision function V_aa_holl_peg(l0)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: l0
  
 double precision :: thr_l0,thr_delta
 double precision :: sqpi3,delta_l0,delta,delta4,inv_delta4,inv_detla4_2
 double precision :: delta_2,delta_3,delta_4,dexp_inv_detla4_2
 thr_l0    = 1.d-20
 thr_delta = 1.d-20
 V_aa_holl_peg = 0.d0
 if(l0.lt.thr_l0)then
  return 
 endif
 sqpi3             = 3.d0 * sqpi  
 delta             = delta_l0(k_holl_peg,l0)
 if(dabs(delta).lt.thr_delta)then
  return 
 endif
 delta4            = 4.d0 * delta
 if(delta4.lt.thr_delta)then
  return 
 endif
 delta_2           = delta * delta
 delta_3           = delta * delta_2
 delta_4           = delta * delta_3
 inv_delta4        = 1.d0/delta4 
 inv_detla4_2      = inv_delta4 * inv_delta4
 dexp_inv_detla4_2 = dexp(-inv_detla4_2)
 ! first term 
 double precision :: big_1,big_21,big_22,big_2,numerator
 big_1             = sqpi3 * (delta * 16.d0 * delta_3) - 40.d0 * delta_2 - 1.d0
 big_1             = big_1 * delta4 * dexp_inv_detla4_2 
 
 big_21            = sqpi3 * (delta * 24.d0 * delta_3) - 48.d0 * (delta_2 + 4.d0 * delta_4) - 1.d0
 big_22            = sqpi * (1.d0 + derf(inv_delta4))
 big_2             = big_21 * big_22  
 numerator         = pi * l0 * (big_2 + big_1)
 double precision :: smal_1, smal_21, smal_22, smal_2, denom
 smal_1            = 4.d0 * (delta + 40.d0 * delta_3) * dexp_inv_detla4_2
 smal_21           = sqpi * (1.d0 + 48.d0 * (delta_2 + 4.d0 * delta_4))
 smal_22           = 1.d0 + derf(inv_delta4)
 smal_2            = smal_22 * smal_21
 denom             = 3.d0 * delta_4 * ( smal_2 + smal_1  )
end
