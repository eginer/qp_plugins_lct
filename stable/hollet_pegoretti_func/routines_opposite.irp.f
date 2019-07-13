BEGIN_PROVIDER [double precision, q_holl_peg ]
 implicit none
 BEGIN_DOC
! q parameter of J. Chem. Phys. 148, 164111 (2018), needed for the opposite spin correlation 
 END_DOC
 q_holl_peg = 2.54d0

END_PROVIDER 

double precision function lambda_rho(q,rho)
 implicit none
 double precision, intent(in) :: rho,q
 lambda_rho = q * rho**0.333333333333d0
end

double precision function V_ab_holl_peg(rho_in,on_top_in)
 implicit none
 BEGIN_DOC 
! equation 25 of J. Chem. Phys. 148, 164111 (2018)
! note : the value of the on-top is such that the two-body density integrated gives N_a * N_b
 END_DOC
 include 'utils/constants.include.F'
 double precision, intent(in) :: rho_in,on_top_in
 double precision             :: rho,on_top,thr_rho,thr_n2
 double precision :: lamb,lamb_2,inv_2lamb,inv_2lamb_2,dexp_inv_2lamb_2,lambda_rho
 double precision :: big_1,big_21,big_22,big_2,num,denom
 double precision :: smal_1,smal_2
 thr_rho = 1.d-10
 thr_n2  = 1.d-30
 V_ab_holl_peg = 0.d0
 rho = max(rho_in,thr_rho) 
 on_top = max(on_top_in,thr_n2) 
 if(rho.lt.thr_rho)then
  return
 endif
 lamb              = lambda_rho(q_holl_peg,rho)
 lamb_2            = lamb * lamb 
 if(lamb.lt.thr_rho)then
  return
 endif
 inv_2lamb         = 0.5d0/lamb
 inv_2lamb_2       = inv_2lamb * inv_2lamb
 dexp_inv_2lamb_2  = dexp(-inv_2lamb_2)
 !!! numerator of equation 25 
 ! first term in equation 25 
 big_1             = 2.d0 * lamb * ( sqpi * lamb - 1.d0 ) * dexp_inv_2lamb_2
 ! second term in 25 
 big_21            = sqpi * (sqpi * lamb - 2.d0 * lamb_2  -1.d0 )
 big_22            = 1.d0 + derf(inv_2lamb)
 big_2             = big_21 * big_22
 ! numerator 
 num               = 2.d0 * pi * on_top * (big_1 + big_2)
 smal_1            = 2.d0 * lamb * dexp_inv_2lamb_2
 smal_2            = sqpi * (1.d0 + 2.d0 * lamb_2) * (1.d0 + derf(inv_2lamb))
 denom             = lamb_2 * (smal_1 + smal_2)
 if(denom.gt.thr_rho)then
  V_ab_holl_peg     = num/denom
 else
  return
 endif
end
