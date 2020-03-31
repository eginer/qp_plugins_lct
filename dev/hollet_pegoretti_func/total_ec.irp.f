subroutine give_ec_holpeg(rho,on_top,lapl_aa,lapl_bb,ec_ab,ec_aa,ec_bb)
 implicit none
 double precision, intent(in) :: rho,on_top,lapl_aa,lapl_bb 
 double precision, intent(out):: ec_ab,ec_aa,ec_bb
 double precision :: V_ab_holl_peg,V_aa_holl_peg
 ec_ab = V_ab_holl_peg(rho,on_top)
 ec_aa = V_aa_holl_peg(lapl_aa)
 ec_bb = V_aa_holl_peg(lapl_bb)

end
