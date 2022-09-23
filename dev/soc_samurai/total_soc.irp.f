BEGIN_PROVIDER [ complex*8, mo_v_soc_tot, (mo_num, mo_num,3)]
 implicit none
 BEGIN_DOC
! mo_v_soc_tot(k,l,mu) = <MO_k | V_soc_1_e ^\mu + V_soc_2_e ^\mu | MO_l>, 
! 
! with mu = 1 ::> V^- = V^x + i V^y, 
!      
!      mu = 2 ::> V^+ = V^x - i V^y,
!
!      mu = 3 ::> V^z 
 END_DOC
 mo_v_soc_tot = mo_one_e_soc

END_PROVIDER 
