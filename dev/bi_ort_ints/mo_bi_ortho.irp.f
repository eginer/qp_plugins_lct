 BEGIN_PROVIDER [ double precision, mo_r_coef, (ao_num, mo_num)]
 implicit none
 mo_r_coef = mo_coef
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, mo_l_coef, (ao_num, mo_num)]
 implicit none
 mo_l_coef = mo_coef
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, mo_r_coef_transp, (mo_num, ao_num)]
 implicit none
 mo_r_coef_transp = mo_coef_transp
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, mo_l_coef_transp, (mo_num, ao_num)]
 implicit none
 mo_l_coef_transp = mo_coef_transp
 END_PROVIDER 
