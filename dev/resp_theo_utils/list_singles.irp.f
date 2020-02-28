BEGIN_PROVIDER [integer, n_occ_resp]
 implicit none
 n_occ_resp = (elec_alpha_num - n_core_orb) 
END_PROVIDER 

BEGIN_PROVIDER [integer, list_occ_resp, (n_occ_resp)]
 implicit none
 integer :: i ,j
 j = 0
 do i = n_core_orb + 1, elec_alpha_num
  j += 1
  list_occ_resp(j) = i
 enddo
END_PROVIDER 

BEGIN_PROVIDER [integer, n_virt_resp ]
 implicit none
 n_virt_resp =  (n_act_orb - elec_alpha_num - n_core_orb)
END_PROVIDER 

BEGIN_PROVIDER [integer, list_virt_resp, (n_virt_resp)]
 implicit none
 integer :: i ,j
 j = 0
 do i = elec_alpha_num + 1, n_act_orb
  j += 1
  list_virt_resp(j) = i
 enddo
END_PROVIDER 

BEGIN_PROVIDER [integer, n_singles_resp]
 implicit none
 integer :: i,j
 n_singles_resp = n_occ_resp * n_virt_resp
END_PROVIDER 

 BEGIN_PROVIDER [integer, list_singles, (n_singles_resp,2)]
&BEGIN_PROVIDER [integer, list_singles_rev, (n_occ_resp,n_virt_resp)]
 implicit none
 integer :: i,j,orbi,orbj,kk
 kk = 0
 do i = 1, n_occ_resp
  orbi = list_occ_resp(i)
  do j = 1, n_virt_resp
   orbj = list_virt_resp(j)
   kk += 1
   list_singles(kk,1) = orbi
   list_singles(kk,2) = orbj
   list_singles_rev(i,j) = kk
  enddo
 enddo

END_PROVIDER 

