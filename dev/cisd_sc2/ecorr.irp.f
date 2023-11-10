BEGIN_PROVIDER [ double precision, e_cor_a, (n_virt_cisd,2)]
 implicit none
 integer :: k,l,c,s1,sigma,a
 BEGIN_DOC
!e_cor_a(a,sigma) = - sum_(klc)sum_(sigma_1) e_{k,sigma l,tau}^{a,sigma c,tau}
 END_DOC
 e_cor_a = 0.d0
 do s1 = 1, 2
  do sigma = 1, 2
   do k = 1, n_inact_cisd
    do l = 1, n_inact_cisd
     do a = 1, n_virt_cisd 
      do c = 1, n_virt_cisd
       e_cor_a(a,sigma) += - ecorr_contrib(c,a,k,l,sigma,s1)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, e_cor_ab, (n_virt_cisd,n_virt_cisd,2)]
 implicit none
 integer :: k,l,b,a,tau,sigma
 BEGIN_DOC
!e_cor_ab(a,b,sigma,tau) =  sum_(kl) e_{k,sigma l,tau}^{a,sigma b,tau}
 END_DOC
 e_cor_ab = 0.d0
 do tau = 1, 2
  do sigma = 1, 2
   do l = 1, n_inact_cisd
    do k = 1, n_inact_cisd
     do b = 1, n_virt_cisd
      do a = 1, n_virt_cisd 
       e_cor_ab(a,b,sigma,tau) += ecorr_contrib(a,b,k,l,sigma,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [ double precision, e_cor_j, (n_inact_cisd,2)]
 implicit none
 integer :: k,c,d,tau,s1,j,a
 BEGIN_DOC
!e_cor_j(j,tau) = - sum_(kcd) e_{k,sigma j,tau}^{c,sigma d,tau}
 END_DOC
 e_cor_j = 0.d0
 do tau = 1, 2
  do s1= 1, 2
   do j = 1, n_inact_cisd
    do k = 1, n_inact_cisd
     do c = 1, n_virt_cisd
      do d = 1, n_virt_cisd 
       e_cor_j(j,tau) += - ecorr_contrib(d,c,j,k,s1,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, e_cor_jb_direct, (n_virt_cisd,n_inact_cisd,2)]
 implicit none
 integer :: k,c,b,j,tau,s1
 BEGIN_DOC
!e_cor_jb(j,tau) =  sum_(kc) e_{k,sigma j,tau}^{c,sigma b,tau}
 END_DOC
 e_cor_jb_direct = 0.d0
 do tau = 1, 2
  do s1= 1, 2
   do j = 1, n_inact_cisd
    do k = 1, n_inact_cisd
     do c = 1, n_virt_cisd
      do b = 1, n_virt_cisd 
       e_cor_jb_direct(b,j,tau) +=  ecorr_contrib(b,c,j,k,s1,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, e_cor_jb_exchange, (n_virt_cisd,n_inact_cisd,2,2)]
 implicit none
 integer :: k,d,a,j,tau,sigma
 BEGIN_DOC
!e_cor_jb(j,tau,sigma) =  sum_(kd) e_{k,sigma j,tau}^{a,sigma d,tau}
 END_DOC
 e_cor_jb_exchange = 0.d0
 do tau = 1, 2
  do sigma= 1, 2
   do k = 1, n_inact_cisd
    do j = 1, n_inact_cisd
     do a = 1, n_virt_cisd
      do d = 1, n_virt_cisd 
       e_cor_jb_exchange(a,j,sigma,tau) +=  ecorr_contrib(d,a,j,k,sigma,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 




BEGIN_PROVIDER [ double precision, e_cor_jab, (n_virt_cisd,n_virt_cisd,n_inact_cisd,2,2)]
 implicit none
 integer :: j,a,b,k,tau,sigma
 BEGIN_DOC
!e_cor_jab(a,b,j,tau,sigma) = - sum_(k) e_{k,sigma j,tau}^{a,sigma b,tau}
 END_DOC
 e_cor_jab = 0.d0
 do tau = 1, 2
  do sigma= 1, 2
   do k = 1, n_inact_cisd
    do j = 1, n_inact_cisd
     do a = 1, n_virt_cisd
      do b = 1, n_virt_cisd 
       e_cor_jab(b,a,j,sigma,tau) +=  - ecorr_contrib(b,a,j,k,sigma,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


