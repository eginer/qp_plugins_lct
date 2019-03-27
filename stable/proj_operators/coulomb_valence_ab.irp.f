 BEGIN_PROVIDER [integer, n_valence_orb_for_hf,(2)]
&BEGIN_PROVIDER [integer, n_max_valence_orb_for_hf]
 implicit none
 integer :: i
 n_valence_orb_for_hf = 0
 do i = 1, elec_alpha_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   n_valence_orb_for_hf(1) +=1 
  endif
 enddo

 do i = 1, elec_beta_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   n_valence_orb_for_hf(2) +=1 
  endif
 enddo
 n_max_valence_orb_for_hf = maxval(n_valence_orb_for_hf)

END_PROVIDER 

BEGIN_PROVIDER [integer, list_valence_orb_for_hf, (n_max_valence_orb_for_hf,2)]
 implicit none
 integer :: i,j
 j = 0
 do i = 1, elec_alpha_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   j +=1 
   list_valence_orb_for_hf(j,1) = i
  endif
 enddo

 j = 0
 do i = 1, elec_beta_num
  if(  trim(mo_class(i))=="Inactive" & 
  .or. trim(mo_class(i))=="Active"   &
  .or. trim(mo_class(i))=="Virtual" )then
   j +=1 
   list_valence_orb_for_hf(j,2) = i
  endif
 enddo

END_PROVIDER 

BEGIN_PROVIDER [integer, n_inact_act_virt_orb]
 implicit none
 integer :: i
 n_inact_act_virt_orb = 0
!do i = 1, mo_num
! if(  trim(mo_class(i))=="Inactive" & 
! .or. trim(mo_class(i))=="Active"   &
! .or. trim(mo_class(i))=="Virtual" )then
!  n_inact_act_virt_orb +=1 
! endif
!enddo
 n_inact_act_virt_orb = mo_num 
END_PROVIDER 

BEGIN_PROVIDER [integer, list_inact_act_virt_orb, (n_inact_act_virt_orb)]
 implicit none
 integer :: i,j
!j = 0
!do i = 1, mo_num
! if(  trim(mo_class(i))=="Inactive" & 
! .or. trim(mo_class(i))=="Active"   &
! .or. trim(mo_class(i))=="Virtual" )then
!  j +=1 
!  list_inact_act_virt_orb (j) = i
! endif
!enddo

 do i = 1, mo_num
   list_inact_act_virt_orb(i) = i
 enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, integrals_for_valence_hf_pot, (n_inact_act_virt_orb,n_inact_act_virt_orb,n_max_valence_orb_for_hf,n_max_valence_orb_for_hf)]
 implicit none
 integer :: i_i,i_j,i,j,i_m,i_n,m,n
 double precision :: get_two_e_integral
 do i_m = 1, n_max_valence_orb_for_hf! electron 1 
  m = list_valence_orb_for_hf(i_m,1)
  do i_n = 1, n_max_valence_orb_for_hf! electron 2 
   n = list_valence_orb_for_hf(i_n,1)
   do i_i = 1, n_inact_act_virt_orb ! electron 1 
    i = list_inact_act_virt_orb(i_i)
    do i_j = 1, n_inact_act_virt_orb ! electron 2 
     j = list_inact_act_virt_orb(i_j)
     !                             2   1   2   1
     integrals_for_valence_hf_pot(i_j,i_i,i_n,i_m) = get_two_e_integral(m,n,i,j,mo_integrals_map) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

subroutine f_HF_valence_ab(r1,r2,integral_psi,two_bod)
 implicit none
 BEGIN_DOC
! f_HF_ab(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha beta spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > =  \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: integral_psi,two_bod
 integer :: i,j,m,n,i_m,i_n
 integer :: i_i,i_j
 double precision :: mo_two_e_integral
 double precision :: mos_array_r1(mo_num)
 double precision :: mos_array_r2(mo_num)
 double precision, allocatable  :: mos_array_valence_r1(:)
 double precision, allocatable  :: mos_array_valence_r2(:)
 double precision, allocatable  :: mos_array_valence_hf_r1(:)
 double precision, allocatable  :: mos_array_valence_hf_r2(:)
 double precision :: get_two_e_integral
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 allocate(mos_array_valence_r1( n_inact_act_virt_orb ) , mos_array_valence_r2( n_inact_act_virt_orb ) )
 allocate(mos_array_valence_hf_r1( n_valence_orb_for_hf(1) ) , mos_array_valence_hf_r2( n_valence_orb_for_hf(2) ) )
 do i_m = 1, n_valence_orb_for_hf(1)
  mos_array_valence_hf_r1(i_m) = mos_array_r1(list_valence_orb_for_hf(i_m,1))
 enddo
 do i_m = 1, n_valence_orb_for_hf(2)
  mos_array_valence_hf_r2(i_m) = mos_array_r2(list_valence_orb_for_hf(i_m,2))
 enddo

 do i_m = 1, n_inact_act_virt_orb 
  mos_array_valence_r1(i_m) = mos_array_r1(list_inact_act_virt_orb(i_m))
  mos_array_valence_r2(i_m) = mos_array_r2(list_inact_act_virt_orb(i_m))
 enddo


 integral_psi = 0.d0
 two_bod = 0.d0
 do m = 1, n_valence_orb_for_hf(1)! electron 1 
  do n = 1, n_valence_orb_for_hf(2)! electron 2 
   two_bod += mos_array_valence_hf_r1(m) * mos_array_valence_hf_r1(m) * mos_array_valence_hf_r2(n) * mos_array_valence_hf_r2(n) 
   do i = 1, n_inact_act_virt_orb
    do j = 1, n_inact_act_virt_orb
     !                                             2 1 2 1
     integral_psi +=  integrals_for_valence_hf_pot(j,i,n,m)  & 
     * mos_array_valence_r1(i) * mos_array_valence_hf_r1(m)  & 
     * mos_array_valence_r2(j) * mos_array_valence_hf_r2(n)    
    enddo
   enddo
  enddo
 enddo
end


subroutine integral_f_HF_valence_ab(r1,integral_psi)
 implicit none
 BEGIN_DOC
! f_HF_ab(X1,X2) = function f_{\Psi^B}(X_1,X_2) of equation (22) of paper J. Chem. Phys. 149, 194301 (2018)
! for alpha beta spins and an HF wave function
! < HF | wee_{\alpha\alpha} | HF > =  \int (X1,X2) f_HF_aa(X1,X2)
 END_DOC
 double precision, intent(in) :: r1(3)
 double precision, intent(out):: integral_psi
 integer :: i,j,m,n,i_m,i_n
 integer :: i_i,i_j
 double precision :: mo_two_e_integral
 double precision :: mos_array_r1(mo_num)
 double precision, allocatable  :: mos_array_valence_r1(:)
 double precision, allocatable  :: mos_array_valence_hf_r1(:)
 double precision :: get_two_e_integral
 call give_all_mos_at_r(r1,mos_array_r1) 
 allocate(mos_array_valence_r1( n_inact_act_virt_orb ))
 allocate(mos_array_valence_hf_r1( n_valence_orb_for_hf(1) ) )
 do i_m = 1, n_valence_orb_for_hf(1)
  mos_array_valence_hf_r1(i_m) = mos_array_r1(list_valence_orb_for_hf(i_m,1))
 enddo
 do i_m = 1, n_inact_act_virt_orb 
  mos_array_valence_r1(i_m) = mos_array_r1(list_inact_act_virt_orb(i_m))
 enddo

 integral_psi = 0.d0
 do m = 1, n_valence_orb_for_hf(1)! electron 1 
  do n = 1, n_valence_orb_for_hf(2)! electron 2 
   do i = 1, n_inact_act_virt_orb
    do j = n,n
     !                                             2 1 2 1
     integral_psi +=  integrals_for_valence_hf_pot(j,i,n,m)  & 
     * mos_array_valence_r1(i) * mos_array_valence_hf_r1(m) 
    enddo
   enddo
  enddo
 enddo
end
