BEGIN_PROVIDER [integer, n_pi_z_aos]
 implicit none
 integer :: i,j,i_ao
 n_pi_z_aos = 0
 do i = 1, nucl_num
  if (dabs(nucl_charge(i)-nucl_target).lt.1.d-10) then
   do j = 1, nucl_n_aos(i)
    i_ao = nucl_aos(i,j)
    if(ao_power(i_ao, 3).eq.1)then
     n_pi_z_aos += 1
    endif
   enddo
  endif
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ integer, list_pi_z_aos, (n_pi_z_aos)]
 implicit none 
 integer :: i,j,i_ao,n_ao
 n_ao = 0
 do i = 1, nucl_num
  if (dabs(nucl_charge(i)-nucl_target).lt.1.d-10) then
   do j = 1,nucl_n_aos(i)
    i_ao = nucl_aos(i,j)
    if(ao_power(i_ao, 3).eq.1)then
     n_ao += 1
     list_pi_z_aos(n_ao) = i_ao
    endif
   enddo
  endif
 enddo

END_PROVIDER

BEGIN_PROVIDER [ integer, n_pi_z_mos ]
 implicit none
 integer :: i_ao,j_mo,i_ao_z
 double precision :: accu
 n_pi_z_mos = 0 
 do j_mo = 1, mo_num
  accu = 0.d0
  do i_ao = 1, n_pi_z_aos
   i_ao_z = list_pi_z_aos(i_ao)
   accu += dabs(mo_coef(i_ao_z,j_mo))
  enddo
  if(accu.gt.1.d-5)then
   n_pi_z_mos += 1
  endif
 enddo
END_PROVIDER

BEGIN_PROVIDER [ integer, list_pi_z_mos, (n_pi_z_mos)]
 integer :: i_ao,j_mo,i_ao_z, n_mos
 double precision :: accu
 n_mos = 0 
 do j_mo = 1, mo_num
  accu = 0.d0
  do i_ao = 1, n_pi_z_aos
   i_ao_z = list_pi_z_aos(i_ao)
   accu += dabs(mo_coef(i_ao_z,j_mo))
  enddo
  if(accu.gt.1.d-5)then
   n_mos += 1
   list_pi_z_mos(n_mos) = j_mo
  endif
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [ integer, n_occ_pi_z_mos ]
&BEGIN_PROVIDER [ integer, n_virt_pi_z_mos ]
 implicit none
 integer :: i, i_mo
 n_occ_pi_z_mos = 0
 n_virt_pi_z_mos = 0
 do i = 1, n_pi_z_mos
  i_mo = list_pi_z_mos(i)
  if(i_mo.le.elec_alpha_num)then
   n_occ_pi_z_mos += 1
  else
   n_virt_pi_z_mos += 1
  endif
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ integer, list_occ_pi_z_mos, (n_occ_pi_z_mos) ]
&BEGIN_PROVIDER [ integer, list_virt_pi_z_mos, (n_virt_pi_z_mos) ]
 implicit none
 integer :: i, i_mo, n_mo_occ, n_mo_virt
 n_mo_occ  = 0
 n_mo_virt = 0
 do i = 1, n_pi_z_mos
  i_mo = list_pi_z_mos(i)
  if(i_mo.le.elec_alpha_num)then
   n_mo_occ += 1
   list_occ_pi_z_mos(n_mo_occ)   = i_mo
  else
   n_mo_virt += 1
   list_virt_pi_z_mos(n_mo_virt) = i_mo
  endif
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ logical, is_pi_z_mos, (mo_num)]
 implicit none
 integer :: i_mo,i
 is_pi_z_mos = .False.
 do i =1, n_pi_z_mos
  i_mo = list_pi_z_mos(i)
  is_pi_z_mos(i_mo) = .True.
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ integer, n_occ_non_pi_z_mos ]
&BEGIN_PROVIDER [ integer, n_virt_non_pi_z_mos]
 implicit none
 integer :: i
 n_occ_non_pi_z_mos = 0
 n_virt_non_pi_z_mos= 0
 do i = 1, mo_num
  if(is_pi_z_mos(i))cycle
  if(i.le.elec_alpha_num)then
   n_occ_non_pi_z_mos += 1
  else 
   n_virt_non_pi_z_mos+= 1
  endif
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [ integer, list_occ_non_pi_z_mos, (n_occ_non_pi_z_mos) ]
&BEGIN_PROVIDER [ integer, list_virt_non_pi_z_mos, (n_virt_non_pi_z_mos) ]
 implicit none
 integer :: i, i_mo, n_occ, n_virt
 n_occ = 0
 n_virt= 0
 do i = 1, mo_num
  if(is_pi_z_mos(i))cycle
  if(i.le.elec_alpha_num)then
   n_occ  += 1
   list_occ_non_pi_z_mos(n_occ) = i
  else 
   n_virt += 1
   list_virt_non_pi_z_mos(n_virt) = i
  endif
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ integer, new_order_mos, (mo_num)]
 implicit none
 integer :: i,i_mo,n_mo
 n_mo = 0
 ! first the non pi z mos occupied
 do i = 1, n_occ_non_pi_z_mos
  n_mo += 1
  i_mo = list_occ_non_pi_z_mos(i) 
  new_order_mos(n_mo) = i_mo
 enddo
 ! then the occupied pi z mos 
 do i = 1, n_occ_pi_z_mos
  n_mo += 1
  i_mo = list_occ_pi_z_mos(i)  
  new_order_mos(n_mo) = i_mo
 enddo
 ! then the virt pi z mos 
 do i = 1, n_occ_pi_z_mos
  n_mo += 1
  i_mo = list_virt_pi_z_mos(i)
  new_order_mos(n_mo) = i_mo
 enddo
 ! then the virt non pi z mos 
 do i = 1, n_occ_non_pi_z_mos
  n_mo += 1
  i_mo = list_virt_non_pi_z_mos(i)
  new_order_mos(n_mo) = i_mo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, new_mo_coef_pi_z_order, (ao_num, mo_num)]
 implicit none
 integer :: i,j, i_mo
 do i = 1, mo_num
  i_mo = new_order_mos(i)
  do j = 1, ao_num
   new_mo_coef_pi_z_order(j,i) = mo_coef(j,i_mo)
  enddo
 enddo

END_PROVIDER 

subroutine save_mos_pi_z
  implicit none
  double precision, allocatable  :: buffer(:,:)
  integer                        :: i,j

  call ezfio_set_mo_basis_mo_num(mo_num)                                                                                                                
  call ezfio_set_mo_basis_mo_label(mo_label)
  call ezfio_set_mo_basis_ao_md5(ao_md5)
  allocate ( buffer(ao_num,mo_num) )
  buffer = 0.d0
  do j = 1, mo_num
    do i = 1, ao_num
      buffer(i,j) = new_mo_coef_pi_z_order(i,j)
    enddo 
  enddo
  call ezfio_set_mo_basis_mo_coef(buffer)
  call ezfio_set_mo_basis_mo_occ(mo_occ)
  call ezfio_set_mo_basis_mo_class(mo_class)
  deallocate (buffer)
  
end

