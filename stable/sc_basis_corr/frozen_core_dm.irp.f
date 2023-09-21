
subroutine routine_save_one_e_dm_fc
 implicit none
 BEGIN_DOC
 ! routine called by :c:func:`save_one_e_dm`
 END_DOC
 double precision, allocatable  :: mo_one_e_alpha_fc(:,:)
 double precision, allocatable  ::  mo_one_e_beta_fc(:,:)
 double precision, allocatable  :: ao_one_e_alpha_fc(:,:)
 double precision, allocatable  ::  ao_one_e_beta_fc(:,:)
 double precision :: mo_alpha, mo_beta
 allocate(mo_one_e_alpha_fc(mo_num, mo_num),mo_one_e_beta_fc(mo_num, mo_num))
 allocate(ao_one_e_alpha_fc(ao_num, ao_num),ao_one_e_beta_fc(ao_num, ao_num))
 mo_one_e_alpha_fc = one_e_dm_mo_alpha(1:mo_num, 1:mo_num, 1)
 mo_one_e_beta_fc = one_e_dm_mo_beta(1:mo_num, 1:mo_num, 1)
 integer :: i,j,k,l,ii,jj
 do ii = 1, n_core
  i = list_core(ii)
  do j = 1, mo_num
   mo_one_e_alpha_fc(j,i) = 0.d0
   mo_one_e_alpha_fc(i,j) = 0.d0
   mo_one_e_beta_fc(j,i) = 0.d0
   mo_one_e_beta_fc(i,j) = 0.d0
  enddo
 enddo

 ao_one_e_alpha_fc = 0.d0
 ao_one_e_beta_fc = 0.d0
 do k = 1, ao_num
   do l = 1, ao_num
     do i = 1, mo_num
       do j = 1, mo_num
         mo_alpha = mo_one_e_alpha_fc(j,i) ! no core density alpha
         mo_beta  = mo_one_e_beta_fc(j,i)  ! no core density beta 
         !    if(dabs(dm_mo).le.1.d-10)cycle
         ao_one_e_alpha_fc(l,k) += mo_coef(k,i) * mo_coef(l,j) *  mo_alpha
         ao_one_e_beta_fc(l,k) += mo_coef(k,i) * mo_coef(l,j)  *  mo_beta
       enddo
     enddo
   enddo
 enddo


 call ezfio_set_aux_quantities_data_one_e_dm_alpha_mo(mo_one_e_alpha_fc)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_mo(mo_one_e_alpha_fc)
 call ezfio_set_aux_quantities_data_one_e_dm_alpha_ao(ao_one_e_alpha_fc)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_ao(ao_one_e_beta_fc)
end
