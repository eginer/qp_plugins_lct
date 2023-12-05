program get_pi_cas
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
 double precision, allocatable :: mo_coef_new(:,:)
 allocate(mo_coef_new(ao_num, mo_num))
 integer :: i,i_mo
! print*,'n_pi_z_aos = ',n_pi_z_aos
! do i = 1, n_pi_z_aos
!  print*,list_pi_z_aos(i)
! enddo
 print*,'n_pi_z_mos',n_pi_z_mos
 do i = 1, n_pi_z_mos
  i_mo = list_pi_z_mos(i)
  print*,i_mo
 enddo
 print*,'new order'
 do i = 1, mo_num
  print*,i, new_order_mos(i)
 enddo
! do i = 1, n_occ_pi_z_mos
!  print*, list_occ_pi_z_mos(i)
! enddo
 call save_mos_pi_z
end
