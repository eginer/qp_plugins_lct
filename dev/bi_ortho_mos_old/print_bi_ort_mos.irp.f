program print_mos_lr_in_r
 implicit none
 double precision :: r(3),xmin,xmax,dx
 double precision, allocatable :: mos_r_array(:),mos_l_array(:)
 integer :: nx,i,imo,kmo
 allocate(mos_r_array(mo_num),mos_l_array(mo_num))
 xmin = -5.d0
 xmax =  5.d0
 nx = 1000
 dx = (xmax-xmin)/dble(nx)
 r = 0.d0
 r(1) = xmin
 imo = 3
 kmo = 3
 do i = 1, nx
  call give_all_mos_l_at_r(r, mos_l_array)
  call give_all_mos_r_at_r(r, mos_r_array)
  write(33,'(100(F16.10,X))')r(1),mos_l_array(kmo),mos_r_array(imo)
  r(1) += dx
 enddo
 print*,''
 do i = 1, mo_num
  print*,i,overlap_mo_l(i,i),overlap_mo_r(i,i)
 enddo
 print*,'left/left'
 do i = 1, mo_num
  print*,i,overlap_mo_l(:,i)
 enddo
 print*,'right/right'
 do i = 1, mo_num
  print*,i,overlap_mo_r(:,i)
 enddo


end
