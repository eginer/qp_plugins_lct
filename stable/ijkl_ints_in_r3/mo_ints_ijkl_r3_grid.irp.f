double precision function integral_bourrin_mo(i,j,k,l)
 implicit none
 integer,intent(in) :: i,j,k,l
 integer :: m,n,q,p
 double precision :: get_ao_bielec_integral_ijkl_r3
 integral_bourrin_mo= 0.d0
 do m = 1, ao_num
  do n = 1, ao_num
   do p = 1, ao_num
    do q = 1, ao_num
     integral_bourrin_mo += mo_coef(m,i) * mo_coef(n,j) *mo_coef(p,k) * mo_coef(q,l) * get_ao_bielec_integral_ijkl_r3(q,p,n,m,ao_integrals_ijkl_r3_map) 
 !   print*,'coef mo  =  ',mo_coef(m,i), mo_coef(n,j), mo_coef(p,k), mo_coef(q,l)
 !   print*,'inte ao map  =', get_ao_bielec_integral_ijkl_r3(q,p,n,m,ao_integrals_ijkl_r3_map)
 !   print*,'inte bourrini  =  ' ,integral_bourrin_mo
    !stop
    enddo
   enddo
  enddo
 enddo


end
