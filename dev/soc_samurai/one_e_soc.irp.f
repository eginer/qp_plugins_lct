BEGIN_PROVIDER [ double precision, constant_soc]
 implicit none
 double precision :: C_soc,alpha
 alpha = 1.d0/137.d0
 constant_soc = alpha**2 / 2.d0
END_PROVIDER 

BEGIN_PROVIDER [ double precision, vec_prod_grad_ao_at_r, (3,n_points_final_grid, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! vec_prod_grad_ao_at_r(coord,ipoint,i,j) =  [grad_r AO_i(r)] X [grad_r AO_j(r)]_coord
!
! where [u] X [v] is the vector product and coord = 1->, coord = 2->y, coord = 3->z.
 END_DOC
 integer :: i,j,ipoint
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1,n_points_final_grid 
    vec_prod_grad_ao_at_r(1,ipoint,j,i) = aos_grad_in_r_array_transp_3(2,ipoint,j)*aos_grad_in_r_array_transp_3(3,ipoint,i) &
                                        - aos_grad_in_r_array_transp_3(3,ipoint,j)*aos_grad_in_r_array_transp_3(2,ipoint,i)
    vec_prod_grad_ao_at_r(2,ipoint,j,i) = aos_grad_in_r_array_transp_3(3,ipoint,j)*aos_grad_in_r_array_transp_3(1,ipoint,i) &
                                        - aos_grad_in_r_array_transp_3(1,ipoint,j)*aos_grad_in_r_array_transp_3(3,ipoint,i)
    vec_prod_grad_ao_at_r(3,ipoint,j,i) = aos_grad_in_r_array_transp_3(1,ipoint,j)*aos_grad_in_r_array_transp_3(2,ipoint,i) &
                                        - aos_grad_in_r_array_transp_3(2,ipoint,j)*aos_grad_in_r_array_transp_3(1,ipoint,i)
   enddo
  enddo
 enddo
vec_prod_grad_ao_at_r = vec_prod_grad_ao_at_r 
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_one_e_soc_cartesian_per_atom, (3,nucl_num,ao_num, ao_num)]
 implicit none
 BEGIN_DOC
 ! ao_one_e_soc_cartesian_per_atom(coord,i_nucl,j,i) = \alpha^2/2 \int dr Z_nucl /|r - R_nucl| [grad_r AO_i(r)] X [grad_r AO_j(r)]_coord
 END_DOC
 integer :: i,j,ipoint,istate,k,l,coord,i_nucl
 double precision :: r_n(3),r_e(3),dis,v_nucl,z_n
 ao_one_e_soc_cartesian_per_atom = 0.d0
 do i_nucl= 1, nucl_num
  r_n(1) = nucl_coord(i_nucl,1)
  r_n(2) = nucl_coord(i_nucl,2)
  r_n(3) = nucl_coord(i_nucl,3)
  z_n = nucl_charge(i_nucl)
  do i= 1, ao_num
   do j = 1, ao_num
    do ipoint = 1, n_points_final_grid
     r_e(1) = final_grid_points (1,ipoint)
     r_e(2) = final_grid_points (2,ipoint)
     r_e(3) = final_grid_points (3,ipoint)
     dis = sqrt((r_e(1)-r_n(1))**2 + (r_e(2)-r_n(2))**2 + (r_e(3)-r_n(3))**2 )
     if (dis < 1.d-8) cycle
     v_nucl = constant_soc * z_n / dis
     do coord = 1,3
      ao_one_e_soc_cartesian_per_atom(coord,i_nucl,j, i) += vec_prod_grad_ao_at_r(coord,ipoint,j,i) * v_nucl * final_weight_at_r_vector(ipoint)
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_one_e_soc_cartesian, (ao_num, ao_num,3)]
 BEGIN_DOC
 ! ao_one_e_soc_cartesian_per_atom(j,i,coord) = \alpha^2/2 \sum_{nucl} \int dr Z_nucl /|r - R_nucl| [grad_r AO_i(r)] X [grad_r AO_j(r)]_coord
 END_DOC
 implicit none 
 integer :: i,j,ipoint,istate,k,l,coord,i_nucl
 ao_one_e_soc_cartesian = 0.d0
 do coord = 1, 3
  do i = 1, ao_num
   do j = 1, ao_num
    do i_nucl = 1, nucl_num
     ao_one_e_soc_cartesian(j,i,coord) += ao_one_e_soc_cartesian_per_atom(coord,i_nucl, j, i)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_one_e_soc_cartesian, (mo_num, mo_num,3)]
 implicit none
 BEGIN_DOC
 ! mo_one_e_soc_cartesian_per_atom(j,i,coord) = \alpha^2/2 \sum_{nucl} \int dr Z_nucl /|r - R_nucl| [grad_r MO_i(r)] X [grad_r MO_j(r)]_coord
 END_DOC
 integer :: i
 mo_one_e_soc_cartesian = 0.d0
 do i = 1, 3
  call ao_to_mo(ao_one_e_soc_cartesian(1,1,i),ao_num,mo_one_e_soc_cartesian(1,1,i),mo_num)
 enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*8, mo_one_e_soc, (mo_num, mo_num, 3)]
 implicit none
 BEGIN_DOC
! mo_one_e_soc(k,l,mu) = <MO_k | V_soc_1_e ^\mu | MO_l>, 
! 
! with mu = 1 ::> V^- = V^x + i V^y, 
!      
!      mu = 2 ::> V^+ = V^x - i V^y,
!
!      mu = 3 ::> V^z 
 END_DOC
 integer :: i,j
 do i = 1, mo_num
  do j = 1, mo_num
   ! V^+ = V^x + i V^y = i(mo_one_e_soc_cartesian(1)) + i (i mo_one_e_soc_cartesian(2))
   !     = -mo_one_e_soc_cartesian(2) + i mo_one_e_soc_cartesian(1)
   mo_one_e_soc(j,i,1) = cmplx(-1.D0 * mo_one_e_soc_cartesian(j,i,2), mo_one_e_soc_cartesian(j,i,1))
   ! V^- = V^x - i V^y = i(mo_one_e_soc_cartesian(1)) - i (i mo_one_e_soc_cartesian(2))
   !     = +mo_one_e_soc_cartesian(2) + i mo_one_e_soc_cartesian(1)
   mo_one_e_soc(j,i,2) = cmplx( 1.D0 * mo_one_e_soc_cartesian(j,i,2), mo_one_e_soc_cartesian(j,i,1))
   ! V^z = i mo_one_e_soc_cartesian(3)
   mo_one_e_soc(j,i,3) = cmplx(0.d0,mo_one_e_soc_cartesian(j,i,3))
  enddo
 enddo

END_PROVIDER 
