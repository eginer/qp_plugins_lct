
BEGIN_PROVIDER [integer, n_good_ao]
 implicit none
 n_good_ao = n_non_metal * 4

END_PROVIDER 

BEGIN_PROVIDER [integer, index_good_ao, (n_non_metal, 4)]
 implicit none
 BEGIN_DOC
 ! index_good_ao(i,a) = index of the good AO for non_metal atom i
 !                 a  = 1 : px
 !                 a  = 2 : py
 !                 a  = 3 : pz
 !                 a  = 4 : s 
 END_DOC
 integer :: i
 do i = 1, n_non_metal
  index_good_ao(i,1) = 138 + (i-1) * 29 + 9
  index_good_ao(i,2) = 138 + (i-1) * 29 + 10
  index_good_ao(i,3) = 138 + (i-1) * 29 + 11
  index_good_ao(i,4) = 138 + (i-1) * 29 + 3  
 enddo

END_PROVIDER 

BEGIN_PROVIDER [integer, n_ref_vec_ext]
 implicit none
 n_ref_vec_ext = n_external_ligand * 4
END_PROVIDER 


BEGIN_PROVIDER [double precision, ref_vec_ext, (ao_num,n_ref_vec_ext)]
 implicit none
 BEGIN_DOC 
! reference vector for the external ligands 
 END_DOC
 integer :: i,index_i,index_final
 integer :: j,index_j,index_final_j
 integer :: n_coord,metal_atom
 double precision :: x,y,z,d
 integer :: n_tmp
 n_tmp   = 0
 ref_vec_ext = 0.d0 
 do i = 1, n_external_ligand
  index_i     = index_external_ligand(i)    ! index of non metal atom 
  index_final = index_non_metal(index_i)    ! index of atom 
  metal_atom  = coord_non_metal(index_i,1)  ! index of connected metal atom 
  x = nucl_coord(metal_atom,1) - nucl_coord(index_final,1)
  y = nucl_coord(metal_atom,2) - nucl_coord(index_final,2)
  z = nucl_coord(metal_atom,3) - nucl_coord(index_final,3)
  d = dsqrt(x**2 + y**2 + z**2)
  x = x / d
  y = y / d
  z = z / d
  if(dabs(d).lt.1.d-10)then
   print*,'n_tmp = ',n_tmp
   print*,'x,y,z'
   print*, x,y,z 
  endif
  n_tmp += 1
  ! ref vector that points toward metal 
  ref_vec_ext(index_good_ao(index_i,1),n_tmp) = x ! x component on px 
  ref_vec_ext(index_good_ao(index_i,2),n_tmp) = y ! y component on py 
  ref_vec_ext(index_good_ao(index_i,3),n_tmp) = z ! z component on pz 
  
  double precision :: x1,y1,z1
  double precision :: x2,y2,z2
  double precision :: ovrp
  j = index_other_ext(i)                      ! other external atom 
  index_j = index_external_ligand(j)          ! index of non metal atom 
  index_final_j = index_non_metal(index_j)    ! index of atom 
  ! direction : from the other external atom toward the metalic atom 
  x1 = nucl_coord(metal_atom,1) - nucl_coord(index_final_j,1)
  y1 = nucl_coord(metal_atom,2) - nucl_coord(index_final_j,2)
  z1 = nucl_coord(metal_atom,3) - nucl_coord(index_final_j,3)
  d = dsqrt(x1**2 + y1**2 + z1**2)
! if(dabs(d).lt.1.d-10)then
!  print*,'n_tmp = ',n_tmp
!  print*,'d = ',d
!  print*,'x1,y1,z1'
!  print*, x1,y1,z1 
! endif
  x1 = x1 / d
  y1 = y1 / d
  z1 = z1 / d
  ! computes the overlap between two directions 
  ovrp = x * x1 + y * y1 + z * z1 
  x1 = x1 - ovrp * x
  y1 = y1 - ovrp * y
  z1 = z1 - ovrp * z
  d = dsqrt(x1**2 + y1**2 + z1**2)
  if(dabs(d).lt.1.d-10)then
   print*,'n_tmp = ',n_tmp
   print*,'d = ',d
   print*,'x2,y2,z2'
   print*, x1,y1,z1 
  endif
  x1 = x1 / d
  y1 = y1 / d
  z1 = z1 / d
  n_tmp += 1
  ! ref vector that points toward metal 
  ref_vec_ext(index_good_ao(index_i,1),n_tmp) = x1 ! x component on px 
  ref_vec_ext(index_good_ao(index_i,2),n_tmp) = y1 ! y component on py 
  ref_vec_ext(index_good_ao(index_i,3),n_tmp) = z1 ! z component on pz 
  
  ! ref vector orthogonal to the two first ones 
  ! vectorial product for normal vector respect to the plane of the two first  
  x2 = y * z1 - z * y1 
  y2 = z * x1 - x * z1 
  z2 = x * y1 - y * z1 
  n_tmp += 1
  ref_vec_ext(index_good_ao(index_i,1),n_tmp) = x2 ! x component on px 
  ref_vec_ext(index_good_ao(index_i,2),n_tmp) = y2 ! y component on py 
  ref_vec_ext(index_good_ao(index_i,3),n_tmp) = z2 ! z component on pz 

  n_tmp += 1
  ref_vec_ext(index_good_ao(index_i,4),n_tmp) = 1.d0 ! s component 
 enddo

END_PROVIDER 

BEGIN_PROVIDER [integer, n_ref_vec_int]
 implicit none
 n_ref_vec_int = n_internal_ligand * 4
END_PROVIDER 


BEGIN_PROVIDER [double precision, ref_vec_int, (ao_num,n_ref_vec_int)]
 implicit none
 BEGIN_DOC
 ! reference vector for the internal ligands 
 END_DOC
 integer :: i,index_i,index_final
 integer :: j,index_j,index_final_j
 integer :: n_coord,metal_atom
 double precision :: x,y,z,d
 integer :: n_tmp


 n_tmp   = 0
 ref_vec_int = 0.d0 
 do i = 1, n_internal_ligand
  index_i     = index_internal_ligand(i)    ! index of non metal atom 
  index_final = index_non_metal(index_i)    ! index of atom 
  x = metal_barycentre(1) - nucl_coord(index_final,1)
  y = metal_barycentre(2) - nucl_coord(index_final,2)
  z = metal_barycentre(3) - nucl_coord(index_final,3)
  d = dsqrt(x**2 + y**2 + z**2)
  x = x / d
  y = y / d
  z = z / d
  n_tmp += 1
  ! ref vector that points toward metal 
  ref_vec_int(index_good_ao(index_i,1),n_tmp) = x ! x component on px 
  ref_vec_int(index_good_ao(index_i,2),n_tmp) = y ! y component on py 
  ref_vec_int(index_good_ao(index_i,3),n_tmp) = z ! z component on pz 
  
  double precision :: x1,y1,z1
  double precision :: x2,y2,z2
  double precision :: ovrp
  ! direction : from the other internal atom toward the metalic atom 
  x1 = metal_difference(1) 
  y1 = metal_difference(2) 
  z1 = metal_difference(3) 
  d = dsqrt(x1**2 + y1**2 + z1**2)
  x1 = x1 / d
  y1 = y1 / d
  z1 = z1 / d
  ! computes the overlap between two directions 
  ovrp = x * x1 + y * y1 + z * z1 
  x1 = x1 - ovrp * x
  y1 = y1 - ovrp * y
  z1 = z1 - ovrp * z
  d = dsqrt(x1**2 + y1**2 + z1**2)
  x1 = x1 / d
  y1 = y1 / d
  z1 = z1 / d
  n_tmp += 1
  ! ref vector that points toward metal 
  ref_vec_int(index_good_ao(index_i,1),n_tmp) = x1 ! x component on px 
  ref_vec_int(index_good_ao(index_i,2),n_tmp) = y1 ! y component on py 
  ref_vec_int(index_good_ao(index_i,3),n_tmp) = z1 ! z component on pz 
  
  ! ref vector orthogonal to the two first ones 
  ! vectorial product for normal vector respect to the plane of the two first  
  x2 = y * z1 - z * y1 
  y2 = z * x1 - x * z1 
  z2 = x * y1 - y * z1 
  n_tmp += 1
  ref_vec_int(index_good_ao(index_i,1),n_tmp) = x2 ! x component on px 
  ref_vec_int(index_good_ao(index_i,2),n_tmp) = y2 ! y component on py 
  ref_vec_int(index_good_ao(index_i,3),n_tmp) = z2 ! z component on pz 

  n_tmp += 1
  ref_vec_int(index_good_ao(index_i,4),n_tmp) = 1.d0 ! s component 
 enddo

END_PROVIDER 

BEGIN_PROVIDER [integer, n_ref_vec_total]
 implicit none
 n_ref_vec_total = (n_ref_vec_ext + n_ref_vec_int) 
 if (n_ref_vec_total .ne. n_non_metal * 4)then
  print*,'PB with n_ref_vec'
  print*,n_ref_vec_total,n_non_metal * 4
  stop
 endif
END_PROVIDER 

BEGIN_PROVIDER [double precision, ref_vec_total, (ao_num, n_ref_vec_total)]
 implicit none
 integer :: i,j,ntmp
 ntmp = 0
 do i = 1, n_ref_vec_ext
  ntmp += 1
  do j = 1, ao_num
   ref_vec_total(j,ntmp) = ref_vec_ext(j,i)
  enddo
 enddo
 ntmp = n_ref_vec_ext

 do i = 1, n_ref_vec_int
  ntmp += 1
  do j = 1, ao_num
   ref_vec_total(j,ntmp) = ref_vec_int(j,i)
  enddo
 enddo
END_PROVIDER 
