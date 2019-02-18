BEGIN_PROVIDER [integer, n_metal]
 implicit none
 BEGIN_DOC
! number of atoms with Z > 18
 END_DOC
 integer :: i
 n_metal = 0
 do i = 1, nucl_num
  if(dabs(nucl_charge(i)).gt.18.d0) then 
   n_metal += 1  
  endif
 enddo

END_PROVIDER 


BEGIN_PROVIDER [integer, index_metal, (n_metal)]
 implicit none
 BEGIN_DOC 
! index of atoms with Z > 18
 END_DOC
 index_metal = 0
 integer :: i,jtmp
 jtmp = 0
 do i = 1, nucl_num
  if(dabs(nucl_charge(i)).gt.18.d0) then 
   jtmp += 1  
   index_metal(jtmp) = i
  endif
 enddo

END_PROVIDER 


BEGIN_PROVIDER [integer, n_non_metal]
 implicit none
 BEGIN_DOC
! number of atoms with Z <= 18
 END_DOC
 integer :: i
 n_non_metal = 0
 do i = 1, nucl_num
  if(dabs(nucl_charge(i)).le.18.d0) then 
   n_non_metal += 1  
  endif
 enddo

END_PROVIDER 


BEGIN_PROVIDER [integer, index_non_metal, (n_non_metal)]
 implicit none
 BEGIN_DOC 
! index of atoms with Z <= 18
 END_DOC
 index_non_metal = 0
 integer :: i,jtmp
 jtmp = 0
 do i = 1, nucl_num
  if(dabs(nucl_charge(i)).le.18.d0) then 
   jtmp += 1  
   index_non_metal(jtmp) = i
  endif
 enddo

END_PROVIDER 



 BEGIN_PROVIDER [integer, n_external_ligand]
&BEGIN_PROVIDER [integer, n_internal_ligand]
 implicit none
 BEGIN_DOC
! n_external_ligand = number of non_metal_atoms with coordination = 1
! n_internal_ligand = number of non_metal_atoms with coordination = 2
 END_DOC
 integer :: i
 n_external_ligand = 0
 n_internal_ligand = 0
 do i = 1, n_non_metal
  if (n_coord_non_metal(i) .eq. 1)then 
   n_external_ligand += 1
  else 
   n_internal_ligand += 1
  endif
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [integer, index_external_ligand, (n_external_ligand)]
&BEGIN_PROVIDER [integer, index_internal_ligand, (n_internal_ligand)]
 implicit none
 integer :: i,int_tmp,ext_tmp
 int_tmp = 0
 ext_tmp = 0
 do i = 1, n_non_metal
  if (n_coord_non_metal(i) .eq. 1)then 
   ext_tmp += 1
   index_external_ligand(ext_tmp) = i
  else 
   int_tmp += 1
   index_internal_ligand(int_tmp) = i
  endif
 enddo

END_PROVIDER 


