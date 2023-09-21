BEGIN_PROVIDER[integer, hf_index]
 BEGIN_DOC
 ! Provide the index of the HF determinant
 END_DOC

 implicit none
 integer :: i,degree

 hf_index=0

 do i=1, N_det
  call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int)
  if (degree==0) then
   hf_index=i
  endif
 enddo

 if(hf_index==0) then
  print*, "WARNING : The Hartree-Fock determinant was not found"
 endif
END_PROVIDER
