BEGIN_PROVIDER [ logical, read_ao_integrals_ijkl_r3 ]
&BEGIN_PROVIDER [ logical, read_mo_integrals_ijkl_r3 ]
&BEGIN_PROVIDER [ logical, write_ao_integrals_ijkl_r3 ]
&BEGIN_PROVIDER [ logical, write_mo_integrals_ijkl_r3 ]

 BEGIN_DOC
! One level of abstraction for disk_access_ao_integrals_ijkl_r3 and disk_access_mo_integrals_ijkl_r3
 END_DOC
implicit none

    if (disk_access_ao_ijkl_r3.EQ.'Read') then
        read_ao_integrals_ijkl_r3 =  .True.
        write_ao_integrals_ijkl_r3 = .False.

    else if  (disk_access_ao_ijkl_r3.EQ.'Write') then
        read_ao_integrals_ijkl_r3 = .False.
        write_ao_integrals_ijkl_r3 =  .True.
    
    else if (disk_access_ao_ijkl_r3.EQ.'None') then
        read_ao_integrals_ijkl_r3 = .False.
        write_ao_integrals_ijkl_r3 = .False.

    else
        print *, 'bielec_integrals_ijkl_r3/disk_access_ao_ijkl_r3 has a wrong type'
        stop 1

    endif

    if (disk_access_mo_ijkl_r3.EQ.'Read') then
        read_mo_integrals_ijkl_r3 =  .True.
        write_mo_integrals_ijkl_r3 = .False.

    else if  (disk_access_mo_ijkl_r3.EQ.'Write') then
        read_mo_integrals_ijkl_r3 = .False.
        write_mo_integrals_ijkl_r3 =  .True.

    else if (disk_access_mo_ijkl_r3.EQ.'None') then
        read_mo_integrals_ijkl_r3 = .False.
        write_mo_integrals_ijkl_r3 = .False.

    else
        print *, 'bielec_integrals_ijkl_r3/disk_access_mo_ijkl_r3 has a wrong type'
        stop 1

    endif

END_PROVIDER
