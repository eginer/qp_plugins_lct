BEGIN_PROVIDER [ logical, read_ao_integrals_mu_of_r ]
&BEGIN_PROVIDER [ logical, read_mo_integrals_mu_of_r ]
&BEGIN_PROVIDER [ logical, write_ao_integrals_mu_of_r ]
&BEGIN_PROVIDER [ logical, write_mo_integrals_mu_of_r ]

 BEGIN_DOC
! One level of abstraction for disk_ao_integrals_mu_of_r and disk_mo_integrals_mu_of_r
 END_DOC
implicit none

    if (disk_ao_integrals_mu_of_r.EQ.'Read') then
        read_ao_integrals_mu_of_r =  .True.
        write_ao_integrals_mu_of_r = .False.

    else if  (disk_ao_integrals_mu_of_r.EQ.'Write') then
        read_ao_integrals_mu_of_r = .False.
        write_ao_integrals_mu_of_r =  .True.
    
    else if (disk_ao_integrals_mu_of_r.EQ.'None') then
        read_ao_integrals_mu_of_r = .False.
        write_ao_integrals_mu_of_r = .False.

    else
        print *, 'bielec_integrals_mu_of_r/disk_ao_integrals_mu_of_r has a wrong type'
        stop 1

    endif

    if (disk_mo_integrals_mu_of_r.EQ.'Read') then
        read_mo_integrals_mu_of_r =  .True.
        write_mo_integrals_mu_of_r = .False.

    else if  (disk_mo_integrals_mu_of_r.EQ.'Write') then
        read_mo_integrals_mu_of_r = .False.
        write_mo_integrals_mu_of_r =  .True.

    else if (disk_mo_integrals_mu_of_r.EQ.'None') then
        read_mo_integrals_mu_of_r = .False.
        write_mo_integrals_mu_of_r = .False.

    else
        print *, 'bielec_integrals_mu_of_r/disk_mo_integrals_mu_of_r has a wrong type'
        stop 1

    endif

END_PROVIDER
