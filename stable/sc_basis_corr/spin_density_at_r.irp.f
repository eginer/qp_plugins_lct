program spin_density_at_r 
 implicit none
 BEGIN_DOC
 !
 END_DOC

 read_wf = .true.
 touch read_wf
 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 double precision :: r(3), dm_a(N_states), dm_b(N_states)
 double precision :: spin_density_at_nucleus
 integer :: istate
 double precision :: mu_bohr, mu_nuclear ! Bohr and nuclear magneton
 double precision :: g_e, g_k ! free electron g factor and nucleus g factor 
 double precision :: A_factor ! A = (4Pi/3)*g_e*g_k*mu_bohr*mu_nuclear*<s_z>*spin_density_at_nucleus = (8Pi/3)*A_Factor*spin_density_at_nucleus
 double precision :: pi, HFCC

 istate = 1
 pi = dacos(-1.d0)

 mu_bohr = 0.5d0 ! atomic unit = e*h_bar/(m_e)
 mu_nuclear = 7.622d0 ! atomic unit , MHz/T 
 g_e = 2.d0 
 g_k = 0.d0 ! 0.566 for Nitrogen from KohMil-AtomicDataAndNuclearDataTables-85
 A_factor = 54.d0 ! for Nitrogen from KohMil-AtomicDataAndNuclearDataTables-85

 
 r(1) = nucl_coord(1, 1) ! # of the nucleus in the .xyz file
 r(2) = nucl_coord(1, 2)
 r(3) = nucl_coord(1, 3)

 call dm_dft_alpha_beta_at_r(r, dm_a, dm_b) 
 
 spin_density_at_nucleus = dm_a(istate) - dm_b(istate)
 HFCC = (4.d0*pi/3.d0)*A_factor*spin_density_at_nucleus

 write(35,*) 'r(1)                                         =', r(1)
 write(35,*) 'r(2)                                         =', r(2)
 write(35,*) 'r(3)                                         =', r(3)
 write(35,*) 'spin_density_at_nucleus                      =', spin_density_at_nucleus
 write(35,*) 'dm_a                                         =', dm_a(istate)
 write(35,*) 'dm_b                                         =', dm_b(istate)
 write(35,*) ''
 write(35,*) 'hyperfine coupling constant (isotropic part) =', HFCC
end program


