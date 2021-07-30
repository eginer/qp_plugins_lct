========================
finite_difference_method
========================
Here, we considere the Hamiltonian : H = T + V_ne + W_ee + \vec{mu}.\vec{E}. 


05/02/20 : We considere \vec{mu}.\vec{E} = - eps.<mu.z> as we make tests on a system were the two atoms are on the axis z. Epsilon will be taken as 0.0001 as proposed in HalKloHel-JCP-99.
To do later : a generalization to any orientation of the system.

Structure of the code (everything here is a copy past of src files, if any need of reference) :

ao_one_e_ints_electric_field.irp.f contains the modified one electron integrals provders in ao basis.
mo_one_e_ints_electric_field.irp.f contains the modified one electron integrals provders in mo basis. (certainly useful to build the guess in the SCF procedure)

rout_write_int_electric_field.irp.f to write the modified one electron integrals in the EZFIO folder.

The SCF calculation with the modified integrals is done in 'scf_electric_field_script'. 
Problems : 
Attribute BH.ezfio/mo_basis/mo_num is not set 
-> Maybe delete the MO writting in a first try. Then think about this mo_integrals meaning (it looks like they are used to build the guess)
