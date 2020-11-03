
source ~/programs/qp2/quantum_package.rc

rm -rf B2.ezfio*
rm two_rdm 
rm rotation_matrix
rm new_mo_coef
rm hcore*
rm pouet 


qp create_ezfio -b 6-31g B2.xyz 
qp run scf 
#qp set_frozen_core -l
qp set_mo_class -c "[1-4]" -a "[5-18]"
qp set davidson threshold_davidson 1.e-18 
qp run cisd | tee ${EZFIO_FILE}.cisd.out
#qp run cisd | tee ${EZFIO_FILE}.cisd.out
qp set mu_of_r mu_of_r_potential cas_ful 
qp set becke_numerical_grid grid_type_sgn 0 
qp run basis_correction | tee ${EZFIO_FILE}.DFT_before.out 
qp run test_2_rdm | tee ${EZFIO_FILE}.cisd_two_e.out 
# Writing in plain text the rotation matrix to move to natural orbitals and transformed 2-RDM 
qp run write_rotation_matrix | tee ${EZFIO_FILE}.rot.out 
nlines=`wc -l two_rdm | cut -d "t" -f 1`
nlines=${nlines%%*( )} #  nlines is the number of lines in the two-rdm file
qp set yuan_plugins dim_two_rdm $nlines

# Applying the rotation matrix to the MOs and writting the TWO-RDM in the EZFIO file 
qp run new_mo_coefs_rot_mat | tee ${EZFIO_FILE}.new_mos.out 
qp reset -d 
#qp run print_hcore 
#cp -r ${EZFIO_FILE}/mo_basis/mo_coef.gz new.gz
# Testing the alpha-beta two-e energy 
qp run test_e_two_rdm | tee ${EZFIO_FILE}.test_two_e
# Testing the new MOs overlap 
# You compute the basis set correction with a frozen core cisd in the natural orbitals basis 
qp run basis_correction | tee ${EZFIO_FILE}.DFT_after.out 
