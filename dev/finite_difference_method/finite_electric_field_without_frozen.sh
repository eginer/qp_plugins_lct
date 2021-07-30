#/bin/qpsh
# specify the QP folder
QP=$QP_ROOT
#sourcing the quantum_package.rc file
. ${QP}/quantum_package.rc


thresh=0.000000000001
PT2_1=0.01
PT2_2=0.001
PT2=0.000001
ezfio=$1
basis=$2
xyz=$3
plus_part=$4
field_abs=0.0001

#if [[ -z $field_abs ]]; then
# echo "Field_abs will be set to 0.0001 a.u. by default."
# field_abs=0.0001
#fi

if [[ -z $plus_part ]]; then
 echo "We perform calc for the positive and negative values for the field strenght."
 echo "Absolute field strenght :  "$field_abs
 echo "SCF thresh :" $thresh
 echo "CIPSI PT2 : " $PT2  

 func=pbe_ueg

# pbe-ueg no_core_density a true frozen_core  10-12 thresh scf

 field=$field_abs
 mkdir plus
 cp ${xyz} plus/
 cd plus
 
 #create ezfio 
 qp create_ezfio -b ${basis} -o ${ezfio}_plus ${xyz} 
 qp set_file ${ezfio}_plus 

 #set ezfio
 qp set finite_difference_method field_strenght  $field 
 qp set scf_utils thresh_scf $thresh 
 qp run write_eff_basis_ints_electric_field

 #run scf and dipole calculations
 qp run scf | tee ${ezfio}.${field}.out
# qp_set_frozen_core ${ezfio}_plus
# qp set density_for_dft no_core_density true

# qp run density_at_r | tee HF_n_n2_r_plus.dat 
# qp run print_z_dipole | tee HF_z_dipole_plus.out

 #basis set correction w n_HF
# qp set sc_basis_corr basis_cor_func $func
# qp_set_frozen_core ${ezfio}_plus
# qp set density_for_dft no_core_density true
 qp run basis_correction | tee ${ezfio}.${field}.BSC_n_HF.out 

# #CIPSI energy calculations
# qp set perturbation pt2_max ${PT2_1}
# qp run fci | tee fci_pt2_${PT2_1}_${basis}_plus.out
# qp run save_natorb

# qp set perturbation pt2_max ${PT2_2}
# qp run fci | tee fci_pt2_${PT2_2}_${basis}_plus.out
# qp run save_natorb

# qp set perturbation pt2_max ${PT2}
# qp run fci | tee fci_pt2_${PT2}_${basis}_plus.out

# qp run density_at_r | tee CIPSI_n_n2_r_plus.dat 

##basis set correction w n_CIPSI
#qp set sc_basis_corr basis_cor_func $func
#qp_set_frozen_core ${ezfio}_plus
#qp set density_for_dft no_core_density true
#qp run basis_correction | tee ${ezfio}.${field}.BSC_n_CIPSI.out 

 qp unset_file
 cd ..

else
 echo "We perform calc only for negative value for the field strenght."
 echo "Absolute field strenght :  "$field_abs
 echo "SCF thresh :" $thresh
 echo "CIPSI PT2 : " $PT2  

 func=pbe_ueg
fi

 
 field=-$field_abs
 mkdir minus
 cp ${xyz} minus
 cd minus

 #create ezfio 
 qp create_ezfio -b ${basis} -o ${ezfio}_minus ${xyz} 
 qp set_file ${ezfio}_minus 

 #set ezfio
 qp set finite_difference_method field_strenght  $field 
 qp set scf_utils thresh_scf $thresh 
 qp run write_eff_basis_ints_electric_field

 #run scf and dipole calculations
 qp run scf | tee ${ezfio}.${field}.out
# qp_set_frozen_core ${ezfio}_minus
# qp set density_for_dft no_core_density true

# qp run density_at_r | tee HF_n_n2_r_minus.dat 
# qp run print_z_dipole | tee HF_z_dipole_minus.out

 #basis set correction w n_HF
# qp set sc_basis_corr basis_cor_func $func
# qp_set_frozen_core ${ezfio}_minus
# qp set density_for_dft no_core_density true
 qp run basis_correction | tee ${ezfio}.${field}.BSC_n_HF.out 

# #CIPSI energy calculations
# qp set perturbation pt2_max ${PT2_1}
# qp run fci
# qp run save_natorb

# qp set perturbation pt2_max ${PT2_2}
# qp run fci
# qp run save_natorb

# qp set perturbation pt2_max ${PT2}
# qp run fci | tee fci_pt2_${PT2}_${basis}_minus.out

# qp run density_at_r | tee CIPSI_n_n2_r_minus.dat 

# #basis set correction w n_CIPSI
# qp set sc_basis_corr basis_cor_func $func
# qp_set_frozen_core ${ezfio}_minus
# qp set density_for_dft no_core_density true
# qp run basis_correction | tee ${ezfio}.${field}.BSC_n_CIPSI.out 

 qp unset_file
cd ..
