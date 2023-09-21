source /home/eginer/qp2/quantum_package.rc
basis=aug-cc-pvtz
sys=water

qp create_ezfio -b $basis ${sys}.xyz -a -o ${sys}_$basis
qp run scf | tee ${EZFIO_FILE}.scf.out
qp set_frozen_core 

#### FIRST SANITY CHECK FOR CCSD CALCULATIONS WITH BARE HAMILTONIAN
qp run ccsd | tee ${EZFIO_FILE}.ccsd_qp.out
qp run fcidump_pyscf | tee ${EZFIO_FILE}.fcidump_pyscf.out 
fcidump=${EZFIO_FILE}.FCIDUMP
python PYSCF_EOMCC_from_fcidump.py $basis ${fcidump} | tee ${EZFIO_FILE}.PySCF_EOM_bare_H.out 


#### PBE-UEG CALCULATIONS 
qp set sc_basis_corr basis_cor_func pbe_ueg
qp run basis_correction | tee ${EZFIO_FILE}.basiscorr.out 
qp run fcidump_for_vbarb | tee ${EZFIO_FILE}.fcidump_for_vbarb.out
fcidump=${EZFIO_FILE}.vbarb_FCIDUMP
python PYSCF_EOMCC_from_fcidump.py $basis ${fcidump} | tee ${EZFIO_FILE}.PySCF_EOM_vbarb.out 
