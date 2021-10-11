source ~/programs/new_qp2/qp2/quantum_package.rc 

basis=aug-cc-pvdz 
mu="hf"
ezfio=He_${basis}_${mu}_rs_ks_scf
rm -rf $ezfio
rm data_$mu
qp create_ezfio -b $basis He.xyz -o $ezfio
qp set_file $ezfio 
qp run scf 
qp set dft_keywords mu_dft_type $mu 
qp run save_mu_of_r_ints 
cd $ezfio/work
for i in *erf*   
do  
 ln -s $i ${i/_erf/} 
done
cd ../../

qp run rs_ks_scf | tee ${ezfio}.out 
