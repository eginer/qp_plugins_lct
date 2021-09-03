source ~/programs/qp2/quantum_package.rc 
ezfio=f2.ezfio
qp set_file  $ezfio 
core=2
a1=4
a2=10
v1=11
v2=28
core=`qp set_mo_class -q | grep Core | cut -d "(" -f 2 | cut -d ")" -f 1 | sed 's/ //g'`
active=`qp set_mo_class -q | grep Active | cut -d "(" -f 2 | cut -d ")" -f 1 | sed 's/ //g'`
virtual=`qp set_mo_class -q | grep Virtual | cut -d "(" -f 2 | cut -d ")" -f 1 | sed 's/ //g'`
echo $core
echo $active
echo $virtual
qp set_mo_class -c "[${core}]" -a "[${active}]" -d "[${virtual}]"
qp set mu_of_r io_mu_of_r Write
qp set mu_of_r mu_of_r_potential cas_truncated 
qp run write_eff_basis_ints > ${ezfio}.DFT_0.out
grep "TOTAL ENERGY        = " ${ezfio}.DFT_0.out
qp set mu_of_r io_mu_of_r Read
for i in `seq 1 10`
#for i in `seq 1 1`
 do 
 qp set_mo_class -c "[${core}]" -a "[${active}]" -d "[${virtual}]"
 qp run write_eff_basis_ints >  ${ezfio}.DFT_${i}.out
 echo "iteration " $i
 grep "TOTAL ENERGY        = " ${ezfio}.DFT_${i}.out
 qp run fci >  ${ezfio}.fci_${i}.out 
 qp set_mo_class -c "[${core}]" -a "[${active}]" -v "[${virtual}]"
 qp run iter_casscf > ${ezfio}.scf_${i}.out 
done
