mutab="1.e-6   0.25  0.5  1.0  2.0  5.0  1e6"
#mutab="1.e-6"
source ~/programs/qp2/quantum_package.rc
ezfio=H2O

for mu in $mutab
do
 newdir=${ezfio}_$mu
 qp set_file $newdir
 qp run plot_density | tee ${newdir}.data
done
