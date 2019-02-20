#!/usr/bin/env qpsh

 ezfio=He_av${1}z
 qp create_ezfio -b aug-cc-pv${1}z -m 3 He.xyz  -o $ezfio
 qp run scf
 qp run cis
 qp run save_natorb
 qp run scf

 x=0.
 y=0.
 z=0.

 qp set firth_order_der x_center_wee $x
 qp set firth_order_der y_center_wee $y
 qp set firth_order_der z_center_wee $z

 qp run print_wee > pouet
 out=${ezfio}_${x}_${y}_${z}
 cp fort.33 $out
 mu=`grep "mu_ana =" pouet | cut -d "=" -f 2`


cat << EOF > ${out}.plt
set term pdf
set output "$out.pdf"
set xlabel "r_{12}"
set ylabel "Interaction"
ymax=2*$mu/sqrt(pi) + 1.
set yrange [0:ymax]

plot "${out}" w l title "W_{ee} paper", "" u 1:3 w l title "Fit mu_{num}" , "" u 1:4 w l title "Fit mu_{ana}", 1/x title "exact"

EOF

gnuplot ${out}.plt
