set terminal pdf enhanced
set output "He_aug-cc-pvdz_0.5.wf.pdf"
set ylabel "Wave functions, {/ Symbol m} = 0.5, AVdz"
set xlabel "r_{12}"
set key bottom right

plot "He_aug-cc-pvdz_0.5.cusp_wf" u 1:3 w l title "exact", "" u 1:4 w l title "FCI", "" u 1:5 w l title "{/ Symbol j}_{/ Symbol m}(r_{12})", "" u 1:7 w l title "e^{j(r_{12},{/ Symbol m)}} {/ Symbol j}_{/ Symbol m}(r_{12})"
