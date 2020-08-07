
set terminal pdf
set output "density-mu.pdf"
set xlabel "distance from O"
set ylabel "density"
plot "H2O_1.e-6.density" u 1:3 w l title "n(r), {/Symbol m}=0", "H2O_0.25.density" u 1:3 w l title "n(r), {/Symbol m}=0.25", "H2O_0.5.density" u 1:3 w l title "n(r), {/Symbol m}=0.5", "H2O_1.0.density" u 1:3 w l title "n(r), {/Symbol m}=1.0", "H2O_1e6.density" u 1:3 w l title "n(r), {/Symbol m}=10^6"
