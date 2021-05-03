
set terminal pdf
set output "on-top-mu.pdf"
set xlabel "distance from O"
set ylabel "on-top pair density"
plot "H2O_1.e-6.density" w l title "n2(r,r), {/Symbol m}=0", "H2O_0.25.density" w l title "n2(r,r), {/Symbol m}=0.25", "H2O_0.5.density" w l title "n2(r,r), {/Symbol m}=0.5", "H2O_1.0.density" w l title "n2(r,r), {/Symbol m}=1.0", "H2O_1e6.density" w l title "n2(r,r), {/Symbol m}=10^6"
