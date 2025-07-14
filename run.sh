rm -f log*.txt
#rm trust.txt ch.txt 
rm -f AP_stats.txt D.txt agents_groups.txt groups.txt average_tau.txt stat_tau.txt 
octave --no-gui main_simulator.m | tee log`date +%s`.txt
