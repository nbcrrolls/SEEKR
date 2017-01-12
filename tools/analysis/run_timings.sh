for anchor in anchor*/md/ens_equil; do
 echo "extracting timings for anchor $anchor"
 cd $anchor
 #python ~/SEEKR/tools/analysis/calc_ens_equil_timings.py 
 cd ../../../
 done
for anchor in anchor*/md/fwd_rev; do
 echo "extracting timings for anchor $anchor"
 cd $anchor
 python ~/SEEKR/tools/analysis/test_calc_fwd_rev_timings.py 
 cd ../../../
 done

