#!/bin/bash

for anchor in anchor*/md/ens_equil; do
 echo "extracting equilibrium distribution for anchor $anchor"
 cd $anchor
 cpptraj -i ~/SEEKR/tools/analysis/align.cpptraj
 echo "calculating equilibrium distribution for anchor $anchor"
 python ~/SEEKR/tools/analysis/equil_md.py '../building/holo.parm7' 'ens_equil_center.dcd' 'BEN'
 echo "extracting ligand angle for $anchor"
 python ~/SEEKR/tools/analysis/measure_angle.py ../building/holo.parm7 ens_equil_center.dcd 'resname BEN and name C4' 'resname BEN and name C' 'bynum 2479 2490 2500 2536 2719 2746 2770 2788 2795 2868 2927'
 cd ../../../
 done
#for anchor in anchor*/md/fwd_rev; do
# echo "extracting FHPD for anchor $anchor"
# cd $anchor
# python ~/SEEKR/tools/analysis/fhpd_md.py
# cd ../../../
# done

