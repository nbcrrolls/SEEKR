#!/bin/bash

for anchor in anchor*/md/ens_equil; do
 echo "extracting equilibrium distribution for anchor $anchor"
 cd $anchor
 cpptraj -i ~/SEEKR/tools/analysis/align.cpptraj
 echo "calculating equilibrium distribution for anchor $anchor"
 python ~/SEEKR/tools/analysis/equil_md.py '../building/holo.parm7' 'ens_equil_center.dcd' 'APN'
 echo "extracting ligand angle for $anchor"
 python ~/SEEKR/tools/analysis/measure_angle.py ../building/holo.parm7 ens_equil_center.dcd 'resname APN and name C47' 'resname APN and name O37' 'resname BCD and name O3 or name O13 or name O23'
 cd ../../../
 done
for anchor in anchor*/md/ens_equil2; do
 echo "extracting equilibrium distribution for anchor $anchor"
 cd $anchor
 cpptraj -i ~/SEEKR/tools/analysis/align.cpptraj
 echo "calculating equilibrium distribution for anchor $anchor"
 python ~/SEEKR/tools/analysis/equil_md.py '../building/holo.parm7' 'ens_equil_center.dcd' 'APN'
 echo "extracting ligand angle for $anchor"
 python ~/SEEKR/tools/analysis/measure_angle.py ../building/holo.parm7 ens_equil_center.dcd 'resname APN and name C47' 'resname APN and name O37' 'resname BCD and name O3 or name O13 or name O23'
 cd ../../../
 done
for anchor in anchor*/md/fwd_rev; do
 echo "extracting FHPD for anchor $anchor"
 cd $anchor
 python ~/SEEKR/tools/analysis/fhpd_md.py '../building/holo.parm7' '../ens_equil/ens_equil_center.dcd' 1 'APN'
 cd ../../../
 done
for anchor in anchor*/md/fwd_rev2; do
 echo "extracting FHPD for anchor $anchor"
 cd $anchor
 python ~/SEEKR/tools/analysis/fhpd_md.py '../building/holo.parm7' '../ens_equil2/ens_equil_center.dcd' 1 'APN'
 cd ../../../
 done


