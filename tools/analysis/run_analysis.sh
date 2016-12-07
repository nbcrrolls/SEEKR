#!/bin/bash

for anchor in anchor*/md/ens_equil; do
 echo "extracting equilibrium distribution for anchor $anchor"
 cd $anchor
 cpptraj -i ~/seekr/tools/analysis/align.cpptraj
 python ~/seekr/tools/analysis/equil_md.py
 echo "extracting ligand angle for $anchor"
 python ~/seekr/tools/analysis/measure_angle.py
 cd ../../../
 done
for anchor in anchor*/md/fwd_rev; do
 echo "extracting FHPD for anchor $anchor"
 cd $anchor
 python ~/seekr/tools/analysis/fhpd_md.py
 cd ../../../
 done

