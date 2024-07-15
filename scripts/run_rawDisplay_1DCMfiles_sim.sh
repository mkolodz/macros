#!/bin/bash
# root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/sim_deleteme.root\"})"

#FILES="sifi_filtered_SystemMatrix_CodedMaskHIT_simv5_Pixel*_0to39.root"

# for i in {0..199}
for i in {0..199}
do
    echo "Pixel${i}"
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch2/gccb/magda/results_sim/sifi_filtered_SystemMatrix_CodedMaskHIT_simv5_Pixel${i}_0to39.root\"})"
#     sifi_dst_4to1_sim ${f} -p ~/project/sifi-framework/macros/params_4to1.txt -o /scratch2/gccb/magda/results_sim/sifi_${f}.root
done
