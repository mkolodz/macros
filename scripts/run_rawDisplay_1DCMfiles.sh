#!/bin/bash
# for i in {567..571}
# do
#     echo "run00${i} being rD-processed..."
#     root -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly.C({\"/scratch1/gccb/data/Jan2023Beam/results/00${i}_sifi.root\"})"
# #     sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_spills.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root
# done
# 
# for i in {575..584}
# do
#     echo "run00${i} being rD-processed..."
#     root -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly.C({\"/scratch1/gccb/data/Jan2023Beam/results/00${i}_sifi.root\"})"
# #     sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_spills.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root
# done

# echo "run00597 being rD-processed..."
# sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/run00597_single_nlc.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00597_sifi.root
# 
# echo "run00596 being rD-processed..."
# sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/run00596_single_nlc.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00596_sifi.root

for i in {1..8}
do
# sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/chopping/run00596_single_nlc_all_chunk00${i}.root -i 0 -e 100000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00596_sifi_chunk00${i}.root
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly.C({\"/scratch1/gccb/data/Jan2023Beam/results/sifi_run00596_chunk00${i}.root\"})"
#     root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly.C({\"/scratch1/gccb/data/Jan2023Beam/results/run00571_sifi.root\"})"

done

for i in {1..5}
do
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly.C({\"/scratch1/gccb/data/Jan2023Beam/results/sifi_run00597_chunk00${i}.root\"})"

# # sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/chopping/run00597_single_nlc_all_chunk00${i}.root -i 0 -e 100000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00597_sifi_chunk00${i}.root
done




