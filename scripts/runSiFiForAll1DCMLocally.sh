#!/bin/bash
# for i in {567..571}
# do
#     echo "run00${i} being sifi-processed..."
# #     sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_spills.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root
#     sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_bg.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi_bg.root
# done
# 
# for i in {575..584}
# do
#     echo "run00${i} being sifi-processed..."
# #     sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_spills.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root
#     sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_bg.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi_bg.root
# done
# 
# # echo "run00597 being sifi-processed..."
# # sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/run00597_single_nlc.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00597_sifi.root
# # 
# # echo "run00596 being sifi-processed..."
# # sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/run00596_single_nlc.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00596_sifi.root


for i in {509..519}
do
    echo "run00${i} being sifi-processed..."
    sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_spills.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root
    sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_bg.root -i 0 -e 1000000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi_bg.root
done
