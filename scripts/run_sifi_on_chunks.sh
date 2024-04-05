#!/bin/bash

for i in {1..8}
do
sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/chopping/run00596_single_nlc_all_chunk00${i}.root -i 0 -e 100000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00596_sifi_chunk00${i}.root
done

for i in {1..5}
do
sifi_dst_4to1 0x1000::/scratch1/gccb/data/Jan2023Beam/root/chopping/run00597_single_nlc_all_chunk00${i}.root -i 0 -e 100000000 -p ~/project/sifi-framework/macros/params.txt -o /scratch1/gccb/data/Jan2023Beam/results/run00597_sifi_chunk00${i}.root
done
