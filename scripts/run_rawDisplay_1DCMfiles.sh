#!/bin/bash
for i in {567..571}
do
    echo "run00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root\"})"
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi_bg.root\"})"
done

for i in {575..584}
do
    echo "run00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root\"})"
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi_bg.root\"})"
done

for i in {1..8}
do
    echo "run00596_chunk00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/sifi_run00596_chunk00${i}.root\"})"
done

for i in {1..5}
do
    echo "run00597_chunk00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/sifi_run00597_chunk00${i}.root\"})"
done

for i in {509..519}
do
    echo "run00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi.root\"})"
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/run00${i}_sifi_bg.root\"})"
done

for i in {1..4}
do
    echo "run00507_chunk00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/sifi_run00507_chunk00${i}.root\"})"
done

for i in {1..6}
do
    echo "run00520_chunk00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/sifi_run00520_chunk00${i}.root\"})"
done

for i in {1..6}
do
    echo "run00523_chunk00${i} being rD-processed..."
    root -b -l -q "/home/magda/project/macros/rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({\"/scratch1/gccb/data/Jan2023Beam/results/sifi_run00523_chunk00${i}.root\"})"
done
