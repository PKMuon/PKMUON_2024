#!/bin/bash

mkdir -p root_file

# 循环运行Geant4模拟，传入不同的i值
for i in $(seq 1 $(lscpu | grep '^CPU(s):' | grep -o '[0-9]\+')); do

  cp SingleEngMu.mac SingleEngMu_$i.mac
  sed -i "s/Mu_1GeV\\.root/Mu_1GeV_${i}.root/g" SingleEngMu_$i.mac

  # 运行Geant4模拟
  ./muPos SingleEngMu_$i.mac &

done

#rm -rf SingleEngMu_*.mac
