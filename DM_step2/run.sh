#!/bin/bash

# 循环运行Geant4模拟，传入不同的i值
for i in {0..19}
do

  cp SingleEngMu.mac SingleEngMu_$i.mac
  sed -i "" "s/mu1GeV_i.root/mu1GeV_${i}.root/g" SingleEngMu_$i.mac

  # 运行Geant4模拟
  ./muPos SingleEngMu_$i.mac &

done

#rm -rf SingleEngMu_*.mac
