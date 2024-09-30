#!/bin/bash

# 循环运行Geant4模拟，传入不同的i值
for i in {0..9}
do

  cp CryMu.mac CryMu_$i.mac
  sed -i "s/muCRY_i.root/muCRY_${i}.root/g" CryMu_$i.mac

  # 运行Geant4模拟
  ./muPos CryMu_$i.mac &

done

#rm -rf CryMu_*.mac
