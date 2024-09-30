#! /bin/bash

ls -l /
echo "--------"
ls -l /usr
echo `pwd`
hostname
cat /etc/redhat-release
cat $_CONDOR_JOB_AD
sleep 2
echo "yeah"
echo $2

source ~/.bashrc
source /data/bond/yuxd/root/install/bin/thisroot.sh

cd /home/pku/yuxd/bond/PKMUON_G4sim/DM_step2/analysis/analysis

#root -b -q -l draw.cc > $1
root -b -q -l SaveHist.cc > $1
