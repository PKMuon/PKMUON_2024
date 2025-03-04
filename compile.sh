#!/bin/bash -ve

source /data/bond/yuxd/root/install/bin/thisroot.sh

g++ -g muPos.cc src/*.cc -o muPos -w -Iinclude -I$HOME/.local/include $(root-config --cflags --libs) $(geant4-config --cflags --libs)
