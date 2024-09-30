bkgtype=("air" "vac")
sigtype=("const" "maxwell")
DMmass=("0p005" "0p05" "0p1" "0p2" "0p5" "1" "10" "100")

mypath=$PWD
echo $mypath

for (( i=0; i<${#sigtype[@]}; i++ ))
do
  for (( j=0; j<${#bkgtype[@]}; j++ ))
  do
    echo signal velocity: ${sigtype[$i]} and background material: ${bkgtype[$j]}
    if [ -d "UL_${sigtype[$i]}_${bkgtype[$j]}" ]; then
      rm -rf UL_${sigtype[$i]}_${bkgtype[$j]}
    fi
    mkdir UL_${sigtype[$i]}_${bkgtype[$j]}
    cd UL_${sigtype[$i]}_${bkgtype[$j]}
    for (( k=0; k<${#DMmass[@]}; k++ ))
    do
      echo Dark Matter mass ${DMmass[$k]}
      cp -r ../datacard.txt datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}.txt
      sed -i "s/SIG/${sigtype[$i]}_${DMmass[$k]}/g" datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}.txt
      sed -i "s/BKG/${bkgtype[$j]}/g" datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}.txt
      combine -M AsymptoticLimits datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}.txt -n ${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]} --run blind > log_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]} &
    done
    cd $mypath
  done
done
