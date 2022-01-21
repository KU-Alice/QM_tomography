#! /usr/bin/env python3

filename="$PWD/ratio.txt" # this is file name where s/b passed our requirement stored in 1 line either 0 or 1
#echo "reading file"
n=2 # setting initial logic n =1 means that s/b ratio is in range of our preference
trial=1
Nparticleinitial=1000 # initial number of particle to produce in starlight
#echo $Nparticleinitial
while [ $n -ne 1 ];do #loop to run until s/b ratio is passed
  #echo $n
  Nparticle=$(($trial*1000)) # everytime increasing Nparticle in starlight by 1000 increment

  slightname="slight.in" #reading strlight input parameter to replace Nparticle
  search=$(printf "N_EVENTS=%i" $Nparticleinitial) # the line to find in starlight

  replace=$(printf "N_EVENTS=%i" $Nparticle) # line to change to increase the number of particles
  #echo $replace
  sed -i  "s/$search/$replace/" $slightname # changing the line
  #../starlight # producing starlight with new parameter
  rm starlightgammatomu.dat # removing old file
  root -l -b -q ConvertStarlightAsciiToTree.C\(\"slight.out\",\"starlightgammatomu.dat\"\) # converting starlight output to px1 py1 pz1 E1 px2 py2 pz2 E2

  while IFS= read -r line
  do
  echo "the line is $line"

  n=$line
  done < "$filename"
  alienv setenv AliPhysics/latest -c python3 $PWD/datamerge.py python3 datamerge.py -input1 starlightjpsimumuetacut.dat -input2 starlightgammatomu.dat -ratio 2 -output output.dat $n
  trial=$(($trial+1))
  Nparticleinitial=$Nparticle
  #echo "this is very interesting $n"
  done

#echo $n
#python3 $PWD/datamerge.py python3 datamerge.py -input1 starlightjpsimumuetacut.dat -input2 starlightgammagammamumufinal.dat -ratio 2 -output output.dat 0
