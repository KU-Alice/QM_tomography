#! /usr/bin/env python3
#######################important###########################
#######################################################
###############################################
###############################
####for this script to work initially slight.in need to be set to produce 1000 particles##########
filename=$PWD/ratio.txt # this is file name where s/b passed our requirement stored in 1 line either 0 or 1

#Some required directory can be changed according to the
Signal_FIle_Dir=$PWD/data/incohjpsi
Bkg_File_Dir=$PWD/data/gammagamma
Mixed_File_Dir=$PWD/data/mixed_data
Starlight_Dir=/home/agautam/Run2Analysis/Starlight/starlight_installed/
#echo "reading file"
n=2 # setting initial logic n =1 means that s/b ratio is in range of our preference
trial=0
Nparticleinitial=1000 # previous number of particle  in slignt.in
Nstart=180000 # initial number
#FileName_Sig="$PWD/Starlight."

#echo $Nparticleinitial
while [ $n -ne 1 ];do #loop to run until s/b ratio is passed
  echo $n
  Nparticle=$(((($trial*1000))+$Nstart)) # everytime increasing Nparticle in starlight by 1000 increment

  slightname="slight.in" #reading strlight input parameter to replace Nparticle
  search=$(printf "N_EVENTS=%i" $Nparticleinitial) # the line to find in starlight

  replace=$(printf "N_EVENTS=%i" $Nparticle) # line to change to increase the number of particles
  #echo $replace
  sed -i  "s/$search/$replace/" $Bkg_File_Dir/$slightname # changing the line
  (cd $Bkg_File_Dir && $Starlight_Dir/starlight) # producing starlight with new parameter
  rm starlightgammatomubig.dat # removing old file
  root -l -b -q ConvertStarlightAsciiToTree.C\(\"$Bkg_File_Dir\/slight.out\",\"$Bkg_File_Dir\/starlightgammatomubig.dat\"\) # converting starlight output to px1 py1 pz1 E1 px2 py2 pz2 E2

  while IFS= read -r line
  do
    echo "the line is $line"

    n=$line
    done < "$filename"
    alienv setenv AliPhysics/latest -c python3 $PWD/datamerge.py -input1 $Signal_FIle_Dir/incohjpsitomumu.dat -input2 $Bkg_File_Dir/starlightgammatomubig.dat -ratio 2 -output $Mixed_File_Dir/outputbig.dat $n
    trial=$(($trial+1))
    Nparticleinitial=$Nparticle
    #echo "this is very interesting $n"
  done

#echo $n
#python3 $PWD/datamerge.py python3 datamerge.py -input1 starlightjpsimumuetacut.dat -input2 starlightgammagammamumufinal.dat -ratio 2 -output output.dat 0
