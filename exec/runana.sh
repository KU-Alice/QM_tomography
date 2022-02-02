input1=${1:-"ang_tot.txt"} #mixed  files as input
input2=${2:-"ang_sig.txt"} # output of the code
input3=${3:-"ang_bkg.txt"} # mass limit to apply
#slicing the inputs name
sliced_input1=$(echo ${input1:0:7})
sliced_input2=$(echo ${input2:0:7})
sliced_input3=$(echo ${input3:0:7})
echo "testing input file name"
#echo ${input1:0:7}
echo $1
#echo $sliced_input1
do_conversion=0
do_fitting=1
do_qtfitting=1
isMC=1 #mc flag

Data_Dir=$PWD/datafiles
Root_Dir=$PWD/root
Python_Dir=$PWD/python
#echo $Rot_Dir

sig_file="incohjpsitomumu.dat"
bkg_file="starlightgammatomubig"
mix_file="incohjpsi_gammagamma_to_mumu.dat"
if [ $do_conversion -eq 1 ]; then
  echo "converting"

  root -l -b -q $Root_Dir/Correlation.C\(\"$Data_Dir\/$mix_file\",\"$Data_Dir\/$sliced_input1\",0\)
  root -l -b -q $Root_Dir/Correlation.C\(\"$Data_Dir\/$mix_file\",\"$Data_Dir\/$sliced_input3\",4\)
  root -l -b -q $Root_Dir/Correlation.C\(\"$Data_Dir\/$sig_file\",\"$Data_Dir\/$sliced_input2\",0\)
else
  echo "not converting"
fi
if [ $do_fitting -eq 1 ]; then
  echo "running fitting "
  if [ $do_qtfitting -eq 1 ]; then
    if [ $isMC -eq 1 ]; then
      echo "this is mc "
      alienv setenv AliPhysics/latest -c root -l  $Root_Dir/totFit.C\(\"$Data_Dir\/$input1\",\"$Data_Dir\/$input2\",\"$Data_Dir\/$input3\"\)
      #alienv setenv AliPhysics/latest -c root -l  testplot.C\(\"$Data_Dir\/$input1\"\)
    else
      echo "running fitting on data"
      alienv setenv AliPhysics/latest -c root -l -b -q $Root_Dir/totFit.C\(\"$1\",\"$2\",\)
    fi
  else
    echo "This is data now running traditional fitting"
    alienv setenv AliPhysics/latest -c python3 $Python_Dir/readhistogram.py
  fi
else
  echo "not runnning fitting "
fi
