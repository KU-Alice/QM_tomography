input1=${1:"ang_tot.txt"} #mixed  files as input
input2=${2:-"ang_sig.txt"} # output of the code
input3=${3:-"ang_bkg.txt"} # mass limit to apply


do_conversion=0
do_fitting=1
do_qtfitting=1
isMC=1 #mc flag

Data_Dir=$PWD/datafiles
Root_Dir=$PWD/root
Python_Dir=$PWD/python
#echo $Rot_Dir

if [ $do_conversion -eq 1 ]; then
  echo "converting"
  root -l -b -q $Root_Dir/Correlation.C\(\"$1\",\"$2\"\)
else
  echo "not converting"
fi
if [ $do_fitting -eq 1 ]; then
  echo "running fitting "
  if [ $do_qtfitting -eq 1 ]; then
    if [ $isMC -eq 1 ]; then
      echo "this is mc "
      alienv setenv AliPhysics/latest -c root -l -b -q $Root_Dir/totFit.C\(\"$Data_Dir/$1\",\"$Data_Dir/$2\",\"$Data_Dir/$3\"\)
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
