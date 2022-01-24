input1=$1 #mixed  files as input
input2=$2 # output of the code
input3=$3 # mass limit to apply


do_conversion=0
do_qtfitting=0
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
if [ $do_qtfitting -eq 1 ]; then
  if [ $isMC -eq 1 ]; then
    echo "this is mc "
    alienv setenv AliPhysics/latest -c root -l -b -q $Root_Dir/totFit.C\(\"$1\",\"$2\",\"$3\"\)
  else
    echo "running fitting on data"
    alienv setenv AliPhysics/latest -c root -l -b -q $Root_Dir/totFit.C\(\"$1\",\"$2\",\)
  fi
else
  echo "not running traditional fitting"
  alienv python3 readhistogram.py alienv setenv AliPhysics/latest -c python3 $Python_Dir/readhistogram.py
fi
