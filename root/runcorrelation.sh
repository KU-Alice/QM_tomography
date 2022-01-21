input1=$1 #mixed  files as input
input2=$2 # output of the code
input3=$3 # mass limit to apply
root -l -b -q Correlation.C\(\"$1\",\"$2\",$3\)
