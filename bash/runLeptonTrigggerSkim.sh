inputFileName=$1
outputFileName=$2
selection=$3

origDir=`pwd`
cmsswDir=$CMSSW_BASE/src
#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $cmsswDir
eval `scramv1 runtime -sh`
xrdcp $inputFileName $(basename $inputFileName)
root -b -l <<EOF
assert(0==gSystem->Load("MitVHBBAnalysis/macros/leptonTriggerSkim_C.so"))
printf("leptonTriggerSkim(\"$(basename $inputFileName)\",\"$outputFileName\",$selection,false)\n")
assert(leptonTriggerSkim("$(basename $inputFileName)","$outputFileName",$selection,false))
.q
EOF
rm $(basename $inputFileName)
cd $origDir
