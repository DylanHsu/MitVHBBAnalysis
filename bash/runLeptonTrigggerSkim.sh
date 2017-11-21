inputFileName=$1
outputFileName=$2
selection=$3

origDir=`pwd`
cmsswDir=$CMSSW_BASE/src
#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $cmsswDir
eval `scramv1 runtime -sh`
root -b -l <<EOF
assert(0==gSystem->Load("MitVHBBAnalysis/macros/leptonTriggerSkim_C.so"))
printf("vhbbPlotSkim(\"$inputFileName\",\"$outputFileName\",$selection,false)\n");
assert(leptonTriggerSkim("$inputFileName","$outputFileName",$selection,false));
.q
EOF
cd $origDir
