inputFileName=$1
outputFileName=$2
selection=$3

origDir=`pwd`
cmsswDir=$CMSSW_BASE/src
#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $cmsswDir
eval `scramv1 runtime -sh`
tempInputFile=$origDir/$(basename $inputFileName)
tempOutputFile=$origDir/triggerSkim_$(basename $inputFileName)
xrdcp $inputFileName $tempInputFile
root -b -l <<EOF
assert(0==gSystem->Load("MitVHBBAnalysis/macros/leptonTriggerSkim_C.so"))
printf("leptonTriggerSkim(\"$tempInputFile\",\"$tempOutputFile\",$selection,false)\n")
assert(leptonTriggerSkim("$tempInputFile","$tempOutputFile",$selection,false))
.q
EOF
cp $tempOutputFile $outputFileName
rm $tempInputFile $tempOutputFile
cd $origDir
