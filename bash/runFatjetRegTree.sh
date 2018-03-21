inputFileName=$1
outputFileName=$2

origDir=`pwd`
cmsswDir=$CMSSW_BASE/src
tempOutput=$origDir/tempOutput.root
#tempInput=$origDir/tempInput.root
#cp -v $inputFileName $tempInput

#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $cmsswDir
mkdir -p $(dirname "${2}")
eval `scramv1 runtime -sh`
root -b -l -q MitVHBBAnalysis/macros/fatjetRegTree.C+\(\"$inputFileName\",\"$tempOutput\",vhbbPlot::kVH\)
cd $origDir
cp -v $tempOutput $outputFileName
rm $tempOutput
