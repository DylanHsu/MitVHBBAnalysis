#inputFileName=$1
#outputFileName=$2
#Parse output of catalogT2Prod
outputBasename=$1
sampleType=$2 # Data | MC
xs=$3
inputFileName=$4

outputFileName="/data/t3home000/$USER/fatjetRegTree/${outputBasename}.root"


origDir=`pwd`
cmsswDir=$CMSSW_BASE/src
tempOutput=$origDir/tempOutput_${1}.root
#tempInput=$origDir/tempInput.root
#cp -v $inputFileName $tempInput

#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $cmsswDir
mkdir -p $(dirname "${outputFileName}")
eval `scramv1 runtime -sh`
root -b -l -q MitVHBBAnalysis/macros/fatjetRegTree.C+\(\"$inputFileName\",\"$tempOutput\"\)
cd $origDir
cp -v $tempOutput $outputFileName
rm $tempOutput
