inputFileName=$1
outputFileName=$2
selection=$3
sample=$4
useBoosted=$5
modSplitIndex=$6
if [ "$useBoosted" = "" ] 
then
    useBoosted=false
fi
if [ "$modSplitIndex" = "" ] 
then
    modSplitIndex=-1
fi

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
root -b -l <<EOF
bool loadedMacro=(0==gSystem->Load("MitVHBBAnalysis/macros/vhbbPlotSkim_C.so"))
if(!loadedMacro) throw std::runtime_error("Could not load macro shared object MitVHBBAnalysis/macros/vhbbPlotSkim_C.so, go compile it in ACLiC");
printf("vhbbPlotSkim(\"$inputFileName\",\"$tempOutput\",$selection,$sample,false,0,$useBoosted,$modSplitIndex)\n");
vhbbPlotSkim("$inputFileName","$tempOutput",$selection,$sample,false,0,$useBoosted,$modSplitIndex);
.q
EOF
cd $origDir
cp -v $tempOutput $outputFileName
rm $tempOutput
