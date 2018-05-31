dataCardDir=$1
selection=$2
useBoosted=$3
MVAVarType=$4
year=$5
batchSampleName=$6
batchSampleType=$7
modSplitIndex=$8
if [ "$modSplitIndex" == '' ]
then
  modSplitIndex="-1"
fi

source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsswDir=$CMSSW_BASE/src
cd $cmsswDir
eval `scramv1 runtime -sh`
root -b -l <<EOF
gSystem->Load("MitVHBBAnalysis/macros/whAnalysis_C.so");
whAnalysis("${dataCardDir}",${selection},${useBoosted},${MVAVarType},${year},0,false,"${batchSampleName}",${batchSampleType},${modSplitIndex})
EOF
#whAnalysis("${dataCardDir}",${selection},${useBoosted},${MVAVarType},${year},0,false,"${batchSampleName}",${batchSampleType})

