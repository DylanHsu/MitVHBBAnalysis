dataCardDir=$1
selection=$2
useBoosted=$3
MVAVarType=$4
year=$5
vzbbMode=$6
batchSampleName=$7
batchSampleType=$8

source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsswDir=$CMSSW_BASE/src
cd $cmsswDir
eval `scramv1 runtime -sh`
root -b -l <<EOF
gSystem->Load("MitVHBBAnalysis/macros/whAnalysis_C.so");
whAnalysis("${dataCardDir}",${selection},${useBoosted},${MVAVarType},${year},0,false,${vzbbMode},"${batchSampleName}",${batchSampleType})
EOF
#whAnalysis("${dataCardDir}",${selection},${useBoosted},${MVAVarType},${year},0,false,"${batchSampleName}",${batchSampleType})


