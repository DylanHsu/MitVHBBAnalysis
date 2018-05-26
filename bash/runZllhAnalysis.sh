dataCardDir=$1
selection=$2
useBoosted=$3
MVAVarType=$4
binZpt=$5
year=$6
batchSampleName=$7
batchSampleType=$8
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsswDir=$CMSSW_BASE/src
cd $cmsswDir
eval `scramv1 runtime -sh`
root -b -l <<EOF
gSystem->Load("MitVHBBAnalysis/macros/zllhAnalysis_C.so");
zllhAnalysis("${dataCardDir}",${selection},${useBoosted},${MVAVarType},${binZpt},${year},0,false,"${batchSampleName}",${batchSampleType})
EOF

