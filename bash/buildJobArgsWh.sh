dataCardDir=$1 # zhbb/testcondor
MVAVarType=$2  # 1|3
year=$3        # 2016|2017
useBoosted=$4  # true|false

resolvedSelections="kWHLightFlavorCR kWHHeavyFlavorCR kWH2TopCR kWHSR"
boostedSelections="kWHLightFlavorFJCR kWHHeavyFlavorFJCR kWHTT1bFJCR kWHTT2bFJCR kWHFJSR"
resolvedSelectionsVZ="kWHLightFlavorCR kWHHeavyFlavorCR kWH2TopCR kWHVZbbCR"
boostedSelectionsVZ="kWHLightFlavorFJCR kWHHeavyFlavorFJCR kWHTT1bFJCR kWHTT2bFJCR kWHVZbbFJCR"
resolvedSelectionsVZ=""
boostedSelectionsVZ=""

mkdir -p "MitVHBBAnalysis/datacards/${dataCardDir}/VZbb"
jobArgsFile="MitVHBBAnalysis/datacards/${dataCardDir}/jobArgs.txt"
rm $jobArgsFile 2>/dev/null
config="MitVHBBAnalysis/config/whSamples${year}.cfg"

PURPLE='\033[0;35m'
RED='\033[0;31m'
NC='\033[0m' # No Color

if [ "$useBoosted" == 'true' ] ; then
  echo -e "${RED}You have chosen to use a boosted category, amazing!${NC}"
elif [ "$useBoosted" == 'false' ] ; then
  echo -e "${RED}You are running without the boosted category${NC}"
else 
  echo -e "${RED}4th argument must be (true|false), exiting${NC}"
  exit
fi
echo -e "${RED}Compiling the macro...${NC}"
root -b -l << EOF
 .L MitVHBBAnalysis/macros/whAnalysis.C+g
EOF
  
#Resolved category
while read -r line
do
  a=( $line )
  batchSampleName=${a[0]}
  batchSampleType=${a[1]}
  for sel in $resolvedSelections
  do
    echo "${dataCardDir} ${sel} ${useBoosted} ${MVAVarType} ${year} false ${batchSampleName} ${batchSampleType} ${idx}" >> $jobArgsFile
  done
  for sel in $resolvedSelectionsVZ
  do
    echo "${dataCardDir} ${sel} ${useBoosted} ${MVAVarType} ${year} true ${batchSampleName} ${batchSampleType} ${idx}" >> $jobArgsFile
  done
done < $config

#Boosted category
if [ "$useBoosted" == 'true' ]
then
  while read -r line
  do
    a=( $line )
    batchSampleName=${a[0]}
    batchSampleType=${a[1]}
    for sel in $boostedSelections
    do
      echo "${dataCardDir} ${sel} ${useBoosted} ${MVAVarType} ${year} false ${batchSampleName} ${batchSampleType} ${idx}" >> $jobArgsFile
    done
    for sel in $boostedSelectionsVZ
    do
      echo "${dataCardDir} ${sel} ${useBoosted} ${MVAVarType} ${year} true ${batchSampleName} ${batchSampleType} ${idx}" >> $jobArgsFile
    done
  done < $config
fi

echo "Done building the job arguments for $dataCardDir area"
echo -e "Remove any jobs you don't want from the file ${PURPLE}${jobArgsFile}${NC}"
echo "When you're ready to submit, do"
echo -e "  ${RED}PandaCore/bin/submit --exec MitVHBBAnalysis/bash/runWhAnalysis.sh --arglist $jobArgsFile --cache /data/t3serv014/$USER/submit/${dataCardDir}${NC} --njobs 1000"
echo "You can check on your jobs using"
echo -e "  ${RED}PandaCore/bin/check --cache /data/t3serv014/$USER/submit/${dataCardDir}${NC}"
echo -e "Plots and histograms output will arrive at ${PURPLE}MitVHBBAnalysis/datacards/${dataCardDir}/split${NC}"
