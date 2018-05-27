dataCardDir=$1 # whbb/testcondor

toMerge=`ls MitVHBBAnalysis/datacards/${dataCardDir}/split/ | grep .root | sed 's/\([[:alnum:]]\+\)\(SR\|CR\)\([[:alnum:]]\?\)_.*/\1\2\3/g' | sort -u`

PURPLE='\033[0;35m'
NC='\033[0m' # No Color
for token in $toMerge
do
  echo -e "Merging MitVHBBAnalysis/datacards/${dataCardDir}/split/${token}_*.root into MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root"
  hadd -ff MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root MitVHBBAnalysis/datacards/${dataCardDir}/split/${token}_*.root
done
