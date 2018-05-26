dataCardDir=$1 # zhbb/testcondor

toMerge=`ls MitVHBBAnalysis/datacards/${dataCardDir}/split/ | grep .root | sed 's/\(ZptBin[[:digit:]]\+\)_.*.root/\1/g' | sed 's/\([[:alnum:]]\+FJ[[:alnum:]]\+\)_.*/\1/g' | sort -u`

PURPLE='\033[0;35m'
NC='\033[0m' # No Color
for token in $toMerge
do
  echo -e "Merging into MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root"
  hadd -v 0 MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root MitVHBBAnalysis/datacards/${dataCardDir}/split/${token}_*.root
done
