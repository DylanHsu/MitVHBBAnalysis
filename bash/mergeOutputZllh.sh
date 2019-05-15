dataCardDir=$1 # zhbb/testcondor

toMerge=`ls MitVHBBAnalysis/datacards/${dataCardDir}/split/ | grep .root | sed 's/\(ZptBin[[:digit:]]\+\)_.*.root/\1/g' | sed 's/\([[:alnum:]]\+FJ[[:alnum:]]\+\)_.*/\1/g' | sort -u`

PURPLE='\033[0;35m'
NC='\033[0m' # No Color
for token in $toMerge
do
  echo -e "Merging into MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root"
  hadd -ff MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root MitVHBBAnalysis/datacards/${dataCardDir}/split/${token}_*.root
done

toMerge=`ls MitVHBBAnalysis/mva/${dataCardDir}/split/ | grep .root | sed 's/\(ZptBin[[:digit:]]\+\)_.*.root/\1/g' | sed 's/\([[:alnum:]]\+FJ[[:alnum:]]\+\)_.*/\1/g' | sort -u`

PURPLE='\033[0;35m'
NC='\033[0m' # No Color
for token in $toMerge
do
  echo -e "Merging into MitVHBBAnalysis/mva/${dataCardDir}/${token}.root"
  hadd -ff MitVHBBAnalysis/mva/${dataCardDir}/${token}.root MitVHBBAnalysis/mva/${dataCardDir}/split/${token}_*.root
done
