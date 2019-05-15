dataCardDir=$1 # whbb/testcondor

toMerge=`ls MitVHBBAnalysis/datacards/${dataCardDir}/split/ | grep .root | sed 's/\([[:alnum:]]\+\)\(SR\|CR\)\([[:alnum:]]\?\)_.*/\1\2\3/g' | sort -u`
toMergeVZ=`ls MitVHBBAnalysis/datacards/${dataCardDir}/VZbb/split/ | grep .root | sed 's/\([[:alnum:]]\+\)\(SR\|CR\)\([[:alnum:]]\?\)_.*/\1\2\3/g' | sort -u`

PURPLE='\033[0;35m'
NC='\033[0m' # No Color
for token in $toMerge
do
  echo -e "Merging MitVHBBAnalysis/datacards/${dataCardDir}/split/${token}_*.root into MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root"
  hadd -ff MitVHBBAnalysis/datacards/${dataCardDir}/${token}.root MitVHBBAnalysis/datacards/${dataCardDir}/split/${token}_*.root
done
for token in $toMergeVZ
do
  echo -e "Merging MitVHBBAnalysis/datacards/${dataCardDir}/VZbb/split/${token}_*.root into MitVHBBAnalysis/datacards/${dataCardDir}/VZbb/${token}.root"
  hadd -ff MitVHBBAnalysis/datacards/${dataCardDir}/VZbb/${token}.root MitVHBBAnalysis/datacards/${dataCardDir}/VZbb/split/${token}_*.root
done

toMerge=`ls MitVHBBAnalysis/mva/${dataCardDir}/split/ | grep .root | sed 's/\([[:alnum:]]\+\)\(SR\|CR\)\([[:alnum:]]\?\)_.*/\1\2\3/g' | sort -u`
toMergeVZ=`ls MitVHBBAnalysis/mva/${dataCardDir}/VZbb/split/ | grep .root | sed 's/\([[:alnum:]]\+\)\(SR\|CR\)\([[:alnum:]]\?\)_.*/\1\2\3/g' | sort -u`

PURPLE='\033[0;35m'
NC='\033[0m' # No Color
for token in $toMerge
do
  echo -e "Merging MitVHBBAnalysis/mva/${dataCardDir}/split/${token}_*.root into MitVHBBAnalysis/mva/${dataCardDir}/${token}.root"
  hadd -ff MitVHBBAnalysis/mva/${dataCardDir}/${token}.root MitVHBBAnalysis/mva/${dataCardDir}/split/${token}_*.root
done
for token in $toMergeVZ
do
  echo -e "Merging MitVHBBAnalysis/mva/${dataCardDir}/VZbb/split/${token}_*.root into MitVHBBAnalysis/mva/${dataCardDir}/VZbb/${token}.root"
  hadd -ff MitVHBBAnalysis/mva/${dataCardDir}/VZbb/${token}.root MitVHBBAnalysis/mva/${dataCardDir}/VZbb/split/${token}_*.root
done
