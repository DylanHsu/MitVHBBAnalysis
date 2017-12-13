argFile=$1
rm $argFile 2>/dev/null
#skimDir=/mnt/hadoop/scratch/dhsu/dylansVHSkims
#skimDir=/mnt/hadoop/scratch/bmaier/panda/008_v1/flat
skimDir=/data/t3home000/dhsu/benediktsVHSkims/split
#for selType in WHLightFlavorCR WHHeavyFlavorCR WH2TopCR WHSR
#do
  mkdir -p /data/t3home000/dhsu/whbbPlots/root
  for path in `ls $skimDir/SingleElectron_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kData kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/SingleMuon_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kData kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/QCD_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kQCD kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/WZ_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kVZ kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/ZZ_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kVZ kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/WW_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kWW kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/TTbar_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kTT kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/SingleTop_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kTop kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/WJets_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kWjets kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/ZJets_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kZjets kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/VH_*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kVH kWHLightFlavorCR >> $argFile; done
#done

