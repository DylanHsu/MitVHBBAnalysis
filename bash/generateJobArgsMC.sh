argFile=$1
rm $argFile 2>/dev/null
#skimDir=/mnt/hadoop/scratch/dhsu/dylansVHSkims
#skimDir=/mnt/hadoop/scratch/bmaier/panda/008_v1/flat
skimDir=/data/t3home000/dhsu/dylansVHSkims/
#for selType in WHLightFlavorCR WHHeavyFlavorCR WH2TopCR WHSR
#do
  mkdir -p /data/t3home000/dhsu/whbbPlots/root
  for path in `ls $skimDir/SingleElectron*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kData kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/SingleMuon*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kData kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/QCD*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kQCD kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/WZ*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kVZ kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/ZZ*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kVZ kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/WW*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kWW kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/TTbar*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kTT kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/SingleTop*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kTop kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/WJets*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kWjets kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/ZJets*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kZjets kWHLightFlavorCR >> $argFile; done
  for path in `ls $skimDir/W*LNuHbb*.root`; do printf "%s\t%s\t%s\t%s\n" $skimDir/${path##*/} /data/t3home000/dhsu/whbbPlots/root/${path##*/} kVH kWHLightFlavorCR >> $argFile; done
#done

