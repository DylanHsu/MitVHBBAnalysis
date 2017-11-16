# MitVHBBAnalysis

Software to run on Panda Express files for the non-hadronic VH(bb) search

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.
You need to produce Panda Express files using PandaAnalysis with analysis setting complicatedLeptons turned on.

### Dependencies

PandaAnalysis
PandaCore
PandaTree

### Making plotting trees from Panda Express

To make a plot skim you can use the ROOT macro `macros/vhbbPlotSkim.C` which compiles via ACLiC.
You can make plotting trees by running the following example command from your $CMSSW\_BASE/src folder:
```
MitVHBBAnalysis/runPlotSkim.sh /mnt/hadoop/scratch/dhsu/dylansVHSkims/WJets_ht1200to2500.root myWJetsSkim.root kWjets kWHHeavyFlavorCR
```
The process types e.g. `kWjets` and selection types e.g. `kWHHeavyFlavorCR` are enumerated in `macros/vhbbPlot.h`.
You can run a bunch of these commands on the cluster but first you need to create a file containing job arguments of the form:
```
<input panda express file> <output skim file> <enumerated process type> <selection type>
```

Once you have this job arguments file, you can run it (again from $CMSSW\_BASE/src) using the following command:
```
PandaCore/bin/submit MitVHBBAnalysis/runPlotSkim.sh /path/tojobArgs.txt
```

Lastly, merge all your plot skims together for a given selection type and then plot them using `macros/finalPlot2018.C`.

