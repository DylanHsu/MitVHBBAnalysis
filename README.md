# MitVHBBAnalysis

Software to run on Panda Express files for the non-hadronic VH(bb) search

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.
You need to produce Panda Express files using PandaAnalysis with analysis setting complicatedLeptons turned on.

### Dependencies

`CMSSW_9_4_6`
PandaAnalysis
PandaCore
PandaTree
MitAnalysisRunII

### Running Z(ll)H analysis on Condor

The lists of samples being run over are in `config/zllhSamples2016.cfg`, `config/zllhSamples2017.cfg`
Run the following command to make a list of job arguments:
`MitVHBBAnalysis/bash/buildJobArgsZllh.sh` `<subdirectory>` `<MVAVarType>` `<year (2016|2017)>` `<useBoosted (true|false)>`
```
MitVHBBAnalysis/bash/buildJobArgsZllh.sh zhbb/testcondor 3 2016 true
```
This will give you a command to submit the jobs. The job output will show up at `MitVHBBAnalysis/datacards/zhbb/testcondor/split`.
Next, you can merge that output using `MitVHBBAnalysis/bash/mergeOutputZllh.sh`. That will arrive at `MitVHBBAnalysis/datacards/zhbb/testcondor`/

Lastly to make the datacards from the output, you need to load zllhAnalysis.C in ROOT and call `datacardsFromHistograms` as follows:
```
root -b
.L MitVHBBAnalysis/macros/zllhAnalysis.C+g
datacardsFromHistograms("zhbb/testcondor", kZllHSR, true, 3, 0, 2016)
datacardsFromHistograms("zhbb/testcondor", kZllHSR, true, 3, 1, 2016)
datacardsFromHistograms("zhbb/testcondor", kZllHLightFlavorCR, true, 3, 0, 2016)
...
et cetera
```

To make plots use `macros/finalPlot2018.C` passing the plots rootfiles
