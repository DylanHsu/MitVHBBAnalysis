#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <mutex>
#include <sstream>
#include <thread>
#include <unistd.h>

#include "PandaAnalysis/Flat/interface/GeneralTree.h"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include <Compression.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeIndex.h>
#include <TVector2.h>

#include "formulas.h"
#include "TMVA/Reader.h"
#include "MitAnalysisRunII/panda/macros/80x/common.h"
#include "vhbbPlot.h"
#include "PandaAnalysis/Flat/interface/Common.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// for sel in kWHSR kWHLightFlavorCR kWHHeavyFlavorCR kWH2TopCR; do for i in 0 1; do root -b -l -q MitVHBBAnalysis/macros/nwhAnalysis.C+\(\"zhbb/test\",${sel},false,3,${i},2016,0,true\) & done; done

const bool useHtBinnedVJetsKFactor=true;
const int NJES = (int)shiftjes::N; // Number of JES variations
const int nLepSel=2; // Number of lepton selections
const int nPlots=50; // Max number of plots
const int nThreads=10;
float sf_training=1.4286;
vector<TString> leptonStrings={"mn","en"};
using namespace vhbbPlot;
std::mutex mvaTreeMutex, jecAk4UncMutex, jecAk8UncMutex;

struct analysisObjects {
  // Physics
  selectionType selection;
  std::map<selectionType,vector<TString>> cuts;
  float isojetBtagCut;
  // Configuration parameters
  bool useBoostedCategory, vzbbMode;
  int MVAVarType; 
  unsigned year;
  unsigned debug;
  double lumi;
  unsigned long whichTriggers;
  // Plotting and histograms parameters
  vector<float> xmin, xmax;
  vector<int> nbins;
  vector<TString> histoNames, histoTitles;
  vector<float> MVAbins;
  TString MVAVarName, shapeType; // title of histos, and simple name
  // Reading/writing trees
  // Histogram pointers
  TH1F *histos[nLepSel][nPlots][nPlotCategories];
  TH1F *histo_Baseline                           [nLepSel][nPlotCategories];  
  TH1F *histo_pileupUp                           [nLepSel][nPlotCategories];
  TH1F *histo_pileupDown                         [nLepSel][nPlotCategories];
  TH1F *histo_VHCorrUp                           [nLepSel][nPlotCategories];
  TH1F *histo_VHCorrDown                         [nLepSel][nPlotCategories];
  TH1F *histo_QCDr1f2                            [nLepSel][nPlotCategories];
  TH1F *histo_QCDr1f5                            [nLepSel][nPlotCategories];
  TH1F *histo_QCDr2f1                            [nLepSel][nPlotCategories];
  TH1F *histo_QCDr2f2                            [nLepSel][nPlotCategories];
  TH1F *histo_QCDr5f1                            [nLepSel][nPlotCategories];
  TH1F *histo_QCDr5f5                            [nLepSel][nPlotCategories];
  TH1F *histo_QCDScaleUp                         [nLepSel][nPlotCategories];
  TH1F *histo_QCDScaleDown                       [nLepSel][nPlotCategories];
  TH1F *histo_eleSFUp                            [nLepSel][nPlotCategories];
  TH1F *histo_eleSFDown                          [nLepSel][nPlotCategories];
  TH1F *histo_muSFUp                             [nLepSel][nPlotCategories];
  TH1F *histo_muSFDown                           [nLepSel][nPlotCategories];
  TH1F *histo_btag[GeneralTree::nCsvShifts][5][3][nLepSel][nPlotCategories];
  TH1F *histo_VGluUp                             [nLepSel][nPlotCategories];
  TH1F *histo_VGluDown                           [nLepSel][nPlotCategories];
  TH1F *histo_jes[(int)shiftjes::N]              [nLepSel][nPlotCategories];
  TH1F *histo_doubleBUp                          [nLepSel][nPlotCategories];
  TH1F *histo_doubleBDown                        [nLepSel][nPlotCategories];
  // Corrections storage
  TH1D *puWeights=0, *puWeightsUp=0, *puWeightsDown=0;
  TH2D *kfactors_ZJets=0, *kfactors_WJets=0;
  BTagCalibrationReader *deepcsvSFs=0; 
  BTagCalibration *deepcsvCalib=0; 
  vector<double> ZjetsEWKCorr, WjetsEWKCorr;
  
  vector<JetCorrectionUncertainty*> jecUncsAK4, jecUncsAK8;
  std::map<TString, JetCorrectionUncertainty*> jecUncSourcesAK4, jecUncSourcesAK8;
  
  TTree *mvaTree=0; TFile *mvaFile=0;
  // Common vars
  float mva_ZBosonPt, mva_ZBosonM, mva_CosThetaCS, mva_CosThetaStar;
  float mva_weight;
  unsigned char mva_category;
  ULong64_t mva_eventNumber;
  float mva_lepton1Pt, mva_lepton1Eta, mva_lepton1Charge;
  float mva_mT, mva_WBosonPt, mva_dPhil1W;
  // Resolved only vars
  float mva_sumEtSoft1, mva_nSoft2, mva_nSoft5, mva_nSoft10;
  float mva_bjet1Pt, mva_bjet2Pt, mva_bjet1btag, mva_bjet2btag;
  float mva_hbbpt, mva_hbbm, mva_dPhiWH, mva_ptBalanceWH;
  float mva_dPhiLep1Met, mva_topMass, mva_pfmet;
  float mva_dRBjets, mva_dEtaBjets, mva_dEtaLep1H;
  float mva_nAddJet;
  // Boosted only vars
  float mva_nIsojet, mva_MSD, mva_Tau21SD, mva_Tau32SD, mva_fjPt;
  float mva_psi022004031003, mva_psi022004022003, mva_psi022004030503;
  float mva_ptBalanceWHFJ, mva_dEtaLep1FJ, mva_dPhiWHFJ;
  float mva_HTTFRec;
  // MVA output
  vector<TMVA::Reader*> reader;
  float mvaInputs[nThreads][16];
  
};

void analyzeSample(pair<TString,vhbbPlot::sampleType> sample, TTree *events, analysisObjects &ao, int split=-1);
void writeDatacards(analysisObjects &ao, TString dataCardDir);

// useBoostedCategory:
//   False means you will use all possible events to do the resolved Z(ll)H(bb)
//   True means the events with Z back to back with 250 GeV fatjet are reserved for boosted ZH
// MVAVarType:
//   1 - simple kinematic variable
//   2 - multiclass BDT (not used)
//   3 - simple BDT

void whAnalysis(
  TString dataCardDir,
  vhbbPlot::selectionType selection,
  bool useBoostedCategory=false,
  int MVAVarType=3,
  unsigned year=2016,
  unsigned debug=0,
  bool multithread=false,
  bool vzbbMode=false,
  TString batchSampleName="",
  vhbbPlot::sampleType batchSampleType=vhbbPlot::kData
) {
  struct analysisObjects ao;
  ao.MVAVarType=MVAVarType;
  ao.year=year;
  ao.debug=debug;
  ao.lumi=(year==2016)? 35900:41500;
  ao.useBoostedCategory=useBoostedCategory;
  ao.selection = selection;
  ao.vzbbMode = vzbbMode;
  TString ntupleDir2016 = multithread?
    "/data/t3home000/dhsu/dylansVHSkims/2016/v_009_vhbb3":
    "/mnt/hadoop/scratch/dhsu/dylansVHSkims/2016/v_009_vhbb3/split";
  TString ntupleDir2017 = multithread?
    "/data/t3home000/dhsu/dylansVHSkims/2017/v_010_vhbb3":    
    "/mnt/hadoop/scratch/dhsu/dylansVHSkims/2017/v_010_vhbb3/split";    
  TString ntupleDir = (year==2016)? ntupleDir2016:ntupleDir2017;

  // Analysis Cuts
  ao.isojetBtagCut = (ao.year==2017)? deepcsvLoose : cmvaLoose;
  ao.cuts[kWHLightFlavorCR      ] ={"boostedVeto","WpT","pTjj",           "dPhiLep1Met","looseBTag","mediumBVeto","metSig"};
  ao.cuts[kWHHeavyFlavorLoMassCR] ={"boostedVeto","WpT","pTjj",           "dPhiLep1Met","2jets"    ,"tightBTag","mjjSBLo","metSig"};
  ao.cuts[kWHHeavyFlavorHiMassCR] ={"boostedVeto","WpT","pTjj",           "dPhiLep1Met","2jets"    ,"tightBTag","mjjSBHi","metSig"};
  ao.cuts[kWH2TopCR             ] ={"boostedVeto","WpT","pTjj",           "dPhiLep1Met","4+jets"   ,"tightBTag","lowMET"};
  ao.cuts[kWHVZbbCR             ] ={"boostedVeto","WpT","pTjj","dPhiWH"  ,"dPhiLep1Met","2-3jets"  ,"tightBTag","looseBTag2","mjjVZ"};
  ao.cuts[kWHSR                 ] ={"boostedVeto","WpT","pTjj","dPhiWH"  ,"dPhiLep1Met","2-3jets"  ,"tightBTag","looseBTag2","mjj"};
  ao.cuts[kWHLightFlavorCR      ] ={"boostedVeto","WpT","pTjj",           "dPhiLep1Met","looseBTag","mediumBVeto","metSig"};
  ao.cuts[kWHLightFlavorFJCR    ] ={"boostedCat" ,"WpTFJ","pTFJ","dPhiWHFJ","bvetoFJ","0ijb"            };
  ao.cuts[kWHHeavyFlavorFJCR    ] ={"boostedCat" ,"WpTFJ","pTFJ","dPhiWHFJ","btagFJ" ,"0ijb","mSD_SB"   };
  ao.cuts[kWHTT2bFJCR           ] ={"boostedCat" ,"WpTFJ","pTFJ","dPhiWHFJ","btagFJ" ,"1ijb"            };
  ao.cuts[kWHTT1bFJCR           ] ={"boostedCat" ,"WpTFJ","pTFJ","dPhiWHFJ","bvetoFJ","1ijb"            };
  ao.cuts[kWHVZbbFJCR           ] ={"boostedCat" ,"WpTFJ","pTFJ","dPhiWHFJ","btagFJ" ,"0ijb","mSDVZ_SR" };
  ao.cuts[kWHFJSR               ] ={"boostedCat" ,"WpTFJ","pTFJ","dPhiWHFJ","btagFJ" ,"0ijb","mSD_SR"   };
  ao.cuts[kWHFJPresel           ] ={"boostedCat" ,"WpTFJ","pTFJ","dPhiWHFJ"                             };
  /////////////////////////////
  // List of Samples
  vector<pair<TString,vhbbPlot::sampleType>> samples;
  bool isBatchMode= (batchSampleName!="");
  TString batchSuffix="";
  if(vzbbMode) 
    dataCardDir = dataCardDir + "/VZbb";
  if(isBatchMode) {
    // Handle batch mode for Condor
    //multithread=false; // force single threading
    dataCardDir = dataCardDir + "/split"; // write output in a subdirectory
    batchSuffix = "_"+batchSampleName; // add a suffix to the output with this sample's name
    samples.emplace_back(batchSampleName, batchSampleType);
  } else if(year==2016) {
    samples.emplace_back("SingleElectron"                 , vhbbPlot::kData   );
    samples.emplace_back("SingleMuon"                     , vhbbPlot::kData   );
    samples.emplace_back("QCD_ht100to200"                 , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht200to300"                 , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht300to500"                 , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht500to700"                 , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht700to1000"                , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht1000to1500"               , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht1500to2000"               , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht2000toinf"                , vhbbPlot::kQCD    );
    samples.emplace_back("WZTo1L1Nu2Q"                    , vhbbPlot::kVZ     );
    samples.emplace_back("WZTo1L3Nu"                      , vhbbPlot::kVZ     );
    samples.emplace_back("WZTo2L2Q"                       , vhbbPlot::kVZ     );
    samples.emplace_back("ZZTo2L2Nu"                      , vhbbPlot::kVZ     );
    samples.emplace_back("ZZTo2L2Q"                       , vhbbPlot::kVZ     );
    samples.emplace_back("ZZTo4Q"                         , vhbbPlot::kVZ     );
    samples.emplace_back("WWToLNuQQ"                      , vhbbPlot::kWW     );
    samples.emplace_back("TTbar_Powheg"                   , vhbbPlot::kTT     );
    samples.emplace_back("SingleTop_tT"                   , vhbbPlot::kTop    );
    samples.emplace_back("SingleTop_tTbar"                , vhbbPlot::kTop    );
    samples.emplace_back("SingleTop_tW"                   , vhbbPlot::kTop    );
    samples.emplace_back("SingleTop_tbarW"                , vhbbPlot::kTop    );
    samples.emplace_back("WJets_bHadrons_pt100to200"      , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_bHadrons_pt200toinf"      , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_bQuarks_pt100to200"       , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_bQuarks_pt200toinf"       , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_ht100to200"               , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_ht200to400"               , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_ht400to600"               , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_ht600to800"               , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_ht800to1200"              , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_ht1200to2500"             , vhbbPlot::kWjets  );
    samples.emplace_back("WJets_ht2500toinf"              , vhbbPlot::kWjets  );
    samples.emplace_back("ZJets_bHadrons_pt100to200"      , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_bHadrons_pt200toinf"      , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_bQuarks_pt100to200"       , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_bQuarks_pt200toinf"       , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_ht100to200"               , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_ht200to400"               , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_ht400to600"               , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_ht600to800"               , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_ht800to1200"              , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_ht1200to2500"             , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_ht2500toinf"              , vhbbPlot::kZjets  );
    samples.emplace_back("WmLNuHbb"                       , vhbbPlot::kWH     );
    samples.emplace_back("WpLNuHbb"                       , vhbbPlot::kWH     );
    samples.emplace_back("ZllHbb_mH125"                   , vhbbPlot::kZH     );
    samples.emplace_back("ggZllHbb_mH125"                 , vhbbPlot::kZH     );
  } else if(year==2017) {
    samples.emplace_back("SingleElectron"                 , vhbbPlot::kData   ); 
    samples.emplace_back("SingleMuon"                     , vhbbPlot::kData   );
    samples.emplace_back("Diboson_ww_CP5"                 , vhbbPlot::kWW     );
    samples.emplace_back("WZTo1L1Nu2Q"                    , vhbbPlot::kVZ     );
    samples.emplace_back("WZTo2L2Q"                       , vhbbPlot::kVZ     );
    samples.emplace_back("ZZTo2L2Q"                       , vhbbPlot::kVZ     );
    samples.emplace_back("QCD_ht1000to1500_CP5"           , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht1500to2000_CP5"           , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht2000toinf_CP5"            , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht200to300_CP5"             , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht300to500_CP5"             , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht500to700_CP5"             , vhbbPlot::kQCD    );
    samples.emplace_back("QCD_ht700to1000_CP5"            , vhbbPlot::kQCD    );
    samples.emplace_back("SingleTop_tT_CP5"               , vhbbPlot::kTop    );
    samples.emplace_back("SingleTop_tTbar_CP5"            , vhbbPlot::kTop    );
    samples.emplace_back("SingleTop_tW_CP5"               , vhbbPlot::kTop    );
    samples.emplace_back("SingleTop_tbarW_CP5"            , vhbbPlot::kTop    );
    samples.emplace_back("TTTo2L2Nu_CP5"                  , vhbbPlot::kTT     );
    samples.emplace_back("TTToSemiLeptonic_CP5"           , vhbbPlot::kTT     );
    samples.emplace_back("TTToHadronic_CP5"               , vhbbPlot::kTT     );
    samples.emplace_back("W1JetsToLNu_WpT100to150_CP5"    , vhbbPlot::kWjets  ); 
    samples.emplace_back("W1JetsToLNu_WpT150to250_CP5"    , vhbbPlot::kWjets  ); 
    samples.emplace_back("W1JetsToLNu_WpT250to400_CP5"    , vhbbPlot::kWjets  ); 
    samples.emplace_back("W1JetsToLNu_WpT400toinf_CP5"    , vhbbPlot::kWjets  ); 
    samples.emplace_back("W2JetsToLNu_WpT50to150_CP5"     , vhbbPlot::kWjets  ); 
    //samples.emplace_back("W2JetsToLNu_WpT100to150_CP5"    , vhbbPlot::kWjets  ); 
    //samples.emplace_back("W2JetsToLNu_WpT150to250_CP5"    , vhbbPlot::kWjets  ); 
    //samples.emplace_back("W2JetsToLNu_WpT250to400_CP5"    , vhbbPlot::kWjets  ); 
    samples.emplace_back("W2JetsToLNu_WpT150to250_CP5_notTruncated"    , vhbbPlot::kWjets  ); 
    samples.emplace_back("W2JetsToLNu_WpT250to400_CP5_notTruncated"    , vhbbPlot::kWjets  ); 
    samples.emplace_back("W2JetsToLNu_WpT400toinf_CP5"    , vhbbPlot::kWjets  ); 
    
    //samples.emplace_back("WJets_ht100to200_CP5"           , vhbbPlot::kWjets  ); 
    //samples.emplace_back("WJets_ht200to400_CP5"           , vhbbPlot::kWjets  ); 
    //samples.emplace_back("WJets_ht400to600_CP5"           , vhbbPlot::kWjets  ); 
    //samples.emplace_back("WJets_ht600to800_CP5"           , vhbbPlot::kWjets  ); 
    //samples.emplace_back("WJets_ht800to1200_CP5"          , vhbbPlot::kWjets  ); 
    //samples.emplace_back("WJets_ht1200to2500_CP5"         , vhbbPlot::kWjets  ); 

    samples.emplace_back("WmLNuHbb"                       , vhbbPlot::kWH     );
    samples.emplace_back("WpLNuHbb"                       , vhbbPlot::kWH     );
    samples.emplace_back("Z1Jets_ZpT150to250_CP5"         , vhbbPlot::kZjets  ); 
    samples.emplace_back("Z1Jets_ZpT250to400_CP5"         , vhbbPlot::kZjets  ); 
    samples.emplace_back("Z1Jets_ZpT400toinf_CP5"         , vhbbPlot::kZjets  ); 
    samples.emplace_back("Z1Jets_ZpT50to150_CP5"          , vhbbPlot::kZjets  ); 
    samples.emplace_back("Z2Jets_ZpT150to250_CP5"         , vhbbPlot::kZjets  ); 
    samples.emplace_back("Z2Jets_ZpT250to400_CP5"         , vhbbPlot::kZjets  ); 
    samples.emplace_back("Z2Jets_ZpT400toinf_CP5"         , vhbbPlot::kZjets  ); 
    samples.emplace_back("Z2Jets_ZpT50to150_CP5"          , vhbbPlot::kZjets  ); 
    samples.emplace_back("ZJets_m4_ht100to200_CP5"        , vhbbPlot::kZjets  ); 
    samples.emplace_back("ZJets_m4_ht200to400_CP5"        , vhbbPlot::kZjets  ); 
    samples.emplace_back("ZJets_m4_ht400to600_CP5"        , vhbbPlot::kZjets  ); 
    samples.emplace_back("ZJets_m4_ht600toinf_CP5"        , vhbbPlot::kZjets  ); 
    samples.emplace_back("ZJets_m4_ht70to100_CP5"         , vhbbPlot::kZjets  ); 
    samples.emplace_back("ZllHbb_mH125"                   , vhbbPlot::kZH     );
    samples.emplace_back("ggZllHbb_mH125"                 , vhbbPlot::kZH     );
  }
  if(multithread) std::random_shuffle(samples.begin(),samples.end());
  // End List of Samples
  /////////////////////////////
  
  // Load Shared Objects for ACLIC
  gSystem->Load("libPandaAnalysisFlat.so");
  
  assert(dataCardDir!="");
  if(dataCardDir!="") system(Form("mkdir -p MitVHBBAnalysis/datacards/%s",dataCardDir.Data()));
  ao.whichTriggers = (1<<pa::kSingleEleTrig) | (1<<pa::kSingleMuTrig);

 // Load Pileup Weights
  TString puPath =
    (year==2016)? "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root" :
    "MitAnalysisRunII/data/90x/puWeights_90x.root"; 
  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data())); assert(fPUFile);
  ao.puWeights     = (TH1D*)(fPUFile->Get("puWeights"    )); assert(ao.puWeights    );
  ao.puWeightsUp   = (TH1D*)(fPUFile->Get("puWeightsUp"  )); assert(ao.puWeightsUp  );
  ao.puWeightsDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(ao.puWeightsDown);
  ao.puWeights    ->SetDirectory(0); 
  ao.puWeightsUp  ->SetDirectory(0); 
  ao.puWeightsDown->SetDirectory(0); 
  delete fPUFile;
  // Done Loading Pileup Weights

  // Load JEC Uncertainties
  TString ak4JecUncPath = (year==2016)? 
    "PandaAnalysis/data/jec/23Sep2016V4/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt":
    "PandaAnalysis/data/jec/17Nov2017_V8//Fall17_17Nov2017_V8_MC_UncertaintySources_AK4PFchs.txt";
  TString ak8JecUncPath = (year==2016)? 
    "PandaAnalysis/data/jec/23Sep2016V4/Summer16_23Sep2016V4_MC_UncertaintySources_AK8PFPuppi.txt":
    "PandaAnalysis/data/jec/17Nov2017_V8//Fall17_17Nov2017_V8_MC_UncertaintySources_AK8PFPuppi.txt";
  ao.jecUncsAK4.reserve(NJES);
  ao.jecUncsAK8.reserve(NJES);
  for(unsigned iJES=1; iJES<NJES-1; iJES++) {
    if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
    TString shiftName(jesName(static_cast<shiftjes>(iJES)).Data());
    if(!shiftName.EndsWith("Up")) continue;
    TString jecUncName = shiftName(3,shiftName.Length()-5);//JES*Up
    JetCorrectionUncertainty *theAK4Unc = new JetCorrectionUncertainty(
      JetCorrectorParameters(ak4JecUncPath.Data(), jecUncName.Data()));
    JetCorrectionUncertainty *theAK8Unc = new JetCorrectionUncertainty(
      JetCorrectorParameters(ak8JecUncPath.Data(), jecUncName.Data()));
    ao.jecUncSourcesAK4[jecUncName] = theAK4Unc;
    ao.jecUncSourcesAK8[jecUncName] = theAK8Unc;
    ao.jecUncsAK4[iJES  ] = theAK4Unc;
    ao.jecUncsAK4[iJES+1] = theAK4Unc;
    ao.jecUncsAK8[iJES  ] = theAK8Unc;
    ao.jecUncsAK8[iJES+1] = theAK8Unc;
  }
  // Done Loading JEC Uncertainties

  // Choice of the MVA variable type, binning, and the name
  // This can be different for each control region
    // 1 - simple pT variable
  if(ao.MVAVarType==1) {
    if(selection==kWHSR || selection==kWHVZbbCR) {
      ao.MVAbins={100,120,140,160,180,200,250,300,350};
      ao.MVAVarName="H(bb) pT";
      ao.shapeType="ptShape";
    } else if(selection==kWHFJSR || selection==kWHVZbbFJCR) {
      ao.MVAbins={250,300,350,400,450,500,550,600};
      ao.MVAVarName="H(bb) pT";
      ao.shapeType="ptShape";
    }
  } else if(ao.MVAVarType==3) {
    if(selection==kWHSR || selection==kWHVZbbCR) {
      ao.MVAbins={-1,0,0.4,0.6,0.8,0.9,1};
      ao.MVAVarName="BDT Output";
      ao.shapeType="singleClassBDTShape"; 
    } else if(selection==kWHFJSR || selection==kWHVZbbFJCR) {
      ao.MVAbins={-0.8,-0.4,0,0.4,0.5,0.6};
      ao.MVAVarName="BDT Output";
      ao.shapeType="singleClassBDTShape"; 
    }
  }
  // 2 - multiclass BDT, not implemented
  if(selection==kWHLightFlavorCR) {
    if(year==2016) {
      ao.MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667};
      ao.MVAVarName="Subleading H(bb) CMVA";
      ao.shapeType="lesserCMVAShape";
    } else {
      ao.MVAbins={0, 0.1, 0.2, 0.3, 0.4, 0.5};
      ao.MVAVarName="Subleading H(bb) DeepCSV";
      ao.shapeType="lesserCSVShape";
    }
  } else if(selection==kWHHeavyFlavorLoMassCR || selection==kWHHeavyFlavorHiMassCR || selection==kWH2TopCR || selection==kWHPresel) {
    if(year==2016) {
      ao.MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000, 0.7333, 0.8667, 1.0000};
      ao.MVAVarName="Subleading H(bb) CMVA";
      ao.shapeType="lesserCMVAShape";
    } else {
      ao.MVAbins={0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1};
      ao.MVAVarName="Subleading H(bb) DeepCSV";
      ao.shapeType="lesserCSVShape";
    }
  } else if((selection>=kWHLightFlavorFJCR && selection<kWHFJSR) || selection==kWHFJPresel) {
    if(selection==kWHHeavyFlavorFJCR) {
      if(ao.vzbbMode) ao.MVAbins={40,45,50,55,60,65,70,75,80};
      else            ao.MVAbins={120,130,140,160,180,200};
    } else {
      ao.MVAbins={40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,160,180,200};
    }
    ao.MVAVarName="Fatjet soft drop mass [GeV]";
    ao.shapeType="softDropMassShape";
  }
  // Done choosing shape variable
  
  // Declare histograms for plotting
  printf("Building plotting histograms, please wait...\n");
  ao.xmin.resize(nPlots);
  ao.xmax.resize(nPlots);
  ao.nbins.resize(nPlots);
  ao.histoNames.resize(nPlots); 
  ao.histoTitles.resize(nPlots);
  { int p=0;
    ao.histoNames[p]="MVAVar"                  ; ao.histoTitles[p]=""                         ;                                                p++;
    ao.histoNames[p]="lepton1Pt"               ; ao.histoTitles[p]="Lepton p_{T} [GeV]"       ; ao.nbins[p]=  33; ao.xmin[p]=    20; ao.xmax[p]=   350; p++; 
    ao.histoNames[p]="lepton1Eta"              ; ao.histoTitles[p]="Lepton #eta"              ; ao.nbins[p]=  25; ao.xmin[p]=  -2.5; ao.xmax[p]=   2.5; p++; 
    ao.histoNames[p]="lepton1Charge"           ; ao.histoTitles[p]="Lepton charge"            ; ao.nbins[p]=   3; ao.xmin[p]=    -1; ao.xmax[p]=     2; p++; 
    ao.histoNames[p]="WBosonPhi"               ; ao.histoTitles[p]="W boson #eta"             ; ao.nbins[p]=  25; ao.xmin[p]=  -2.5; ao.xmax[p]=   2.5; p++; 
    ao.histoNames[p]="mT"                      ; ao.histoTitles[p]="W boson m_{T} [GeV]"      ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=  200.; p++;
    ao.histoNames[p]="pfmet"                   ; ao.histoTitles[p]="E_{T}^{miss} [GeV]"       ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=  500.; p++;
    ao.histoNames[p]="pfmetsig"                ; ao.histoTitles[p]="E_{T}^{miss} significance"; ao.nbins[p]=  30; ao.xmin[p]=     0; ao.xmax[p]=    30; p++; 
    ao.histoNames[p]="bdtValue"                ; ao.histoTitles[p]="BDT Output"               ; ao.nbins[p]=  40; ao.xmin[p]=    -1; ao.xmax[p]=    1.; p++; 
    if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
      // fatjet only plots
      ao.histoNames[p]="WBosonPt"           ; ao.histoTitles[p]="W boson pT [GeV]"      ; ao.nbins[p]=  55; ao.xmin[p]=   250; ao.xmax[p]=   700; p++; 
      ao.histoNames[p]="mSD"                ; ao.histoTitles[p]="Fatjet mSD [GeV]"      ; ao.nbins[p]=  32; ao.xmin[p]=    40; ao.xmax[p]=   200; p++; 
      ao.histoNames[p]="pTFJ"               ; ao.histoTitles[p]="Fatjet pT [GeV]"       ; ao.nbins[p]=  25; ao.xmin[p]=   250; ao.xmax[p]=   600; p++; 
      ao.histoNames[p]="Tau21SD"            ; ao.histoTitles[p]="#tau_{2}/#tau_{1} SD"  ; ao.nbins[p]=  20; ao.xmin[p]=     0; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="Tau32SD"            ; ao.histoTitles[p]="#tau_{3}/#tau_{2} SD"  ; ao.nbins[p]=  20; ao.xmin[p]=     0; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="doubleB"            ; ao.histoTitles[p]="Fatjet double b-tag"   ; ao.nbins[p]=  40; ao.xmin[p]=   -1.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="deltaEtaLep1FJ"     ; ao.histoTitles[p]="#Delta#eta(lep,FJ)"    ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="deltaPhiWHFJ"       ; ao.histoTitles[p]="#Delta#phi(W,FJ) [Rad]"; ao.nbins[p]=  20; ao.xmin[p]= 1.571; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="ptBalanceWHFJ"      ; ao.histoTitles[p]="|Fatjet pT / W pT|"    ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=    3.; p++; 
      ao.histoNames[p]="nIsojet"            ; ao.histoTitles[p]="N isojets"             ; ao.nbins[p]=   8; ao.xmin[p]=    0.; ao.xmax[p]=    8.; p++; 
      ao.histoNames[p]="isojetNBtags"       ; ao.histoTitles[p]="N isojet b-tags"       ; ao.nbins[p]=   4; ao.xmin[p]=    0.; ao.xmax[p]=    4.; p++; 
      ao.histoNames[p]="fjHTTFRec"          ; ao.histoTitles[p]="FJ HTT f_{rec}"        ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=   0.4; p++;
      ao.histoNames[p]="fjHTTMass"          ; ao.histoTitles[p]="FJ HTT mass [GeV]"     ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=  200.; p++;
      ao.histoNames[p]="psi022004031003"    ; ao.histoTitles[p]="#psi(2,2.0,4,3,1.0,3)" ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]= 0.05 ; p++;
      ao.histoNames[p]="psi022004022003"    ; ao.histoTitles[p]="#psi(2,2.0,4,2,2.0,3)" ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]= 0.04 ; p++;
      ao.histoNames[p]="psi022004030503"    ; ao.histoTitles[p]="#psi(2,2.0,4,3,0.5,3)" ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]= 2    ; p++;
      ao.histoNames[p]="dPhil1W"            ; ao.histoTitles[p]="#Delta#phi(lep,W)"       ; ao.nbins[p]=  32; ao.xmin[p]=    0.; ao.xmax[p]= 3.142; p++; 
    } else {
      ao.histoNames[p]="WBosonPt"           ; ao.histoTitles[p]="W boson pT [GeV]"        ; ao.nbins[p]=  50; ao.xmin[p]=     0; ao.xmax[p]=   500; p++; 
      ao.histoNames[p]="Mjj"                ; ao.histoTitles[p]="Dijet mass [GeV]"        ; ao.nbins[p]=  25; ao.xmin[p]=     0; ao.xmax[p]=   250; p++; 
      ao.histoNames[p]="pTjj"               ; ao.histoTitles[p]="Dijet pT [GeV]"          ; ao.nbins[p]=  18; ao.xmin[p]=    50; ao.xmax[p]=   350; p++; 
      ao.histoNames[p]="bjet1Pt"            ; ao.histoTitles[p]="B-jet 1 pT [GeV]"        ; ao.nbins[p]=  38; ao.xmin[p]=    20; ao.xmax[p]=   400; p++; 
      ao.histoNames[p]="bjet2Pt"            ; ao.histoTitles[p]="B-jet 2 pT [GeV]"        ; ao.nbins[p]=  38; ao.xmin[p]=    20; ao.xmax[p]=   400; p++; 
      ao.histoNames[p]="bjet1btag"          ; ao.histoTitles[p]="B-jet 1 btag"            ; ao.nbins[p]=  50; ao.xmin[p]=    0.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="bjet2btag"          ; ao.histoTitles[p]="B-jet 2 btag"            ; ao.nbins[p]=  50; ao.xmin[p]=    0.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="nJet"               ; ao.histoTitles[p]="N central AK4CHS jets"   ; ao.nbins[p]=   8; ao.xmin[p]=    0.; ao.xmax[p]=    8.; p++; 
      ao.histoNames[p]="deltaPhiWH"         ; ao.histoTitles[p]="#Delta#phi(W,H) [Rad]"   ; ao.nbins[p]=  20; ao.xmin[p]= 1.571; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="ptBalanceWH"        ; ao.histoTitles[p]="|H pT / W pT|"           ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=    3.; p++; 
      ao.histoNames[p]="sumEtSoft1"         ; ao.histoTitles[p]="#sum E_{T}(soft 1)"      ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=  300.; p++; 
      ao.histoNames[p]="nSoft2"             ; ao.histoTitles[p]="N^{soft}_{2}"            ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=   25.; p++; 
      ao.histoNames[p]="nSoft5"             ; ao.histoTitles[p]="N^{soft}_{5}"            ; ao.nbins[p]=  12; ao.xmin[p]=    0.; ao.xmax[p]=   12.; p++; 
      ao.histoNames[p]="nSoft10"            ; ao.histoTitles[p]="N^{soft}_{10}"           ; ao.nbins[p]=   8; ao.xmin[p]=    0.; ao.xmax[p]=    8.; p++; 
      ao.histoNames[p]="topWBosonCosThetaCS"; ao.histoTitles[p]="W boson cos #theta^{CS}" ; ao.nbins[p]=  40; ao.xmin[p]=   -1.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="hbbCosThetaJJ"      ; ao.histoTitles[p]="cos #theta(bb)"          ; ao.nbins[p]=  40; ao.xmin[p]=   -1.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="hbbCosThetaCSJ1"    ; ao.histoTitles[p]="cos #theta^{CS} bjet1"   ; ao.nbins[p]=  40; ao.xmin[p]=   -1.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="deltaEtaLep1H"      ; ao.histoTitles[p]="#Delta#eta(lep,jj)"      ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="dPhil1W"            ; ao.histoTitles[p]="#Delta#phi(lep,W)"       ; ao.nbins[p]=  32; ao.xmin[p]=    0.; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="dPhil1b1"           ; ao.histoTitles[p]="#Delta#phi(lep,b1)"      ; ao.nbins[p]=  32; ao.xmin[p]=    0.; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="dPhil1b2"           ; ao.histoTitles[p]="#Delta#phi(lep,b2)"      ; ao.nbins[p]=  32; ao.xmin[p]=    0.; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="dPhiWb1"            ; ao.histoTitles[p]="#Delta#phi(W,b1)"        ; ao.nbins[p]=  32; ao.xmin[p]=    0.; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="dPhiWb2"            ; ao.histoTitles[p]="#Delta#phi(W,b2)"        ; ao.nbins[p]=  32; ao.xmin[p]=    0.; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="dPhib1b2"           ; ao.histoTitles[p]="#Delta#phi(b1,b2)"       ; ao.nbins[p]=  32; ao.xmin[p]=    0.; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="dEtal1W"            ; ao.histoTitles[p]="|#Delta#eta(lep,W)|"     ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="dEtal1b1"           ; ao.histoTitles[p]="|#Delta#eta(lep,b1)|"    ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="dEtal1b2"           ; ao.histoTitles[p]="|#Delta#eta(lep,b2)|"    ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="dEtaWb1"            ; ao.histoTitles[p]="|#Delta#eta(W,b1)|"      ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="dEtaWb2"            ; ao.histoTitles[p]="|#Delta#eta(W,b2)|"      ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="dEtab1b2"           ; ao.histoTitles[p]="|#Delta#eta(b1,b2)|"     ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
    }
  }
  
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(int p=0; p<nPlots; p++) {
    if(ao.histoNames[p]=="") continue;
    for(unsigned char ic=0; ic<nPlotCategories; ic++) {
      if(ao.histoNames[p]=="MVAVar") {
        ao.histos[lep][p][ic] = new TH1F(
          Form("histo%d_W%sH%s_%s",
            ic,
            leptonStrings[lep].Data(),
            selectionNames[(int)selection].Data(),
            ao.histoNames[p].Data()
          ), ao.MVAVarName, ao.MVAbins.size()-1, ao.MVAbins.data());
      } else {
        ao.histos[lep][p][ic] = new TH1F(
          Form("histo%d_W%sH%s_%s",
            ic,
            leptonStrings[lep].Data(),
            selectionNames[(int)selection].Data(),
            ao.histoNames[p].Data()
          ), ao.histoTitles[p], ao.nbins[p], ao.xmin[p], ao.xmax[p]);
      }
      if(debug) printf("\tMade plot %s\n", ao.histos[lep][p][ic]->GetName());
      ao.histos[lep][p][ic]->Sumw2();
      ao.histos[lep][p][ic]->SetDirectory(0);
    }
  }

  printf("Done building plotting histograms\n");
  GeneralTree gt; 
  printf("Building uncertainty histograms, please wait...\n");
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(unsigned ic=kPlotData; ic!=nPlotCategories; ic++) {
    ao.histo_Baseline     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s"                  , plotBaseNames[ic].Data()));
    if(ic<kPlotVZbb) continue;
    ao.histo_pileupUp     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_pileupUp"         , plotBaseNames[ic].Data()));
    ao.histo_pileupDown   [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_pileupDown"       , plotBaseNames[ic].Data()));
    ao.histo_VHCorrUp     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VHCorrUp"         , plotBaseNames[ic].Data()));
    ao.histo_VHCorrDown   [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VHCorrDown"       , plotBaseNames[ic].Data()));
    ao.histo_QCDr1f2      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr1f2"          , plotBaseNames[ic].Data()));
    ao.histo_QCDr1f5      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr1f5"          , plotBaseNames[ic].Data()));
    ao.histo_QCDr2f1      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr2f1"          , plotBaseNames[ic].Data()));
    ao.histo_QCDr2f2      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr2f2"          , plotBaseNames[ic].Data()));
    ao.histo_QCDr5f1      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr5f1"          , plotBaseNames[ic].Data()));
    ao.histo_QCDr5f5      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr5f5"          , plotBaseNames[ic].Data()));
    ao.histo_QCDScaleUp   [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDScale%sUp"     , plotBaseNames[ic].Data(),plotBaseNames[ic].Data()));
    ao.histo_QCDScaleDown [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDScale%sDown"   , plotBaseNames[ic].Data(),plotBaseNames[ic].Data()));
    ao.histo_eleSFUp      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_eleSFUp"          , plotBaseNames[ic].Data()));
    ao.histo_eleSFDown    [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_eleSFDown"        , plotBaseNames[ic].Data()));
    ao.histo_muSFUp       [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_muSFUp"           , plotBaseNames[ic].Data()));
    ao.histo_muSFDown     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_muSFDown"         , plotBaseNames[ic].Data()));
    ao.histo_VGluUp       [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VjetsGluFracUp"   , plotBaseNames[ic].Data()));
    ao.histo_VGluDown     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VjetsGluFracDown" , plotBaseNames[ic].Data()));
    ao.histo_doubleBUp    [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_CMS_doubleBUp"    , plotBaseNames[ic].Data()));
    ao.histo_doubleBDown  [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_CMS_doubleBDown"  , plotBaseNames[ic].Data()));
    for(unsigned iJES=0; iJES<NJES; iJES++)
      ao.histo_jes[iJES][lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(
        Form("histo_%s_%s", plotBaseNames[ic].Data(), jesName(static_cast<shiftjes>(iJES)).Data())
      );
    for(unsigned iPt=0; iPt<5; iPt++)
    for(unsigned iEta=0; iEta<3; iEta++)
    for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
      GeneralTree::csvShift shift = gt.csvShifts[iShift];
      if (shift==GeneralTree::csvCent) continue;
      ao.histo_btag[iShift][iPt][iEta][lep][ic] =
        (TH1F*)ao.histos[lep][0][ic]->Clone(
          Form("histo_%s_CMS_VH_btag%d_pt%d_eta%d_%s",
          plotBaseNames[ic].Data(),ao.year,iPt,iEta,btagShiftName(shift))
        );
    }
  }
  printf("Done building uncertainty histograms\n");


  ////////////////////////////////////////////////////////////////////////
  // Load corrections to apply offline
  std::vector<std::string> btagSystNames;
  // get the official syst name strings
  btagSystNames.reserve(GeneralTree::nCsvShifts);
  for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
    GeneralTree::csvShift shift = gt.csvShifts[iShift];
    if (shift==GeneralTree::csvCent) continue;
    btagSystNames.push_back(GeneralTree::csvShiftName(shift).Data());
  }    
  ao.deepcsvSFs = new BTagCalibrationReader(
    BTagEntry::OP_RESHAPING,
    GeneralTree::csvShiftName(GeneralTree::csvCent).Data(), 
    btagSystNames);
  ao.deepcsvCalib = new BTagCalibration(
    "DeepCSV", 
    // no iterative fit, waiting on BTV
    //ao.year==2016?
    //  "PandaAnalysis/data/csv/DeepCSV_Moriond17_B_H.csv":
      "PandaAnalysis/data/csv/DeepCSV_94XSF_V2_B_F.csv"
  );
  ao.deepcsvSFs->load(*(ao.deepcsvCalib), BTagEntry::FLAV_B, "iterativeFit");
  ao.deepcsvSFs->load(*(ao.deepcsvCalib), BTagEntry::FLAV_C, "iterativeFit");
  ao.deepcsvSFs->load(*(ao.deepcsvCalib), BTagEntry::FLAV_UDSG, "iterativeFit");
  
  ao.ZjetsEWKCorr = EWKCorrPars(kZjets);
  ao.WjetsEWKCorr = EWKCorrPars(kWjets);
  if(!useHtBinnedVJetsKFactor) {
    TFile *kfactorsFile = TFile::Open("PandaAnalysis/data/higgs/hbb_kfactors.root","read");
    ao.kfactors_ZJets = (TH2D*) kfactorsFile->Get("h_ZJets");
    ao.kfactors_ZJets->SetDirectory(0);
    ao.kfactors_WJets = (TH2D*) kfactorsFile->Get("h_WJets");
    ao.kfactors_WJets->SetDirectory(0);
    kfactorsFile->Close();
  }

  // Done loading offline corrections
  ////////////////////////////////////////////////////////////////////////
  // Setup MVA training tree if applicable (Not implemented yet for boosted)
  if(selection==kWHSR||selection==kWHVZbbCR) {
    system(Form("mkdir -p MitVHBBAnalysis/mva/%s",dataCardDir.Data()));
    ao.mvaFile = new TFile(Form("MitVHBBAnalysis/mva/%s/WH%s_mvaTree%s.root",
      dataCardDir.Data(),selectionNames[selection].Data(),batchSuffix.Data()),"recreate");
    ao.mvaTree = new TTree("mvaTree","mvaTree");
    ao.mvaTree->Branch("weight"              , &ao.mva_weight             ); 
    ao.mvaTree->Branch("category"            , &ao.mva_category           ); 
    ao.mvaTree->Branch("eventNumber"         , &ao.mva_eventNumber        ); 
    ao.mvaTree->Branch("mT"                  , &ao.mva_mT                 ); 
    ao.mvaTree->Branch("dPhil1W"             , &ao.mva_dPhil1W            ); 
    ao.mvaTree->Branch("WBosonPt"            , &ao.mva_WBosonPt           ); 
    ao.mvaTree->Branch("dPhiLep1Met"         , &ao.mva_dPhiLep1Met        ); 
    ao.mvaTree->Branch("lepton1Pt"           , &ao.mva_lepton1Pt          ); 
    ao.mvaTree->Branch("lepton1Eta"          , &ao.mva_lepton1Eta         ); 
    ao.mvaTree->Branch("lepton1Charge"       , &ao.mva_lepton1Charge      ); 
    ao.mvaTree->Branch("pfmet"               , &ao.mva_pfmet              ); 
    
    ao.mvaTree->Branch("sumEtSoft1"          , &ao.mva_sumEtSoft1         ); 
    ao.mvaTree->Branch("nSoft2"              , &ao.mva_nSoft2             ); 
    ao.mvaTree->Branch("nSoft5"              , &ao.mva_nSoft5             ); 
    ao.mvaTree->Branch("nSoft10"             , &ao.mva_nSoft10            ); 
    ao.mvaTree->Branch("bjet1Pt"             , &ao.mva_bjet1Pt            ); 
    ao.mvaTree->Branch("bjet1btag"           , &ao.mva_bjet1btag          ); 
    ao.mvaTree->Branch("bjet2Pt"             , &ao.mva_bjet2Pt            ); 
    ao.mvaTree->Branch("bjet2btag"           , &ao.mva_bjet2btag          ); 
    ao.mvaTree->Branch("hbbpt"               , &ao.mva_hbbpt              ); 
    ao.mvaTree->Branch("hbbm"                , &ao.mva_hbbm               ); 
    ao.mvaTree->Branch("dPhiWH"              , &ao.mva_dPhiWH             ); 
    ao.mvaTree->Branch("ptBalanceWH"         , &ao.mva_ptBalanceWH        ); 
    ao.mvaTree->Branch("topMass"             , &ao.mva_topMass            ); 
    ao.mvaTree->Branch("dRBjets"             , &ao.mva_dRBjets            ); 
    ao.mvaTree->Branch("dEtaLep1H"           , &ao.mva_dEtaLep1H          ); 
    ao.mvaTree->Branch("nAddJet"             , &ao.mva_nAddJet            ); 
  } else if(selection==kWHFJSR||selection==kWHVZbbFJCR) {
    system(Form("mkdir -p MitVHBBAnalysis/mva/%s",dataCardDir.Data()));
    ao.mvaFile = new TFile(Form("MitVHBBAnalysis/mva/%s/WH%s_mvaTree%s.root",
      dataCardDir.Data(),selectionNames[selection].Data(),batchSuffix.Data()),"recreate");
    ao.mvaTree = new TTree("mvaTree","mvaTree");
    ao.mvaTree->Branch("weight"              , &ao.mva_weight             ); 
    ao.mvaTree->Branch("category"            , &ao.mva_category           ); 
    ao.mvaTree->Branch("eventNumber"         , &ao.mva_eventNumber        ); 
    ao.mvaTree->Branch("mT"                  , &ao.mva_mT                 ); 
    ao.mvaTree->Branch("dPhil1W"             , &ao.mva_dPhil1W            ); 
    ao.mvaTree->Branch("WBosonPt"            , &ao.mva_WBosonPt           ); 
    ao.mvaTree->Branch("dPhiLep1Met"         , &ao.mva_dPhiLep1Met        ); 
    ao.mvaTree->Branch("lepton1Pt"           , &ao.mva_lepton1Pt          ); 
    ao.mvaTree->Branch("lepton1Eta"          , &ao.mva_lepton1Eta         ); 
    ao.mvaTree->Branch("lepton1Charge"       , &ao.mva_lepton1Charge      ); 
    ao.mvaTree->Branch("pfmet"               , &ao.mva_pfmet              ); 
    ao.mvaTree->Branch("nIsojet"             , &ao.mva_nIsojet            ); 
    ao.mvaTree->Branch("MSD"                 , &ao.mva_MSD                ); 
    ao.mvaTree->Branch("Tau21SD"             , &ao.mva_Tau21SD            ); 
    ao.mvaTree->Branch("Tau32SD"             , &ao.mva_Tau32SD            ); 
    ao.mvaTree->Branch("fjPt"                , &ao.mva_fjPt               ); 
    ao.mvaTree->Branch("psi022004031003"     , &ao.mva_psi022004031003    ); 
    ao.mvaTree->Branch("psi022004022003"     , &ao.mva_psi022004022003    ); 
    ao.mvaTree->Branch("psi022004030503"     , &ao.mva_psi022004030503    ); 
    ao.mvaTree->Branch("ptBalanceWHFJ"       , &ao.mva_ptBalanceWHFJ      ); 
    ao.mvaTree->Branch("dEtaLep1FJ"          , &ao.mva_dEtaLep1FJ         ); 
    ao.mvaTree->Branch("dPhiWHFJ"            , &ao.mva_dPhiWHFJ           ); 
    ao.mvaTree->Branch("HTTFRec"             , &ao.mva_HTTFRec            ); 
  }
  
  // Instantiate TMVA reader
  if(MVAVarType>1) {
    TString bdtWeights="";
    if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) 
      if(ao.year==2016)
        bdtWeights = "MitVHBBAnalysis/weights/bdt_BDT_singleClass_boosted_wh2016.weights.xml";
      else
        bdtWeights = "MitVHBBAnalysis/weights/bdt_BDT_singleClass_boosted_wh2017.weights.xml";
    else
      if(ao.year==2016)
        bdtWeights = "MitVHBBAnalysis/weights/bdt_BDT_singleClass_resolved_wh2016.weights.xml";
      else
        bdtWeights = "MitVHBBAnalysis/weights/bdt_BDT_singleClass_resolved_wh2017.weights.xml";
    if(bdtWeights!="") for(unsigned nThread=0; nThread < (multithread? nThreads:1); nThread++) {
      TMVA::Reader *theReader = new TMVA::Reader("Silent");
      // This object is never deleted, which is a small memory leak,
      // but the TMVA Reader destructor has issues
      
      // Vars are hardcoded for now, could make it more general if we care
      if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) {
        theReader->AddVariable("dPhil1W"          , &ao.mvaInputs[nThread][ 0]);
        theReader->AddVariable("WBosonPt"         , &ao.mvaInputs[nThread][ 1]);
        theReader->AddVariable("lepton1Charge"    , &ao.mvaInputs[nThread][ 2]);
        theReader->AddVariable("nIsojet"          , &ao.mvaInputs[nThread][ 3]);
        theReader->AddVariable("MSD"              , &ao.mvaInputs[nThread][ 4]);
        theReader->AddVariable("Tau21SD"          , &ao.mvaInputs[nThread][ 5]);
        theReader->AddVariable("Tau32SD"          , &ao.mvaInputs[nThread][ 6]);
        theReader->AddVariable("fjPt"             , &ao.mvaInputs[nThread][ 7]);
        theReader->AddVariable("psi022004031003"  , &ao.mvaInputs[nThread][ 8]);
        theReader->AddVariable("psi022004022003"  , &ao.mvaInputs[nThread][ 9]);
        theReader->AddVariable("psi022004030503"  , &ao.mvaInputs[nThread][10]);
        theReader->AddVariable("ptBalanceWHFJ"    , &ao.mvaInputs[nThread][11]);
        theReader->AddVariable("dEtaLep1FJ"       , &ao.mvaInputs[nThread][12]);
        theReader->AddVariable("dPhiWHFJ"         , &ao.mvaInputs[nThread][13]);
        theReader->AddVariable("HTTFRec"          , &ao.mvaInputs[nThread][14]);
      } else {
        theReader->AddVariable("WBosonPt"         , &ao.mvaInputs[nThread][ 0]);
        theReader->AddVariable("dPhil1W"          , &ao.mvaInputs[nThread][ 1]);
        theReader->AddVariable("nSoft5"           , &ao.mvaInputs[nThread][ 2]);
        theReader->AddVariable("bjet1Pt"          , &ao.mvaInputs[nThread][ 3]);
        theReader->AddVariable("bjet2Pt"          , &ao.mvaInputs[nThread][ 4]);
        theReader->AddVariable("bjet2btag"        , &ao.mvaInputs[nThread][ 5]);
        theReader->AddVariable("hbbpt"            , &ao.mvaInputs[nThread][ 6]);
        theReader->AddVariable("hbbm"             , &ao.mvaInputs[nThread][ 7]);
        theReader->AddVariable("dPhiWH"           , &ao.mvaInputs[nThread][ 8]);
        theReader->AddVariable("ptBalanceWH"      , &ao.mvaInputs[nThread][ 9]);
        theReader->AddVariable("topMass"          , &ao.mvaInputs[nThread][10]);
        theReader->AddVariable("dRBjets"          , &ao.mvaInputs[nThread][11]);
        theReader->AddVariable("dEtaLep1H"        , &ao.mvaInputs[nThread][12]);
        theReader->AddVariable("nAddJet"          , &ao.mvaInputs[nThread][13]);
      }
      theReader->BookMVA("BDT", bdtWeights.Data());
      ao.reader.push_back(theReader);
    } else {
      printf("error with BDT weight configuration\n"); return;
    }
  }
  ////////////////////////////////////////////////////////////////////////
  // Begin Chain Loop
  vector<thread> threads;
  vector<TTree*> trees;
  vector<TFile*> files;
  for(auto const &sample: samples) {
    printf("### Opening sample %s ###\n",sample.first.Data());
    TString inputFileName = Form("%s/%s.root",ntupleDir.Data(),sample.first.Data());
    
    if(multithread) {
      // Spawn nThreads threads to process that many identical pointers to this file
      for(int split=0; split<nThreads; split++) {
        TFile *inputFile = TFile::Open(inputFileName,"READ"); assert(inputFile);
        files.push_back(inputFile);
        TTree *events = (TTree*)inputFile->Get("events"); assert(events);
        trees.push_back(events);
        threads.emplace_back(analyzeSample, sample, events, std::ref(ao), split);
      }
      for(auto &thread: threads)
        thread.join();
      threads.clear();
      for(auto file: files)
        file->Close();
      files.clear();
    } else {
      TFile *inputFile = TFile::Open(inputFileName,"READ"); assert(inputFile);
      TTree *events = (TTree*)inputFile->Get("events"); assert(events);
      analyzeSample(sample, events, std::ref(ao), -1);
      inputFile->Close();
    }
  } // End Chain Loop
  // Finish writing the MVA tree (if applicable)
  if(selection==kWHSR || selection==kWHVZbbCR || selection==kWHFJSR || selection==kWHVZbbFJCR) {
    ao.mvaFile->cd();
    ao.mvaTree->Write();
    ao.mvaFile->Close();
  }
  if(ao.deepcsvSFs) delete ao.deepcsvSFs;
  if(ao.deepcsvCalib) delete ao.deepcsvCalib;
  
  // Clean up JES
  for(unsigned iJES=0; iJES<NJES; iJES++) {
    if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
    TString shiftName(jesName(static_cast<shiftjes>(iJES)).Data());
    if(!shiftName.EndsWith("Up")) continue;
    delete ao.jecUncSourcesAK4[shiftName];
    delete ao.jecUncSourcesAK8[shiftName];
  }
 
  // renormalize V+Glu shapes to the nominal norm
  if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel)
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights()>0) {
      ao.histo_VGluUp  [lep][ic]->Scale(ao.histo_Baseline[lep][ic]->Integral()/ao.histo_VGluUp  [lep][ic]->Integral());
      ao.histo_VGluDown[lep][ic]->Scale(ao.histo_Baseline[lep][ic]->Integral()/ao.histo_VGluDown[lep][ic]->Integral());
    }
  }
 
  // Compute some uncertainties once event by events weights are filled
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
    for(int nb=1; nb<=ao.histo_Baseline[lep][ic]->GetNbinsX(); nb++){
      // compute QCD scale uncertainties bin-by-bin
      double diffQCDScale[6] = {
       TMath::Abs(ao.histo_QCDr1f2[lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb)),
       TMath::Abs(ao.histo_QCDr1f5[lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb)),
       TMath::Abs(ao.histo_QCDr2f1[lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb)),
       TMath::Abs(ao.histo_QCDr2f2[lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb)),
       TMath::Abs(ao.histo_QCDr5f1[lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb)),
       TMath::Abs(ao.histo_QCDr5f5[lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb))};

      double systQCDScale = diffQCDScale[0];
      for(int nqcd=0; nqcd<6; nqcd++) {
        if(diffQCDScale[nqcd] > systQCDScale) systQCDScale = diffQCDScale[nqcd];
      }

      if(ao.histo_Baseline[lep][ic]->GetBinContent(nb) > 0) 
        systQCDScale = 1.0+systQCDScale/ao.histo_Baseline[lep][ic]->GetBinContent(nb);
      else systQCDScale = 1;

      ao.histo_QCDScaleUp  [lep][ic]->SetBinContent(nb, ao.histo_Baseline[lep][ic]->GetBinContent(nb)*systQCDScale);
      ao.histo_QCDScaleDown[lep][ic]->SetBinContent(nb, ao.histo_Baseline[lep][ic]->GetBinContent(nb)/systQCDScale);
      
      // Force positive bin yields
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
        ao.histo_Baseline    [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_Baseline[lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_pileupUp    [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_pileupUp    [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_pileupDown  [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_pileupDown  [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_VHCorrUp    [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_VHCorrUp    [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_VHCorrDown  [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_VHCorrDown  [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_QCDScaleUp  [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_QCDScaleUp  [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_QCDScaleDown[lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_QCDScaleDown[lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_eleSFUp     [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_eleSFUp     [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_eleSFDown   [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_eleSFDown   [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_muSFUp      [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_muSFUp      [lep][ic]->GetBinContent(nb),1e-7f));
      ao.histo_muSFDown    [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_muSFDown    [lep][ic]->GetBinContent(nb),1e-7f));
      if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) {
        ao.histo_VGluUp       [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_VGluUp       [lep][ic]->GetBinContent(nb),1e-7f));
        ao.histo_VGluDown     [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_VGluDown     [lep][ic]->GetBinContent(nb),1e-7f));
        ao.histo_doubleBUp    [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_doubleBUp    [lep][ic]->GetBinContent(nb),1e-7f));
        ao.histo_doubleBDown  [lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_doubleBDown  [lep][ic]->GetBinContent(nb),1e-7f));
      }
      for(unsigned iJES=0; iJES<NJES; iJES++) { 
        if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
        ao.histo_jes[iJES][lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_jes[iJES][lep][ic]->GetBinContent(nb),1e-7f));
      }
      for(unsigned iPt=0; iPt<5; iPt++)
      for(unsigned iEta=0; iEta<3; iEta++)
      for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
        GeneralTree::csvShift shift = gt.csvShifts[iShift];
        if (shift==GeneralTree::csvCent) continue;
        ao.histo_btag[iShift][iPt][iEta][lep][ic]->SetBinContent(nb, TMath::Max((float)ao.histo_btag[iShift][iPt][iEta][lep][ic]->GetBinContent(nb),1e-7f));
      }
    } // all bins in histograms
    // Renormalize QCD scale uncertainties
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0 &&
       (ic!=kPlotVZbb&&ic!=kPlotVVLF&&ic!=kPlotTop)) {
      ao.histo_QCDScaleUp  [lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_QCDScaleUp  [lep][ic]->GetSumOfWeights());
      ao.histo_QCDScaleDown[lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_QCDScaleDown[lep][ic]->GetSumOfWeights());
    }
    // Renormalize lepton SF and pileup uncertainties
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0 &&
       (ic==kPlotTop||ic==kPlotTT||ic==kPlotZbb||ic==kPlotZb||ic==kPlotZLF)) {
      ao.histo_pileupUp  [lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_pileupUp  [lep][ic]->GetSumOfWeights());
      ao.histo_pileupDown[lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_pileupDown[lep][ic]->GetSumOfWeights());
      ao.histo_eleSFUp   [lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_eleSFUp   [lep][ic]->GetSumOfWeights());
      ao.histo_eleSFDown [lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_eleSFDown [lep][ic]->GetSumOfWeights());
      ao.histo_muSFUp	 [lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_muSFUp    [lep][ic]->GetSumOfWeights());
      ao.histo_muSFDown  [lep][ic]->Scale(ao.histo_Baseline[lep][ic]->GetSumOfWeights()/ao.histo_muSFDown  [lep][ic]->GetSumOfWeights());
    }
  } // all final states and categories

  // Write shape histograms to file
  for(unsigned lep=0; lep<nLepSel; lep++) {
    char outFileDatacardsName[200];
    sprintf(
      outFileDatacardsName,
      "MitVHBBAnalysis/datacards/%s/datacard_W%sH%s%s.root",
      dataCardDir.Data(),
      leptonStrings[lep].Data(),
      selectionNames[selection].Data(),
      batchSuffix.Data()
    );
    TFile* outFileDatacards = new TFile(outFileDatacardsName,"recreate");
    outFileDatacards->cd();

    for(unsigned ic=kPlotData; ic!=nPlotCategories; ic++) {
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() <= 0 && ic!=kPlotData) continue;
      ao.histo_Baseline    [lep][ic]->Write();
      if(ic<kPlotVZbb) continue;
      ao.histo_pileupUp    [lep][ic]->Write();
      ao.histo_pileupDown  [lep][ic]->Write();
      ao.histo_VHCorrUp    [lep][ic]->Write();
      ao.histo_VHCorrDown  [lep][ic]->Write();
      ao.histo_QCDScaleUp  [lep][ic]->Write();
      ao.histo_QCDScaleDown[lep][ic]->Write();
      ao.histo_eleSFUp     [lep][ic]->Write();
      ao.histo_eleSFDown   [lep][ic]->Write();
      ao.histo_muSFUp      [lep][ic]->Write();
      ao.histo_muSFDown    [lep][ic]->Write();
      for(unsigned iJES=0; iJES<NJES; iJES++) { 
        if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
        ao.histo_jes[iJES][lep][ic]->Write();
      }
      for(unsigned iPt=0; iPt<5; iPt++)
      for(unsigned iEta=0; iEta<3; iEta++)
      for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
        GeneralTree::csvShift shift = gt.csvShifts[iShift];
        if (shift==GeneralTree::csvCent) continue;
        ao.histo_btag[iShift][iPt][iEta][lep][ic]->Write();
      }
      if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) {
        ao.histo_VGluUp     [lep][ic]->Write();
        ao.histo_VGluDown   [lep][ic]->Write();
        ao.histo_doubleBUp  [lep][ic]->Write();
        ao.histo_doubleBDown[lep][ic]->Write();
      }
    }
    outFileDatacards->Close();
  }
  
  // Writing datacards
  if(!isBatchMode)
    writeDatacards(ao, dataCardDir);
  
  // Write plots
  char regionName[128];
  for(unsigned lep=0; lep<nLepSel; lep++) {
    sprintf(regionName, "W%sH%s",leptonStrings[lep].Data(),selectionNames[selection].Data());
    TString plotFileName = Form("MitVHBBAnalysis/datacards/%s/plots_%s%s.root",dataCardDir.Data(),regionName,batchSuffix.Data());
    TFile *outputPlots = new TFile(plotFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
    for(int p=0; p<nPlots; p++) {
      if(ao.histoNames[p]=="") continue;
      TDirectory *plotDir = outputPlots->mkdir(ao.histoNames[p]); plotDir->cd();
      for(unsigned ic=kPlotData; ic!=nPlotCategories; ic++) 
        ao.histos[lep][p][ic]->Write(Form("histo%d",ic));
      
    }
    outputPlots->Close();
  }


}

void analyzeSample(
  pair<TString,vhbbPlot::sampleType> sample,
  TTree *events,
  analysisObjects &ao,
  int split
) {
  TString sampleName = sample.first; sampleType type = sample.second;
  // Mass windows  
  float mjjLo=90, mjjHi=150, mSDLo=80, mSDHi=150;
  if(ao.vzbbMode) {
    mjjLo=60; mjjHi=120;
    mSDLo=50; mSDHi=120;
  }
  // Sample properties
  // Only use events with HT<100 for NLO pt binned samples in 2016
  bool isLowMassZjets = sampleName.Contains("ZJets_m10") || sampleName.Contains("ZJets_m4");
  bool isBQuarkEnriched = sampleName.Contains("bQuarks");
  bool isBHadronEnriched = sampleName.Contains("bHadrons");
  bool isInclusiveZjets = type==vhbbPlot::kZjets && !sampleName.Contains("_ht") && !sampleName.Contains("_pt");
  bool isV12jets = sampleName.Contains("Z1Jets") || sampleName.Contains("Z2Jets");
  bool isW2jets = sampleName.Contains("W2Jets"); 
  bool isNLOWjets = sampleName.Contains("W2Jets") || sampleName.Contains("W1Jets");
  bool isNLOZjets = sampleName.Contains("ZJets_pt") || sampleName.Contains("ZJets_m10") || isV12jets || sampleName=="ZJets_inclNLO_CP5"; 

  unsigned nThread = split>=0? split:0;
  // End sample properties
  ////////////////////////////////////////////////////////////////////////
  // Declare local analysis variables
  // Isojets: For boosted categories, the number of 30 GeV AK4 jets not in the fat jet
  vector<int> nIsojet(NJES);
  vector<vector<unsigned char>> isojets(NJES);
  vector<unsigned char> isojetNBtags(NJES);
  // Selection bits
  unsigned selectionBits[NJES], nMinusOneBits;
  unsigned char typeLepSel;
  unsigned char category; int countB;
  // Declare variables for the variations of the weights
  float weight;
  float weight_QCDr1f2 = 1, weight_QCDr1f5 = 1, weight_QCDr2f1 = 1, weight_QCDr2f2 = 1, weight_QCDr5f1 = 1, weight_QCDr5f5 = 1;
  float weight_muSF = 1, weight_elSF = 1;
  float weight_btag[GeneralTree::nCsvShifts][5][3];
  float weight_pileupUp = 1, weight_pileupDown = 1;
  float weight_VHCorrUp = 1, weight_VHCorrDown = 1;
  float weight_VGluUp = 1, weight_VGluDown = 1;
  float weight_doubleBUp = 1, weight_doubleBDown = 1;
  // Done declaring weight variables
  ////////////////////////////////////////////////////////////////////////
    
  // CMVA jet kinematic decorrelated weight nuisances
  vector<double> jetPts[5][3], jetPtsUp[5][3], jetPtsDown[5][3], jetEtas[5][3], jetBtags[5][3];
  vector<int> jetFlavors[5][3];
  for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
    jetPts    [iPt][iEta].reserve(20);
    jetPtsUp  [iPt][iEta].reserve(20);
    jetPtsDown[iPt][iEta].reserve(20);
    jetEtas   [iPt][iEta].reserve(20);
    jetBtags  [iPt][iEta].reserve(20);
    jetFlavors[iPt][iEta].reserve(20);
  }
  ////////////////////////////////////////////////////////////////////////
  // Instantiate GeneralTree object for reading this ntuple
  GeneralTree gt;
  gt.is_hbb         = true;
  gt.is_fatjet      = true;
  gt.is_leptonic    = true;
  gt.btagWeights    = true;
  gt.is_breg        = false;
  gt.useCMVA        = false;
  // Branches not in GeneralTree;
  std::map<TString, void*> extraAddresses;
  float normalizedWeight; unsigned char npnlo;
  extraAddresses["normalizedWeight"] = &normalizedWeight;
    
  // Map of the branches
  std::map<TString, TBranch*> b;
  // Dummy tree for a list of branch names
  TTree *dummyTree = new TTree("dummyTree","dummyTree");
  gt.WriteTree(dummyTree);
  // Done making GeneralTree object
  ////////////////////////////////////////////////////////////////////////

  // Nasty hack code to get the tree branches and their addresses, have to do it for each file
  TObjArray *listOfBranches = events->GetListOfBranches();
  for(unsigned iB=0; iB<(unsigned)listOfBranches->GetEntries(); iB++) {
    TBranch *branch = (TBranch*)listOfBranches->At(iB);
    TString branchName = branch->GetName();
    if(type==kData) { 
      if(branchName.Contains("_JES")) continue;
    }
    bool isExtraBranch=false;
    auto x = extraAddresses.find(branchName);
    isExtraBranch= (x!=extraAddresses.end());
    if(isExtraBranch) {
      b[x->first] = events->GetBranch(x->first);
      if(!b[x->first]) { throw std::runtime_error(Form("Extra branch %s could not be found in events tree\n", x->first.Data())); }
      b[x->first]->SetAddress(x->second);
      if(ao.debug>=2) printf("\tBooking extra branch \"%s\"\n", x->first.Data());
      continue;
    }
    TBranch *dummyBranch = dummyTree->GetBranch(branchName);
    if(!dummyBranch) { 
      printf("WARNING: Couldn't find branch \"%s\" in the dummyTree, skipping (You may be using ntuples produced with an older version of GeneralTree)\n", branchName.Data()); 
      continue;
    }
    void* address=dummyBranch->GetAddress();
    b[branchName] = events->GetBranch(branchName);
    if(!b[branchName]) { throw std::runtime_error(Form("Branch \"%s\" could not be found in events tree", x->first.Data())); }
    b[branchName]->SetAddress(address);
    if(ao.debug>=2) printf("\tBooking GeneralTree branch \"%s\"\n", branchName.Data());
  }
  delete dummyTree; // no longer needed
  // End Book Branches

  Long64_t nentries = events->GetEntries();
  //nentries = TMath::Min(Long64_t(1e5),nentries);
  // Begin Event Loop
  for (Long64_t ientry=0; ientry<nentries; ientry++) {
    if(ao.debug && ientry!=0) usleep(2e5);
    if(ao.debug || ientry%100000==0) printf("> Reading entry %lld/%lld of %s (thread #%d)...\n",ientry,nentries, sampleName.Data(), split);
    if(split!=-1 && (ientry%nThreads)!=split) continue;
    
    // Stitching Cuts/Weights
    float stitchWeight=1;
    if(ao.year==2016) {

      // B-enriched sample stitching
      // Let the b-enriched samples eat 90% of the XS where there are b quarks or status 2 b hadrons
      if(type==kWjets||type==kZjets) {
        bLoad(b["nStatus2BHadrons"],ientry); // number of B hadrons at matrix element level
        bLoad(b["nB"],ientry); // number of B quarks
        bLoad(b["trueGenBosonPt"],ientry);
        if(isInclusiveZjets && (isBQuarkEnriched||isBHadronEnriched) && gt.trueGenBosonPt>=100) continue;
        bool hasBQuarks  = gt.nB > 0;
        bool hasBHadrons = gt.nStatus2BHadrons>0 && gt.nB==0;
        // Orthogonalize        
        if(isBQuarkEnriched && !hasBQuarks) continue;
        if(isBHadronEnriched && !hasBHadrons) continue;
        // Downweight
        if(gt.trueGenBosonPt<100) {
           stitchWeight = 1;
        } else if(hasBQuarks) {
          if(isBQuarkEnriched) stitchWeight = 0.9;
          else                 stitchWeight = 0.1;
        } else if(hasBHadrons) {
          if(isBHadronEnriched) stitchWeight = 0.9;
          else                  stitchWeight = 0.1;
        }
      }
    } else if(ao.year==2017) {
      //if(isW2jets) {
      //  // Handle the overlap of the samples W2JetsToLNu_WpT100to150_CP5, W2JetsToLNu_WpT50to150_CP5
      //  bLoad(b["trueGenBosonPt"],ientry);
      //  if(gt.trueGenBosonPt>=100 && gt.trueGenBosonPt<150)
      //    stitchWeight = 0.5;
      //}
    }

    //////////////////////
    // Clear variables
    typeLepSel=99; // 0: W(mn), 1: W(en)

    // Clear the jet pts for the b tag decorrelation
    for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
      jetPts    [iPt][iEta].clear();
      jetPtsUp  [iPt][iEta].clear();
      jetPtsDown[iPt][iEta].clear();
      jetEtas   [iPt][iEta].clear();
      jetBtags  [iPt][iEta].clear();
      jetFlavors[iPt][iEta].clear();
    }
    for(unsigned iJES=0; iJES<NJES; iJES++) {
      nIsojet[iJES]=0;
      isojetNBtags[iJES]=0;
      isojets[iJES].clear();
    }
    // Done clearing variables
    //////////////////////
    
    // Basic cuts
    bLoad(b["nLooseLep"],ientry);
    if(gt.nLooseLep!=1) continue; 
    if(ao.debug) printf("  Passed lepton multiplicity\n");

    // Trigger
    if(type==kData) {
      bLoad(b["trigger"],ientry);
      bool passTrigger = (gt.trigger & ao.whichTriggers) !=0;
      if(!passTrigger) continue;
      if(ao.debug) printf("  Passed trigger\n");
    }
    bLoad(b["metFilter"],ientry);
    if(gt.metFilter!=1) continue;
    if(ao.debug) printf("  Passed MET filters\n");

    // Lepton ID and isolation
    bLoad(b["nLooseElectron"],ientry);
    bLoad(b["nTightElectron"],ientry);
    bLoad(b["nLooseMuon"],ientry);
    bLoad(b["nTightMuon"],ientry);
    bLoad(b["muonSelBit"],ientry);
    bLoad(b["muonPt"],ientry);
    bLoad(b["muonEta"],ientry);
    bLoad(b["muonPdgId"],ientry);
    bLoad(b["muonCombIso"],ientry);
    bLoad(b["electronSelBit"],ientry);
    bLoad(b["electronPt"],ientry);
    bLoad(b["electronEta"],ientry);
    bLoad(b["electronPdgId"],ientry);
    bLoad(b["electronCombIso"],ientry);
    bLoad(b["muonD0"],ientry);
    bLoad(b["muonDZ"],ientry);
    bLoad(b["electronD0"],ientry);
    bLoad(b["electronDZ"],ientry);

    if(
      gt.nLooseMuon==1                       && 
      gt.muonPt[0]>25                        && 
      (gt.muonSelBit[0] & pa::kTight)!=0     &&
      gt.muonCombIso[0]/gt.muonPt[0] < 0.06  &&
      gt.muonD0[0]<0.20 && gt.muonDZ[0]<0.50
    ) typeLepSel=0;
    else if(
      gt.nLooseElectron==1 &&
      gt.electronPt[0]>30                           && 
      (gt.electronSelBit[0] & pa::kEleMvaWP80)!=0   &&
      gt.electronCombIso[0]/gt.electronPt[0] < 0.06 &&
      (
        (fabs(gt.electronEta[0])<1.479  && gt.electronD0[0]<0.05 && gt.electronDZ[0]<0.10) || 
        (fabs(gt.electronEta[0])>=1.479 && gt.electronD0[0]<0.10 && gt.electronDZ[0]<0.20)
      )
    ) typeLepSel=1;
    if(typeLepSel!=0 && typeLepSel!=1) continue;
    if(ao.debug) printf("  Passed lepton ID/iso multiplicity\n");
    
    // Z boson basics
    bLoad(b["topWBosonPt"],ientry);
    if(gt.topWBosonPt<90) continue;
    if(ao.debug) printf("  Passed W boson reconstruction\n");

    // Lepton kinematics
    float lepton1Pt,lepton1Eta,lepton1Phi,lepton1RelIso,lepton1D0,lepton1DZ;
    int lepton1Charge;
    if (typeLepSel==0) {
      bLoad(b["muonEta"],ientry);
      bLoad(b["muonPhi"],ientry);
      bLoad(b["muonCombIso"],ientry);
      bLoad(b["muonPdgId"],ientry);
      lepton1Pt     = gt.muonPt[0];
      lepton1Eta    = gt.muonEta[0];
      lepton1Phi    = gt.muonPhi[0]; 
      lepton1D0     = gt.muonD0[0];
      lepton1DZ     = gt.muonDZ[0];
      lepton1Charge = gt.muonPdgId[0]>0? -1:1;
      lepton1RelIso = gt.muonCombIso[0]/gt.muonPt[0];
    } else if(typeLepSel==1) {
      bLoad(b["electronEta"],ientry);
      bLoad(b["electronPhi"],ientry);
      bLoad(b["electronCombIso"],ientry);
      bLoad(b["electronPdgId"],ientry);
      lepton1Pt     = gt.electronPt[0]; 
      lepton1Eta    = gt.electronEta[0];
      lepton1Phi    = gt.electronPhi[0];
      lepton1D0     = gt.electronD0[0];
      lepton1DZ     = gt.electronDZ[0];
      lepton1Charge = gt.electronPdgId[0]>0? -1:1;
      lepton1RelIso = gt.electronCombIso[0]/gt.electronPt[0];
    }
    if(ao.debug) printf("  Passed lepton kinematics\n");

    // Jet multiplicity
    bool isBoostedCategory=false;
    float deltaPhiWHFJ=-1;
    if(ao.useBoostedCategory) { 
      bLoad(b["fjMSD_corr"],ientry);
      bLoad(b["fjPt"],ientry);
      bLoad(b["fjEta"],ientry);
      bLoad(b["fjPhi"],ientry);
      bLoad(b["topWBosonPt"],ientry);
      bLoad(b["topWBosonPhi"],ientry);
      if(gt.fjPt[0]>0) // protection against NaN values of fjPhi
        deltaPhiWHFJ = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.fjPhi));
      if(
        gt.fjPt[0] >= 250 && 
        gt.fjMSD_corr[0] >= 40 &&
        fabs(gt.fjEta) < 2.4 &&
        gt.topWBosonPt >= 250 &&
        deltaPhiWHFJ >= 2.5
      ) isBoostedCategory=true;
      // If we consider a boosted category splitting,
      // only put boosted (resolved) events in boosted (resolved) regions
      //printf("isBoostedCategory=%d (fjPt=%.1f, fjMSD=%.1f, topWBosonPt=%.1f, deltaPhiWHFJ=%.2f)\n", isBoostedCategory,gt.fjPt[0], gt.fjMSD[0], gt.topWBosonPt, deltaPhiWHFJ);
      if(isBoostedCategory ^ (ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel))
        continue;
    }
    // Category Assignment for Plotting and Datacards
    if(type!=kData) {
      if(isBoostedCategory) {
        bLoad(b["fjGenNumB"],ientry);
        countB = gt.fjGenNumB;
      } else {
        bLoad(b["nBGenJets"],ientry);
        countB = gt.nBGenJets;
      }
    }
    if     (type==vhbbPlot::kData ) category=kPlotData ;
    else if(type==vhbbPlot::kQCD  ) category=kPlotQCD  ;
    else if(type==vhbbPlot::kWW   ) category=kPlotVVLF ;
    else if(type==vhbbPlot::kTT   ) category=kPlotTT   ; 
    else if(type==vhbbPlot::kTop  ) category=kPlotTop  ;
    else if(type==vhbbPlot::kWH   ) category=kPlotWH   ;
    else if(type==vhbbPlot::kZH   ) {
      if(sampleName.Contains("gg")) 
        category=kPlotGGZH;
      else
        category=kPlotZH;
    } else if(type==vhbbPlot::kWjets) {
      if(countB>1) category=kPlotWbb;
      else if(countB>0) category=kPlotWb;
      else category=kPlotWLF;
    } else if(type==vhbbPlot::kZjets) {
      if(countB>1) category=kPlotZbb;
      else if(countB>0) category=kPlotZb;
      else category=kPlotZLF;
    } else if(type==kVZ) {
      if(countB>1) category=kPlotVZbb;
      else category=kPlotVVLF;
    } else throw std::runtime_error("category problem!");

    // Jet kinematics
    if(isBoostedCategory) {
      // No checks here? 
    } else { 
      bLoad(b["jotPt"],ientry);
      bLoad(b["nJet"],ientry);
      bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
      bLoad(b["hbbpt"],ientry);
      bLoad(b["hbbm_reg"],ientry);
      bLoad(b["jotPt"],ientry);
      if(
        gt.nJet[0]<2 || 
        gt.hbbpt[0]<50 || 
        gt.hbbm_reg[0]<0 ||
        gt.hbbm_reg[0]>250 ||
        gt.jotPt[0][gt.hbbjtidx[0][0]]<25 ||
        gt.jotPt[0][gt.hbbjtidx[0][1]]<25   
      ) continue;
    }
    if(ao.debug) printf("  Passed jet kinematics\n");
    if(ao.debug) printf("Passed preselection!\n");
    
    // Met
    bLoad(b["pfmet"],ientry);
    //bLoad(b["pfmetphi"],ientry);
    
    // Jet B-tagging
    bool bjet1IsLoose=false, bjet1IsMedium=false, bjet1IsTight=false, bjet2IsLoose=false;
    float bjet1btag=-2, bjet2btag=-2;
    bLoad(b["jotCMVA"],ientry);
    bLoad(b["jotCSV"],ientry);
    if(ao.year==2016) {
      bjet1btag = TMath::Max(gt.jotCSV[gt.hbbjtidx[0][0]],gt.jotCSV[gt.hbbjtidx[0][1]]);
      bjet2btag = TMath::Min(gt.jotCSV[gt.hbbjtidx[0][0]],gt.jotCSV[gt.hbbjtidx[0][1]]);
      bjet1IsLoose  = bjet1btag > deepcsv16Loose;
      bjet1IsMedium = bjet1btag > deepcsv16Medium;
      bjet1IsTight  = bjet1btag > deepcsv16Tight;
      bjet2IsLoose  = bjet2btag > deepcsv16Loose;
    } else if(ao.year==2017) {
      bjet1btag = TMath::Max(gt.jotCSV[gt.hbbjtidx[0][0]],gt.jotCSV[gt.hbbjtidx[0][1]]);
      bjet2btag = TMath::Min(gt.jotCSV[gt.hbbjtidx[0][0]],gt.jotCSV[gt.hbbjtidx[0][1]]);
      bjet1IsLoose  = bjet1btag > deepcsvLoose;
      bjet1IsMedium = bjet1btag > deepcsvMedium;
      bjet1IsTight  = bjet1btag > deepcsvTight;
      bjet2IsLoose  = bjet2btag > deepcsvLoose;
    }

    bLoad(b["nJot"],ientry);
    bLoad(b["nJotMax"],ientry);
    bLoad(b["jotPt"],ientry);
    bLoad(b["jotEta"],ientry);
    bLoad(b["jotPhi"],ientry);
    bLoad(b["jotM"],ientry);
    bLoad(b["jotFlav"],ientry);
    float bjet1Pt  = gt.jotPt[0][gt.hbbjtidx[0][0]];
    float bjet2Pt  = gt.jotPt[0][gt.hbbjtidx[0][1]];
    float bjet1Eta = gt.jotEta[gt.hbbjtidx[0][0]];
    float bjet2Eta = gt.jotEta[gt.hbbjtidx[0][1]];
    float bjet1Phi = gt.jotPhi[gt.hbbjtidx[0][0]];
    float bjet2Phi = gt.jotPhi[gt.hbbjtidx[0][1]];
    
    if(isBoostedCategory) {
      // Isojets for boosted category
      bLoad(b["fjEta"],ientry);
      bLoad(b["fjPhi"],ientry);
      for(unsigned iJES=0; iJES<NJES; iJES++) 
      for(unsigned char iJ=0; iJ<gt.nJotMax; iJ++) {
        if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
        if(fabs(gt.jotEta[iJ])>2.4) continue;
        float dR2JetFatjet=pow(gt.jotEta[iJ]-gt.fjEta,2)+pow(TVector2::Phi_mpi_pi(gt.jotPhi[iJ]-gt.fjPhi),2);
        if(dR2JetFatjet<0.64) continue;
        
        float isojetBtag = (ao.year==2016)? gt.jotCMVA[iJ] : gt.jotCSV[iJ];
        if(iJES!=0) {
          jecAk4UncMutex.lock();
          bool isUp = !(iJES%2==0);
          ao.jecUncsAK4[iJES]->setJetPt (gt.jotPt[0][iJ]);
          ao.jecUncsAK4[iJES]->setJetEta(gt.jotEta  [iJ]);
          float relUnc = ao.jecUncsAK4[iJES]->getUncertainty(isUp);
          jecAk4UncMutex.unlock();
          if(!isUp) relUnc*=-1;
          gt.jotPt[iJES][iJ] = gt.jotPt[0][iJ]*(1+relUnc);
        }
        if(gt.jotPt[iJES][iJ]>30) {
          isojets[iJES].push_back(iJ);
          if(isojetBtag>ao.isojetBtagCut) isojetNBtags[iJES]++;
        }
        
        // Decorrelated b-tag nuisances for the isojets
        // Nominal JES scenario only
        if(iJES!=0) continue;
        int iPt=-1, iEta=-1;
        double jetAbsEta=fabs(gt.jotEta[iJ]);
        //if      (gt.jotPt[0][iJ] >= 19.99 && gt.jotPt[0][iJ] < 30 ) iPt = 0;
        //else if (gt.jotPt[0][iJ] >= 30    && gt.jotPt[0][iJ] < 40 ) iPt = 1;
        if      (gt.jotPt[0][iJ] >= 30    && gt.jotPt[0][iJ] < 40 ) iPt = 1;
        else if (gt.jotPt[0][iJ] >= 40    && gt.jotPt[0][iJ] < 60 ) iPt = 2;
        else if (gt.jotPt[0][iJ] >= 60    && gt.jotPt[0][iJ] < 100) iPt = 3;
        else if (gt.jotPt[0][iJ] >= 100                           ) iPt = 4;
        if      (jetAbsEta >= 0   && jetAbsEta < 0.8  ) iEta = 0;
        else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6  ) iEta = 1;
        else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41 ) iEta = 2;
        
        if(iPt>=0 && iEta>=0) {
          jetPts    [iPt][iEta].push_back(gt.jotPt[iJES][iJ]);
          jetEtas   [iPt][iEta].push_back(gt.jotEta[iJ]);
          jetFlavors[iPt][iEta].push_back(gt.jotFlav[iJ]);
          jetBtags  [iPt][iEta].push_back(isojetBtag);
        }
      }
      for(unsigned iJES=0; iJES<NJES; iJES++)
        nIsojet[iJES] = isojets[iJES].size();
    } else {
      // In resolved case, need decorrelated b-tag nuisances for all the ak4 jets
      // Nominal JES scenario only
      for(unsigned char iJ=0; iJ<gt.nJot[0]; iJ++) {
        int iPt=-1, iEta=-1;
        double jetAbsEta=fabs(gt.jotEta[iJ]);
        if      (gt.jotPt[0][iJ] >= 19.99 && gt.jotPt[0][iJ] < 30 ) iPt = 0;
        else if (gt.jotPt[0][iJ] >= 30    && gt.jotPt[0][iJ] < 40 ) iPt = 1;
        else if (gt.jotPt[0][iJ] >= 40    && gt.jotPt[0][iJ] < 60 ) iPt = 2;
        else if (gt.jotPt[0][iJ] >= 60    && gt.jotPt[0][iJ] < 100) iPt = 3;
        else if (gt.jotPt[0][iJ] >= 100                           ) iPt = 4;
        if      (jetAbsEta >= 0   && jetAbsEta < 0.8  ) iEta = 0;
        else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6  ) iEta = 1;
        else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41 ) iEta = 2;
        float btag = gt.jotCSV[iJ];
        if(ao.debug>=3) printf("jet with (pt,|eta|,flav,btag)=(%.2f,%.3f,%d,%.4f) => (iPt,iEta) = (%d,%d)\n",gt.jotPt[0][iJ],jetAbsEta,gt.jotFlav[iJ],btag,iPt,iEta);
        if(iPt>=0 && iEta>=0) {
          jetPts    [iPt][iEta].push_back(gt.jotPt[0][iJ]);
          jetEtas   [iPt][iEta].push_back(gt.jotEta[iJ]);
          jetFlavors[iPt][iEta].push_back(gt.jotFlav[iJ]);
          jetBtags  [iPt][iEta].push_back(btag);
        }
      }
    }
    
    // Load branches and calculate stuff for the cuts
    bLoad(b["eventNumber"],ientry);
    bLoad(b["topWBosonPhi"],ientry);
    if(isBoostedCategory) {
      bLoad(b["fjPt"],ientry);
      bLoad(b["fjPhi"],ientry);
      bLoad(b["fjEta"],ientry);
      bLoad(b["fjMSD_corr"],ientry);
      bLoad(b["fjDoubleCSV"],ientry);
      bLoad(b["fjTau21SD"],ientry);
      bLoad(b["fjTau32SD"],ientry);
      if(type!=vhbbPlot::kData) for(unsigned iJES=1; iJES<NJES; iJES++) {
        if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
        jecAk8UncMutex.lock();
        bool isUp = !(iJES%2==0);
        ao.jecUncsAK8[iJES]->setJetPt (gt.fjPt[0] );
        ao.jecUncsAK8[iJES]->setJetEta(gt.fjEta   );
        float relUnc = ao.jecUncsAK8[iJES]->getUncertainty(isUp);
        jecAk8UncMutex.unlock();
        if(!isUp) relUnc*=-1;
        gt.fjPt[iJES] = gt.fjPt[0]*(1+relUnc);
        gt.fjMSD_corr[iJES] = gt.fjPt[0]*(1+relUnc);
      }

    } else {
      bLoad(b["hbbpt"],ientry);
      bLoad(b["hbbpt_reg"],ientry);
      bLoad(b["hbbeta"],ientry);
      bLoad(b["hbbphi"],ientry);
      bLoad(b["hbbm_reg"],ientry);
      bLoad(b["hbbm"],ientry);
      bLoad(b["pfmet"],ientry);
      bLoad(b["pfmetsig"],ientry);
      bLoad(b["jotM"],ientry);
      if(type!=vhbbPlot::kData) for(unsigned iJES=1; iJES<NJES; iJES++) {
        if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
        jecAk4UncMutex.lock();
        bool isUp = !(iJES%2==0);
        ao.jecUncsAK4[iJES]->setJetPt (gt.jotPt [0][gt.hbbjtidx[0][0]]);
        ao.jecUncsAK4[iJES]->setJetEta(gt.jotEta   [gt.hbbjtidx[0][0]]);
        float relUnc1 = ao.jecUncsAK4[iJES]->getUncertainty(isUp);
        ao.jecUncsAK4[iJES]->setJetPt (gt.jotPt [0][gt.hbbjtidx[0][1]]);
        ao.jecUncsAK4[iJES]->setJetEta(gt.jotEta   [gt.hbbjtidx[0][1]]);
        float relUnc2 = ao.jecUncsAK4[iJES]->getUncertainty(isUp);
        jecAk4UncMutex.unlock();
        if(!isUp) { relUnc1*=-1; relUnc2*=-1; }
        gt.jotPt[iJES][gt.hbbjtidx[0][0]] = gt.jotPt[0][gt.hbbjtidx[0][0]]*(1+relUnc1);
        gt.jotPt[iJES][gt.hbbjtidx[0][1]] = gt.jotPt[0][gt.hbbjtidx[0][1]]*(1+relUnc2);
        TLorentzVector hbbjt1,hbbjt2,hbbsystem;
        hbbjt1.SetPtEtaPhiM(
          gt.jotPt [iJES][gt.hbbjtidx[0][0]],
          gt.jotEta      [gt.hbbjtidx[0][0]],
          gt.jotPhi      [gt.hbbjtidx[0][0]],
          gt.jotM        [gt.hbbjtidx[0][0]]);
        hbbjt2.SetPtEtaPhiM(
          gt.jotPt [iJES][gt.hbbjtidx[0][1]],
          gt.jotEta      [gt.hbbjtidx[0][1]],
          gt.jotPhi      [gt.hbbjtidx[0][1]],
          gt.jotM        [gt.hbbjtidx[0][1]]);
        hbbsystem=hbbjt1+hbbjt2;
        // Assume the regression is conformal...
        gt.hbbpt_reg[iJES] = gt.hbbpt_reg[0] * hbbsystem.Pt()/gt.hbbpt[0];
        gt.hbbm_reg[iJES]  = gt.hbbm_reg[0]  * hbbsystem.M() /gt.hbbm[0];
        gt.hbbphi[iJES]    = hbbsystem.Phi();
      }
      // Handle the NJET variations
      for(unsigned iJES=1; iJES<NJES; iJES++) {
        if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
        gt.nJet[iJES]=0;
        gt.nJot[iJES]=0;
        for(unsigned char iJ=0; iJ<gt.nJotMax; iJ++) {
          jecAk4UncMutex.lock();
          bool isUp = !(iJES%2==0);
          ao.jecUncsAK4[iJES]->setJetPt (gt.jotPt[0][iJ]);
          ao.jecUncsAK4[iJES]->setJetEta(gt.jotEta  [iJ]);
          float relUnc = ao.jecUncsAK4[iJES]->getUncertainty(isUp);
          jecAk4UncMutex.unlock();
          if(!isUp) relUnc*=-1;
          gt.jotPt[iJES][iJ] = gt.jotPt[0][iJ]*(1+relUnc);
          if(gt.jotPt[iJES][iJ] < 20) continue;
          gt.nJot[iJES]++;
          if(fabs(gt.jotEta[iJ])<2.4)
            gt.nJet[iJES]++;
        }
      }
    }
    float deltaPhiWH    = -1;
    float deltaPhiLep1Met = -1;
    float ptBalanceWH   = -1;
    float ptBalanceWHFJ = -1;
    float dEtaBjets     = -1;
    float dPhiBjets     = -1;
    float dRBjets       = -1;
    float dEtaLep1FJ    = -1;
    float dEtaLep1H     = -1;
    float fjECFN_2_4_20 = -1;
    float fjECFN_3_3_10 = -1;
    float fjECFN_2_3_20 = -1;
    float fjECFN_3_3_05 = -1;
    float psi022004031003 = -1;
    float psi022004022003 = -1;
    float psi022004030503 = -1;
    float dPhil1W  = -1;
    float dPhil1b1 = -1;
    float dPhil1b2 = -1;
    float dPhiWb1  = -1;
    float dPhiWb2  = -1;
    float dPhib1b2 = -1;
    float dEtal1W  = -1;
    float dEtal1b1 = -1;
    float dEtal1b2 = -1;
    float dEtaWb1  = -1;
    float dEtaWb2  = -1;
    float dEtab1b2 = -1;
    deltaPhiLep1Met = fabs(TVector2::Phi_mpi_pi(gt.pfmet[0] - lepton1Phi));
    if(isBoostedCategory) {
      ptBalanceWHFJ = gt.fjPt[0] / gt.topWBosonPt;
      dEtaLep1FJ = fabs(lepton1Eta - gt.fjEta);
      bLoad(b["fjHTTFRec"],ientry);
      bLoad(b["fjHTTMass"],ientry);
      // Need to handle these ECFs better for a more general case, long term to-do list
      bLoad(b["fjECFN_2_4_20"],ientry);
      bLoad(b["fjECFN_3_3_10"],ientry);
      bLoad(b["fjECFN_2_3_20"],ientry);
      bLoad(b["fjECFN_3_3_05"],ientry);
      fjECFN_2_4_20 = *((float*)b["fjECFN_2_4_20"]->GetAddress());
      fjECFN_3_3_10 = *((float*)b["fjECFN_3_3_10"]->GetAddress());
      fjECFN_2_3_20 = *((float*)b["fjECFN_2_3_20"]->GetAddress());
      fjECFN_3_3_05 = *((float*)b["fjECFN_3_3_05"]->GetAddress());
      psi022004031003 = fjECFN_2_4_20/pow(TMath::Max(0.0f,fjECFN_3_3_10),1.33);
      psi022004022003 = fjECFN_2_4_20/pow(TMath::Max(0.0f,fjECFN_2_3_20),1.00);
      psi022004030503 = fjECFN_2_4_20/pow(TMath::Max(0.0f,fjECFN_3_3_05),2.67);
      dPhil1W  = fabs(TVector2::Phi_mpi_pi(lepton1Phi      - gt.topWBosonPhi));
    } else {
      deltaPhiWH    = fabs(TVector2::Phi_mpi_pi(gt.hbbphi[0] - gt.topWBosonPhi));
      ptBalanceWH   = gt.hbbpt_reg[0] /  gt.topWBosonPt;
      dEtaBjets     = fabs(gt.jotEta[gt.hbbjtidx[0][0]]-gt.jotEta[gt.hbbjtidx[0][1]]);
      dPhiBjets     = fabs(TVector2::Phi_mpi_pi(gt.jotPhi[gt.hbbjtidx[0][0]]-gt.jotPhi[gt.hbbjtidx[0][1]]));
      dRBjets       = sqrt(dEtaBjets*dEtaBjets + dPhiBjets*dPhiBjets);
      dEtaLep1H = fabs(lepton1Eta - gt.hbbeta[0]);
      bLoad(b["sumEtSoft1"],ientry);
      bLoad(b["nSoft2"],ientry);
      bLoad(b["nSoft5"],ientry);
      bLoad(b["nSoft10"],ientry);
      bLoad(b["topMassLep1Met"],ientry);
      bLoad(b["topWBosonEta"],ientry);
      bLoad(b["topWBosonCosThetaCS"],ientry);
      bLoad(b["hbbCosThetaJJ"],ientry);
      bLoad(b["hbbCosThetaCSJ1"],ientry);
      dPhil1W  = fabs(TVector2::Phi_mpi_pi(lepton1Phi      - gt.topWBosonPhi));
      dPhil1b1 = fabs(TVector2::Phi_mpi_pi(lepton1Phi      - bjet1Phi       ));
      dPhil1b2 = fabs(TVector2::Phi_mpi_pi(lepton1Phi      - bjet2Phi       ));
      dPhiWb1  = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi - bjet1Phi       ));
      dPhiWb2  = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi - bjet2Phi       ));
      dPhib1b2 = fabs(TVector2::Phi_mpi_pi(bjet1Phi        - bjet2Phi       ));
      dEtal1W  = fabs(lepton1Eta      - gt.topWBosonEta);
      dEtal1b1 = fabs(lepton1Eta      - bjet1Eta       );
      dEtal1b2 = fabs(lepton1Eta      - bjet2Eta       );
      dEtaWb1  = fabs(gt.topWBosonEta - bjet1Eta       );
      dEtaWb2  = fabs(gt.topWBosonEta - bjet2Eta       );
      dEtab1b2 = fabs(bjet1Eta        - bjet2Eta       );
    }
    // deltaPhiWHFJ computed already in Jet multiplicity section
    // Set Selection Bits
    std::map<TString, bool> cut;
    for(unsigned iJES=0; iJES<NJES; iJES++) {
      if(iJES==0) {
        if(isBoostedCategory) {
          cut["WpTFJ"      ] = gt.topWBosonPt > 250;
          cut["dPhiWHFJ"   ] = fabs(gt.fjPhi - gt.topWBosonPhi) > 2.5;
          cut["btagFJ"     ] = gt.fjDoubleCSV > doubleBCut;
          cut["bvetoFJ"    ] = gt.fjDoubleCSV < doubleBCut;
        } else {
          cut["WpT"        ] = gt.topWBosonPt > 100;
          cut["dPhiWH"     ] = fabs(gt.fjPhi - gt.topWBosonPhi) > 2.5;
          cut["looseBTag"  ] = bjet1IsLoose;
          cut["tightBTag"  ] = bjet1IsTight;
          cut["mediumBVeto"] = !bjet1IsMedium;
          cut["looseBTag2" ] = bjet2IsLoose;
        }
        cut["lowMET"     ] = gt.pfmet[0] < 170;
        cut["boostedCat" ] = isBoostedCategory;
        cut["boostedVeto"] = !isBoostedCategory;
      }

      if(isBoostedCategory) {
        cut["mSD"     ] = gt.fjMSD_corr[iJES] >= 40;
        cut["mSD_SR"  ] = gt.fjMSD_corr[iJES] >= 80 && gt.fjMSD_corr[iJES]<150;
        cut["mSD_SB"  ] = cut["mSD"] && gt.fjMSD_corr[iJES]<80;
        cut["mSDVZ_SR"] = gt.fjMSD_corr[iJES] >= 40 && gt.fjMSD_corr[iJES]<120;
        cut["pTFJ"    ] = gt.fjPt[iJES] > 250;
        cut["0ijb"    ] = isojetNBtags[iJES]==0;
        cut["1ijb"    ] = !cut["0ijb"];
      } else {
        cut["dPhiWH"     ] = fabs(gt.hbbphi[iJES] - gt.topWBosonPhi) > 2.5;
        cut["pTjj"       ] = gt.hbbpt_reg[iJES] > 100; 
        cut["dPhiLep1Met"] = deltaPhiLep1Met < 2; 
        // Hardcoded Higgs(bb) mass window
        cut["mjj"        ] = gt.hbbm_reg[iJES] >= 90 && gt.hbbm_reg[iJES] < 150; 
        // Hardcoded Z(bb) mass window
        cut["mjjVZ"      ] = gt.hbbm_reg[iJES] >= 60 && gt.hbbm_reg[iJES] < 120; 
        // Sideband mass windows, can be changed with the ao.vzbbMode switch
        cut["mjjSBLo"    ] = gt.hbbm_reg[iJES] < mjjLo; 
        cut["mjjSBHi"    ] = gt.hbbm_reg[iJES] >= mjjHi && gt.hbbm_reg[iJES] < 250;
        cut["2jets"      ] = gt.nJet[iJES]==2;
        cut["2-3jets"    ] = gt.nJet[iJES]<4;
        cut["4+jets"     ] = gt.nJet[iJES]>=4;
        cut["metSig"     ] = gt.pfmetsig>2;
      } 
      selectionBits[iJES]=0; nMinusOneBits=0;
      if(passAllCuts(cut, ao.cuts[ao.selection])) {
        selectionBits[iJES] |= ao.selection;
        if(ao.debug>=1) printf("  Passed cuts for region %s\n",vhbbPlot::selectionNames[ao.selection].Data());
      }
      if(iJES==0 && passNMinusOne(cut, ao.cuts[ao.selection]))
        nMinusOneBits |= ao.selection;
    }

    // Begin Event weighting
    if(type==kData) {
      weight=1;
    } else {
      bLoad(b["normalizedWeight"],ientry);
      bLoad(b["pu"],ientry);
      float puWeight = nPUScaleFactor(ao.puWeights, gt.pu);
      weight_pileupUp = nPUScaleFactor(ao.puWeightsUp, gt.pu)/puWeight;
      weight_pileupDown = nPUScaleFactor(ao.puWeightsDown, gt.pu)/puWeight;
      
      weight = normalizedWeight * ao.lumi * puWeight * stitchWeight; 
      
      if(type==kWjets || type==kZjets) {
        bLoad(b["trueGenBosonPt"],ientry);
        if(type==kZjets)
          gt.sf_ewkV=ao.ZjetsEWKCorr[0]+ao.ZjetsEWKCorr[1]*
            (TMath::Power((gt.trueGenBosonPt+ao.ZjetsEWKCorr[2]),ao.ZjetsEWKCorr[3]));
        else
          gt.sf_ewkV=ao.WjetsEWKCorr[0]+ao.WjetsEWKCorr[1]*
            (TMath::Power((gt.trueGenBosonPt+ao.WjetsEWKCorr[2]),ao.WjetsEWKCorr[3]));
        if(isNLOWjets || isNLOZjets || isLowMassZjets) {
          gt.sf_qcdV=1;
        } else if(useHtBinnedVJetsKFactor) {
          bLoad(b["lheHT"],ientry);
          gt.sf_qcdV=qcdKFactor(type, gt.lheHT);
        } else {
          bLoad(b["lheHT"],ientry);
          bLoad(b["trueGenBosonPt"],ientry);
          if(type==kZjets)
            gt.sf_qcdV=ao.kfactors_ZJets->GetBinContent(ao.kfactors_ZJets->FindBin(
              TMath::Min(1.39,gt.trueGenBosonPt/1000.),
              TMath::Min(1.39,gt.lheHT         /1000.)
            ));
          else if(type==kWjets)
            gt.sf_qcdV=ao.kfactors_WJets->GetBinContent(ao.kfactors_WJets->FindBin(
              TMath::Min(1.39,gt.trueGenBosonPt/1000.),
              TMath::Min(1.39,gt.lheHT         /1000.)
            ));
        }
        if(ao.year==2017) { gt.sf_ewkV=1; } // hack
        weight *= gt.sf_qcdV * gt.sf_ewkV;
      }
      bLoad(b["scale"],ientry);
      double avgQCDScale = (gt.scale[0]+gt.scale[1]+gt.scale[2]+
                            gt.scale[3]+gt.scale[4]+gt.scale[5])/6.;
      if(avgQCDScale > 0) {
        weight_QCDr1f2 = gt.scale[0]/avgQCDScale;
        weight_QCDr1f5 = gt.scale[1]/avgQCDScale;
        weight_QCDr2f1 = gt.scale[2]/avgQCDScale;
        weight_QCDr2f2 = gt.scale[3]/avgQCDScale;
        weight_QCDr5f1 = gt.scale[4]/avgQCDScale;
        weight_QCDr5f5 = gt.scale[5]/avgQCDScale;
      } else {
        weight_QCDr1f2 = 1.0;
        weight_QCDr1f5 = 1.0;
        weight_QCDr2f1 = 1.0;
        weight_QCDr2f2 = 1.0;
        weight_QCDr5f1 = 1.0;
        weight_QCDr5f5 = 1.0;
      }

      if (typeLepSel==0) {
        bLoad(b["muonSfReco"],ientry);
        bLoad(b["muonSfTight"],ientry);
        bLoad(b["muonSfUnc"],ientry);
        bLoad(b["sf_muTrig"],ientry);
        weight *= gt.sf_muTrig; 
        weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
        weight_muSF = (1+gt.muonSfUnc[0]);
      } else if(typeLepSel==1) {
        bLoad(b["electronSfReco"],ientry);
        bLoad(b["electronSfMvaWP80"],ientry);
        bLoad(b["electronSfUnc"],ientry);
        bLoad(b["sf_eleTrig"],ientry);
        weight *= gt.sf_eleTrig; 
        weight *= gt.electronSfReco[0] * gt.electronSfMvaWP80[0];
        weight_elSF = (1+gt.electronSfUnc[0]);
      }

      float recorrect_vhEWK=1, recorrect_vhEWKUp=1, recorrect_vhEWKDown=1;
      if(type==kVZ) {
        bLoad(b["sf_wz"],ientry);
        bLoad(b["sf_zz"],ientry);
        weight *= gt.sf_wz * gt.sf_zz;
      } else if(type==kZH||type==kWH) {
        bLoad(b["sf_vh"],ientry);
        bLoad(b["sf_vhUp"],ientry);
        bLoad(b["sf_vhDown"],ientry);
        TString vhChannel= type==kWH? (lepton1Charge>0?"WplusH":"WminusH") : "ZllH";
        recorrect_vhEWK     = vhEWKCorr(vhChannel,gt.sf_vh    );
        recorrect_vhEWKUp   = vhEWKCorr(vhChannel,gt.sf_vhUp  );
        recorrect_vhEWKDown = vhEWKCorr(vhChannel,gt.sf_vhDown);
        weight *= recorrect_vhEWK;
        weight_VHCorrUp = recorrect_vhEWKUp/recorrect_vhEWK;
        weight_VHCorrDown = recorrect_vhEWKDown/recorrect_vhEWK;
      }
      
      // Hack for the central Btag weights 
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        unsigned iShift=0;
        GeneralTree::csvShift theShift = gt.csvShifts[iShift];
        for(unsigned iJ=0; iJ<jetPts[iPt][iEta].size(); iJ++) {
          unsigned absid = abs(jetFlavors[iPt][iEta][iJ]);
          BTagEntry::JetFlavor flav = absid == 5 ? BTagEntry::FLAV_B : 
            (absid == 4 ? BTagEntry::FLAV_C : BTagEntry::FLAV_UDSG);
          float reshapeFactor = ao.deepcsvSFs->eval_auto_bounds(
            GeneralTree::csvShiftName(theShift).Data(),
            flav,
            jetEtas[iPt][iEta][iJ], 
            jetPts[iPt][iEta][iJ],
            jetBtags[iPt][iEta][iJ]
          );
          if(reshapeFactor<0.001) reshapeFactor=1;
          weight *= reshapeFactor;
          weight_btag[0][iPt][iEta] = reshapeFactor;
        }
        
      }
      if((ao.selection==kWHSR || ao.selection==kWHVZbbCR || ao.selection==kWHFJSR || ao.selection==kWHVZbbFJCR) && ao.MVAVarType>1)
        weight *= sf_training;
    
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        // in 2017, we have to calculate the reshape factor for each jet in each kinematic bin
        for(unsigned iShift=1; iShift<GeneralTree::nCsvShifts; iShift++) {
          GeneralTree::csvShift theShift = gt.csvShifts[iShift];
          weight_btag[iShift][iPt][iEta] = 1.0;
          for(unsigned iJ=0; iJ<jetPts[iPt][iEta].size(); iJ++) {
            unsigned absid = abs(jetFlavors[iPt][iEta][iJ]);
            BTagEntry::JetFlavor flav = absid == 5 ? BTagEntry::FLAV_B : 
              (absid == 4 ? BTagEntry::FLAV_C : BTagEntry::FLAV_UDSG);
            float reshapeFactor = ao.deepcsvSFs->eval_auto_bounds(
              GeneralTree::csvShiftName(theShift).Data(),
              flav,
              jetEtas[iPt][iEta][iJ], 
              jetPts[iPt][iEta][iJ],
              jetBtags[iPt][iEta][iJ]
            );
            if(reshapeFactor>0.001) weight_btag[iShift][iPt][iEta] *= reshapeFactor/weight_btag[0][iPt][iEta];
            else                    weight_btag[iShift][iPt][iEta] *= 1;
          }
        }
      }
      if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) { // Boosted only weighting
      bLoad(b["fjHighestPtGen"],ientry);
        if(gt.fjHighestPtGen==21 && (
          category==kPlotWbb || category==kPlotWb || category==kPlotWLF ||
          category==kPlotZbb || category==kPlotZb || category==kPlotZLF)
        ) { weight_VGluUp = 1.2; weight_VGluDown = 0.8; }
        
        bLoad(b["fjPt"],ientry);
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Boosted_event_topologies
        if(category==kPlotZH || category==kPlotGGZH || category==kPlotWH || category==kPlotWbb) {
          bool passBB = ao.selection==kWHFJSR || ao.selection==kWHHeavyFlavorFJCR || ao.selection==kWHTT2bFJCR;
          float sf,sfUp,sfDown;
          float eff = 0.5; // empirical 
          if     (gt.fjPt[0]<350) { sf=0.92; sfUp=0.95; sfDown=0.89;}
          else if(gt.fjPt[0]<430) { sf=1.01; sfUp=1.04; sfDown=0.97;}
          else                    { sf=0.92; sfUp=0.95; sfDown=0.87;}
          float kBB,kBBUp,kBBDown; // effect of BB SF on the yield (depends on efficiency and pass/fail)
          if(passBB) {
            kBB     = sf;
            kBBUp   = sfUp;
            kBBDown = sfDown;
          } else {
            kBB     = (1.-eff*sf    )/(1.-eff);
            kBBUp   = (1.-eff*sfUp  )/(1.-eff);
            kBBDown = (1.-eff*sfDown)/(1.-eff);
          }
          weight *= kBB;
          weight_doubleBUp   = kBBUp  /kBB;
          weight_doubleBDown = kBBDown/kBB;
        }
      }
    }

    // Calculate the shape variable in all JES scenarios for all the regions
    float MVAVar[(int)shiftjes::N], bdtValue[(int)shiftjes::N];
    for(unsigned iJES=0; iJES<NJES; iJES++) {
      if(ao.MVAVarType>1) {
        // Only calculate it if we need to, it's expensive
        if((selectionBits[iJES] & ao.selection) != 0 && (
          (iJES==0 && ao.selection>=kWHLightFlavorCR && ao.selection<=kWHPresel) || 
          ao.selection==kWHSR || ao.selection==kWHVZbbCR
        )) {
          ao.mvaInputs[nThread][ 0] = gt.topWBosonPt                   ; // "WBosonPt"   
          ao.mvaInputs[nThread][ 1] = dPhil1W                          ; // "dPhil1W"    
          ao.mvaInputs[nThread][ 2] = gt.nSoft5                        ; // "nSoft5"     
          ao.mvaInputs[nThread][ 3] = gt.jotPt[iJES][gt.hbbjtidx[0][0]]; // "bjet1Pt"    
          ao.mvaInputs[nThread][ 4] = gt.jotPt[iJES][gt.hbbjtidx[0][1]]; // "bjet2Pt"    
          ao.mvaInputs[nThread][ 5] = bjet2btag                        ; // "bjet2btag"  
          ao.mvaInputs[nThread][ 6] = gt.hbbpt_reg[iJES]               ; // "hbbpt"      
          ao.mvaInputs[nThread][ 7] = gt.hbbm_reg[iJES]                ; // "hbbm"       
          ao.mvaInputs[nThread][ 8] = deltaPhiWH                       ; // "dPhiWH"     
          ao.mvaInputs[nThread][ 9] = ptBalanceWH                      ; // "ptBalanceWH"
          ao.mvaInputs[nThread][10] = gt.topMassLep1Met[0]             ; // "topMass"    
          ao.mvaInputs[nThread][11] = dRBjets                          ; // "dRBjets"    
          ao.mvaInputs[nThread][12] = dEtaLep1H                        ; // "dEtaLep1H"  
          ao.mvaInputs[nThread][13] = gt.nJet[0]-2                     ; // "nAddJet"    
          bdtValue[iJES] = ao.reader[nThread]->EvaluateMVA("BDT");
        } else if((selectionBits[iJES] & ao.selection) != 0 && (
          (iJES==0 && ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) || 
          ao.selection==kWHFJSR || ao.selection==kWHVZbbFJCR
        )) { 
          ao.mvaInputs[nThread][ 0] = dPhil1W                        ; //"dPhil1W"        
          ao.mvaInputs[nThread][ 1] = gt.topWBosonPt                 ; //"WBosonPt"       
          ao.mvaInputs[nThread][ 2] = lepton1Charge                  ; //"lepton1Charge"  
          ao.mvaInputs[nThread][ 3] = nIsojet[iJES]                  ; //"nIsojet"        
          ao.mvaInputs[nThread][ 4] = gt.fjMSD[iJES]                 ; //"MSD"            
          ao.mvaInputs[nThread][ 5] = gt.fjTau21SD                   ; //"Tau21SD"        
          ao.mvaInputs[nThread][ 6] = gt.fjTau32SD                   ; //"Tau32SD"        
          ao.mvaInputs[nThread][ 7] = gt.fjPt[iJES]                  ; //"fjPt"           
          ao.mvaInputs[nThread][ 8] = psi022004031003                ; //"psi022004031003"
          ao.mvaInputs[nThread][ 9] = psi022004022003                ; //"psi022004022003"
          ao.mvaInputs[nThread][10] = psi022004030503                ; //"psi022004030503"
          ao.mvaInputs[nThread][11] = gt.fjPt[iJES]/gt.topWBosonPt   ; //"ptBalanceWHFJ"  
          ao.mvaInputs[nThread][12] = dEtaLep1FJ                     ; //"dEtaLep1FJ"     
          ao.mvaInputs[nThread][13] = deltaPhiWHFJ                   ; //"dPhiWHFJ"       
          ao.mvaInputs[nThread][14] = gt.fjHTTFRec                   ; //"HTTFRec"        
          bdtValue[iJES] = ao.reader[nThread]->EvaluateMVA("BDT");
        }
      }
      switch(ao.MVAVarType) {
        case 1:
        default:
          if(ao.selection==kWHSR || ao.selection==kWHVZbbCR)
            MVAVar[iJES]=gt.hbbpt_reg[iJES];
          else if(ao.selection==kWHFJSR || ao.selection==kWHVZbbFJCR)
            MVAVar[iJES]=gt.fjPt[iJES];
          else if(ao.selection==kWHLightFlavorCR ||
            ao.selection==kWHHeavyFlavorLoMassCR || 
            ao.selection==kWHHeavyFlavorHiMassCR || 
            ao.selection==kWH2TopCR)
            MVAVar[iJES]=bjet2btag;
          else if(ao.selection==kWHLightFlavorFJCR ||
            ao.selection==kWHHeavyFlavorFJCR       ||
            ao.selection==kWHTT2bFJCR              ||
            ao.selection==kWHTT1bFJCR)
            MVAVar[iJES]=gt.fjMSD_corr[iJES];
          break;
        case 3:
          if(ao.selection==kWHSR || ao.selection==kWHVZbbCR || ao.selection==kWHFJSR || ao.selection==kWHVZbbFJCR)
            MVAVar[iJES]=bdtValue[iJES];
          else if(ao.selection==kWHLightFlavorCR ||
            ao.selection==kWHHeavyFlavorLoMassCR || 
            ao.selection==kWHHeavyFlavorHiMassCR || 
            ao.selection==kWH2TopCR)
            MVAVar[iJES]=bjet2btag;
          else if(ao.selection==kWHLightFlavorFJCR ||
            ao.selection==kWHHeavyFlavorFJCR       ||
            ao.selection==kWHTT2bFJCR              ||
            ao.selection==kWHTT1bFJCR)
            MVAVar[iJES]=gt.fjMSD_corr[iJES];
          break;
      }
    }
    
    bool trainingVeto = (
      ao.MVAVarType==1 || 
      (gt.eventNumber%10)>=3 || 
      category==kPlotData || 
      !(ao.selection==kWHVZbbCR || ao.selection==kWHSR || ao.selection==kWHVZbbFJCR || ao.selection==kWHFJSR)); 
    bool passFullSelNoTrainVeto = (selectionBits[0] & ao.selection) != 0;
    bool passFullSel = passFullSelNoTrainVeto && trainingVeto;

    // Fill the nominal histo and the shape uncertainty histos
    if(category==kData)  {
      if(passFullSel) {
        ao.histo_Baseline    [typeLepSel][category]->Fill(MVAVar[0], weight);
      }
    } else if(category>=kPlotVZbb)  {
      if(passFullSel) {
        ao.histo_Baseline    [typeLepSel][category]->Fill(MVAVar[0], weight);
        ao.histo_pileupUp    [typeLepSel][category]->Fill(MVAVar[0], weight*weight_pileupUp);
        ao.histo_pileupDown  [typeLepSel][category]->Fill(MVAVar[0], weight*weight_pileupDown);
        ao.histo_VHCorrUp    [typeLepSel][category]->Fill(MVAVar[0], weight*weight_VHCorrUp);
        ao.histo_VHCorrDown  [typeLepSel][category]->Fill(MVAVar[0], weight*weight_VHCorrDown);
        ao.histo_QCDr1f2     [typeLepSel][category]->Fill(MVAVar[0], weight*weight_QCDr1f2);
        ao.histo_QCDr1f5     [typeLepSel][category]->Fill(MVAVar[0], weight*weight_QCDr1f5);
        ao.histo_QCDr2f1     [typeLepSel][category]->Fill(MVAVar[0], weight*weight_QCDr2f1);
        ao.histo_QCDr2f2     [typeLepSel][category]->Fill(MVAVar[0], weight*weight_QCDr2f2);
        ao.histo_QCDr5f1     [typeLepSel][category]->Fill(MVAVar[0], weight*weight_QCDr5f1);
        ao.histo_QCDr5f5     [typeLepSel][category]->Fill(MVAVar[0], weight*weight_QCDr5f5);
        ao.histo_eleSFUp     [typeLepSel][category]->Fill(MVAVar[0], weight*weight_elSF);
        ao.histo_eleSFDown   [typeLepSel][category]->Fill(MVAVar[0], weight/weight_elSF);
        ao.histo_muSFUp      [typeLepSel][category]->Fill(MVAVar[0], weight*weight_muSF);
        ao.histo_muSFDown    [typeLepSel][category]->Fill(MVAVar[0], weight/weight_muSF);
        for(unsigned iPt=0; iPt<5; iPt++)
        for(unsigned iEta=0; iEta<3; iEta++)
        for(unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
          GeneralTree::csvShift shift = gt.csvShifts[iShift];
          if (shift==GeneralTree::csvCent) continue;
          ao.histo_btag[iShift][iPt][iEta][typeLepSel][category]->Fill(MVAVar[0], weight*weight_btag[iShift][iPt][iEta]);
        }
        if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) { // Boosted only weighting
          ao.histo_VGluUp      [typeLepSel][category]->Fill(MVAVar[0], weight*weight_VGluUp);
          ao.histo_VGluDown    [typeLepSel][category]->Fill(MVAVar[0], weight*weight_VGluDown);
          ao.histo_doubleBUp   [typeLepSel][category]->Fill(MVAVar[0], weight*weight_doubleBUp);
          ao.histo_doubleBDown [typeLepSel][category]->Fill(MVAVar[0], weight*weight_doubleBDown);
        }
      }
      for(unsigned iJES=0; iJES<NJES; iJES++) {
        bool passFullSelJES = (selectionBits[iJES] & ao.selection) != 0 && trainingVeto;
        if(!passFullSelJES) continue;
        ao.histo_jes[iJES][typeLepSel][category]->Fill(MVAVar[iJES], weight);
      }
    }

    // Fill the plotting histograms and MVA tree (if applicable)
    bLoad(b["mT"],ientry);
    bLoad(b["pfmet"],ientry);
    bLoad(b["pfmetsig"],ientry);
    bLoad(b["topWBosonPhi"],ientry);
    bLoad(b["fjMSD_corr"],ientry);
    bLoad(b["fjPt"],ientry);
    bLoad(b["fjDoubleCSV"],ientry);

    if(passFullSelNoTrainVeto) {
      if(ao.debug>=3) printf("\tPassed this sel\n");
      // Lock the mutex and fill the MVA tree
      if((ao.selection==kWHSR || ao.selection==kWHVZbbCR) && category!=kPlotData) {
        mvaTreeMutex.lock();
        ao.mva_weight        = weight                   ; 
        ao.mva_category      = category                 ;
        ao.mva_eventNumber   = gt.eventNumber           ;
        ao.mva_mT            = gt.mT[0]                 ;
        ao.mva_dPhil1W       = dPhil1W                  ;
        ao.mva_WBosonPt      = gt.topWBosonPt           ;
        ao.mva_dPhiLep1Met   = deltaPhiLep1Met          ;
        ao.mva_lepton1Pt     = lepton1Pt                ;
        ao.mva_lepton1Eta    = lepton1Eta               ;
        ao.mva_lepton1Charge = lepton1Charge            ;
        ao.mva_pfmet         = gt.pfmet[0]              ;
        ao.mva_sumEtSoft1    = gt.sumEtSoft1            ;
        ao.mva_nSoft2        = gt.nSoft2                ;
        ao.mva_nSoft5        = gt.nSoft5                ;
        ao.mva_nSoft10       = gt.nSoft10               ;
        ao.mva_bjet1Pt       = bjet1Pt                  ;
        ao.mva_bjet1btag     = bjet1btag                ;
        ao.mva_bjet2Pt       = bjet2Pt                  ;
        ao.mva_bjet2btag     = bjet2btag                ;
        ao.mva_hbbpt         = gt.hbbpt_reg[0]          ;
        ao.mva_hbbm          = gt.hbbm_reg[0]           ;
        ao.mva_dPhiWH        = deltaPhiWH               ;
        ao.mva_ptBalanceWH   = ptBalanceWH              ;
        ao.mva_topMass       = gt.topMassLep1Met[0]     ;
        ao.mva_dRBjets       = dRBjets                  ;
        ao.mva_dEtaLep1H     = dEtaLep1H                ;
        ao.mva_nAddJet       = gt.nJet[0]-2             ;
        ao.mvaTree->Fill();
        mvaTreeMutex.unlock();
      } else if((ao.selection==kWHFJSR || ao.selection==kWHVZbbFJCR) && category!=kPlotData) {
        mvaTreeMutex.lock();
        ao.mva_weight           = weight                   ; 
        ao.mva_category         = category                 ;
        ao.mva_eventNumber      = gt.eventNumber           ;
        ao.mva_mT               = gt.mT[0]                 ;
        ao.mva_dPhil1W          = dPhil1W                  ;
        ao.mva_WBosonPt         = gt.topWBosonPt           ;
        ao.mva_dPhiLep1Met      = deltaPhiLep1Met          ;
        ao.mva_lepton1Pt        = lepton1Pt                ;
        ao.mva_lepton1Eta       = lepton1Eta               ;
        ao.mva_lepton1Charge    = lepton1Charge            ;
        ao.mva_pfmet            = gt.pfmet[0]              ;
        ao.mva_MSD              = gt.fjMSD_corr[0]         ;
        ao.mva_nIsojet          = nIsojet[0]               ; 
        ao.mva_Tau21SD          = gt.fjTau21SD             ; 
        ao.mva_Tau32SD          = gt.fjTau32SD             ; 
        ao.mva_fjPt             = gt.fjPt[0]               ; 
        ao.mva_psi022004031003  = psi022004031003          ; 
        ao.mva_psi022004022003  = psi022004022003          ; 
        ao.mva_psi022004030503  = psi022004030503          ; 
        ao.mva_ptBalanceWHFJ    = ptBalanceWHFJ            ; 
        ao.mva_dEtaLep1FJ       = dEtaLep1FJ               ; 
        ao.mva_dPhiWHFJ         = deltaPhiWHFJ             ; 
        ao.mva_HTTFRec          = gt.fjHTTFRec             ; 
        ao.mvaTree->Fill();
        mvaTreeMutex.unlock();
      }
    }
    float theVar;
    for(int p=0; p<nPlots; p++) { 
      bool makePlot=false;
      // Variables -- change the makePlot for n-1 later
      // common
      if      (ao.histoNames[p]=="MVAVar"                  ) { theVar = MVAVar[0]                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="lepton1Pt"               ) { theVar = lepton1Pt                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="lepton1Eta"              ) { theVar = lepton1Eta                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="lepton1Charge"           ) { theVar = lepton1Charge              ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="WBosonPt"                ) { theVar = gt.topWBosonPt             ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="WBosonPhi"               ) { theVar = gt.topWBosonPhi            ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="mT"                      ) { theVar = gt.mT[0]                   ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="pfmet"                   ) { theVar = gt.pfmet[0]                ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="pfmetsig"                ) { theVar = gt.pfmetsig                ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bdtValue"                ) { theVar = bdtValue[0]                ; makePlot = passFullSel; }
      // boosted
      else if (ao.histoNames[p]=="mSD"                     ) { theVar = gt.fjMSD_corr[0]           ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="pTFJ"                    ) { theVar = gt.fjPt[0]                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="Tau21SD"                 ) { theVar = gt.fjTau21SD               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="Tau32SD"                 ) { theVar = gt.fjTau32SD               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="doubleB"                 ) { theVar = gt.fjDoubleCSV             ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="deltaEtaLep1FJ"          ) { theVar = dEtaLep1FJ                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="deltaPhiWHFJ"            ) { theVar = deltaPhiWHFJ               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ptBalanceWHFJ"           ) { theVar = ptBalanceWHFJ              ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nIsojet"                 ) { theVar = nIsojet[0]                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="isojetNBtags"            ) { theVar = isojetNBtags[0]            ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="fjHTTFRec"               ) { theVar = gt.fjHTTFRec               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="fjHTTMass"               ) { theVar = gt.fjHTTMass               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="psi022004031003"         ) { theVar = psi022004031003            ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="psi022004022003"         ) { theVar = psi022004022003            ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="psi022004030503"         ) { theVar = psi022004030503            ; makePlot = passFullSel; }
      // resolved
      else if (ao.histoNames[p]=="Mjj"                     ) { theVar = gt.hbbm_reg[0]             ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="pTjj"                    ) { theVar = gt.hbbpt_reg[0]            ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet1Pt"                 ) { theVar = bjet1Pt                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet2Pt"                 ) { theVar = bjet2Pt                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet1btag"               ) { theVar = bjet1btag                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet2btag"               ) { theVar = bjet2btag                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="deltaPhiWH"              ) { theVar = deltaPhiWH                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ptBalanceWH"             ) { theVar = ptBalanceWH                ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nJet"                    ) { theVar = gt.nJet[0]                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="sumEtSoft1"              ) { theVar = gt.sumEtSoft1              ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nSoft2"                  ) { theVar = gt.nSoft2                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nSoft5"                  ) { theVar = gt.nSoft5                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nSoft10"                 ) { theVar = gt.nSoft10                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="topWBosonCosThetaCS"     ) { theVar = gt.topWBosonCosThetaCS[0]  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="hbbCosThetaJJ"           ) { theVar = gt.hbbCosThetaJJ[0]        ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="hbbCosThetaCSJ1"         ) { theVar = gt.hbbCosThetaCSJ1[0]      ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="deltaEtaLep1H"           ) { theVar = dEtaLep1H                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dPhil1W"                 ) { theVar = dPhil1W                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dPhil1b1"                ) { theVar = dPhil1b1                   ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dPhil1b2"                ) { theVar = dPhil1b2                   ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dPhiWb1"                 ) { theVar = dPhiWb1                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dPhiWb2"                 ) { theVar = dPhiWb2                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dPhib1b2"                ) { theVar = dPhib1b2                   ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dEtal1W"                 ) { theVar = dEtal1W                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dEtal1b1"                ) { theVar = dEtal1b1                   ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dEtal1b2"                ) { theVar = dEtal1b2                   ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dEtaWb1"                 ) { theVar = dEtaWb1                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dEtaWb2"                 ) { theVar = dEtaWb2                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dEtab1b2"                ) { theVar = dEtab1b2                   ; makePlot = passFullSel; }
      if(!makePlot) continue;
      if(ao.histoNames[p]=="MVAVar")
        theVar = TMath::Min((float)(ao.MVAbins[ao.MVAbins.size()-1]-0.00001), (float)theVar);
      else
        theVar=TMath::Min((float)(ao.xmax[p]-0.00001), (float)theVar);
      if(ao.debug>=3) printf("\t\tFilling %s\n",ao.histoNames[p].Data());
      ao.histos[typeLepSel][p][category]->Fill(theVar, weight);
    }
       
  } // End Event Loop
}

void writeDatacards(analysisObjects &ao, TString dataCardDir) {  
  for(unsigned lep=0; lep<nLepSel; lep++) {
    ofstream newcardShape;
    newcardShape.open(Form("MitVHBBAnalysis/datacards/%s/datacard_W%sH%s.txt",
                      dataCardDir.Data(),leptonStrings[lep].Data(),selectionNames[ao.selection].Data()));
    newcardShape << Form("imax * number of channels\n");
    newcardShape << Form("jmax * number of background minus 1\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");

    newcardShape << Form("shapes      *   * datacard_W%sH%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",leptonStrings[lep].Data(),selectionNames[ao.selection].Data());
    newcardShape << Form("shapes data_obs * datacard_W%sH%s.root  histo_%s\n",leptonStrings[lep].Data(),selectionNames[ao.selection].Data(),plotBaseNames[kPlotData].Data());

    newcardShape << Form("Observation %f\n",ao.histo_Baseline[lep][kPlotData]->GetSumOfWeights());

    newcardShape << Form("bin   ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
      newcardShape << Form("chW%sH%s  ",leptonStrings[lep].Data(),selectionNames[ao.selection].Data());
    newcardShape << Form("\n");

    newcardShape << Form("process   ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
      newcardShape << Form("%s  ",plotBaseNames[ic].Data());
    newcardShape << Form("\n");
 
    newcardShape << Form("process  ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() <= 0)
        continue;
      if(ic==kPlotWH)
        newcardShape << Form("%d  ",0);
      else if(ic==kPlotZH)
        newcardShape << Form("%d  ",-1);
      else if(ic==kPlotGGZH)
        newcardShape << Form("%d  ",-2);
      else
        newcardShape << Form("%d  ",ic);
    }
    newcardShape << Form("\n");

    newcardShape << Form("rate  ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
      newcardShape << Form("%f  ",ao.histo_Baseline[lep][ic]->GetSumOfWeights());
    newcardShape << Form("\n");

    newcardShape << Form("lumi_13TeV    lnN     ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
       newcardShape << Form("%6.3f ",1.025);
    newcardShape << Form("\n");

    newcardShape << Form("pileup    shape   ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
      newcardShape << Form("1.0  ");
    newcardShape << Form("\n");

    newcardShape << Form("VHCorr    shape   ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++){
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() <= 0)
        continue;
      if(ic!=kPlotZH && ic!=kPlotGGZH) 
        newcardShape << Form("-  ");
      else
        newcardShape << Form("1.0  ");
    }
    newcardShape << Form("\n");

    newcardShape << Form("eleSF    shape   ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
      newcardShape << Form("1.0  ");
    newcardShape << Form("\n");

    newcardShape << Form("muSF    shape   ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
      newcardShape << Form("1.0  ");
    newcardShape << Form("\n");

    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0) {
      newcardShape << Form("QCDScale%s    shape   ",plotBaseNames[ic].Data());
      for(unsigned ic2=kPlotVZbb; ic2!=nPlotCategories; ic2++) {
        if(ao.histo_Baseline[lep][ic2]->GetSumOfWeights() > 0) {
          if(ic==ic2) newcardShape << Form("1.0  ");
          else        newcardShape << Form("-  ");
        }
      }
      newcardShape << Form("\n");
    }

    for(unsigned iJES=0; iJES<NJES; iJES++) {
      if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
      TString shiftName(jesName(static_cast<shiftjes>(iJES)).Data());
      if(!shiftName.EndsWith("Up")) continue;
      TString nuisanceName = shiftName(0,shiftName.Length()-2);
      newcardShape << Form("%s    shape   ",nuisanceName.Data());
      for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
        newcardShape << Form("1.0  ");
      newcardShape << Form("\n");
    }
    
    GeneralTree gt;
    for(unsigned iPt=0; iPt<5; iPt++)
    for(unsigned iEta=0; iEta<3; iEta++)
    for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
      GeneralTree::csvShift shift = gt.csvShifts[iShift];
      if (shift==GeneralTree::csvCent) continue;
      // Get the name of the nuisance only from the Up variations
      TString shiftName(btagShiftName(shift));
      if(!shiftName.EndsWith("Up")) continue;
      TString nuisanceName = shiftName(0,shiftName.Length()-2);
      newcardShape << Form("CMS_VH_btag%d_pt%d_eta%d_%s    shape   ",ao.year,   iPt,iEta,nuisanceName.Data());
      for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0) {
        if((ic==kPlotWbb||ic==kPlotWb||ic==kPlotWLF) && ao.selection==kWH2TopCR) 
          newcardShape << Form("-  ");
        else
          newcardShape << Form("1.0  ");
      }
      newcardShape << Form("\n");
    }

    newcardShape<<"pdf_qqbar lnN "; 
    for(int ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights()<=0)
        continue;
      if(ic==kPlotVZbb||ic==kPlotVVLF||ic==kPlotWbb||ic==kPlotWb||ic==kPlotWLF||
         ic==kPlotZbb||ic==kPlotZb||ic==kPlotZLF||ic==kPlotZH)
        newcardShape<<pdfAcceptUncs[ic]<<" ";
      else newcardShape<<"- ";
    } 
    newcardShape<<std::endl;

    newcardShape<<"pdf_gg lnN ";
    for(int ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights()<=0)
        continue;
      if(ic==kPlotTop||ic==kPlotTT|ic==kPlotGGZH)
        newcardShape<<pdfAcceptUncs[ic]<<" "; 
      else
        newcardShape<<"- ";
    }
    newcardShape<<std::endl;

    //newcardShape<<Form("CMS_VH_TopNorm lnN ");
    //for(int ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
    //  if(ao.histo_Baseline[lep][ic]->GetSumOfWeights()<=0)
    //    continue;
    //  newcardShape<< (ic==kPlotTop? "1.15 ":"- ");
    //}
    //newcardShape<<std::endl;
    //
    //newcardShape<<Form("CMS_VH_VVNorm lnN ");
    //for(int ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
    //  if(ao.histo_Baseline[lep][ic]->GetSumOfWeights()<=0)
    //    continue;
    //  newcardShape<< ((ic==kPlotVZbb||ic==kPlotVVLF)? "1.15 ":"- ");
    //}
    //newcardShape<<std::endl;
    
    // Normalization and double B scale factors
    if(ao.selection>=kWHLightFlavorFJCR && ao.selection<=kWHFJPresel) {
      // Boosted norms
      newcardShape << Form("SF%d_WHFFJ_Wln rateParam * %s 1 [0.2,5]\n",ao.year,plotBaseNames[kPlotWbb].Data());
      newcardShape << Form("SF%d_WHFFJ_Wln rateParam * %s 1 [0.2,5]\n",ao.year,plotBaseNames[kPlotWb].Data());
      newcardShape << Form("SF%d_WLFFJ_Wln rateParam * %s 1 [0.2,5]\n",ao.year,plotBaseNames[kPlotWLF].Data());
      newcardShape << Form("SF%d_TTFJ_Wln rateParam * %s 1 [0.2,5]\n",ao.year,plotBaseNames[kPlotTop].Data());
      newcardShape << Form("SF%d_TTFJ_Wln rateParam * %s 1 [0.2,5]\n",ao.year,plotBaseNames[kPlotTT].Data());
      // In situ measurement of doubleB SF for WLF, Wb, eff checked in kWHFJPresel
      newcardShape << Form("eff%dDoubleB_WLF  param 0.024 0.0060\n",ao.year);
      newcardShape << Form("eff%dDoubleB_Wb   param 0.20  0.05\n",ao.year);
      newcardShape << Form("eff%dSFDoubleB_WLF param 1.0 0.5\n",ao.year);
      newcardShape << Form("eff%dSFDoubleB_Wb  param 1.0 0.5\n",ao.year);
      //newcardShape << Form("eff%dSFDoubleB_WLF extArg 1.0 [0.1,10]\n",ao.year);
      //newcardShape << Form("eff%dSFDoubleB_Wb  extArg 1.0 [0.1,10]\n",ao.year);
      if(ao.selection==kWHHeavyFlavorFJCR || ao.selection==kWHTT2bFJCR || ao.selection==kWHFJSR) {
        newcardShape << Form("passBB%d_WLF rateParam * %s (@0*1.0) eff%dSFDoubleB_WLF\n",ao.year,plotBaseNames[kPlotWLF].Data(),ao.year);
        newcardShape << Form("passBB%d_Wb  rateParam * %s (@0*1.0) eff%dSFDoubleB_Wb \n",ao.year,plotBaseNames[kPlotWb].Data(),ao.year);
      } else if(ao.selection==kWHLightFlavorFJCR || ao.selection==kWHTT1bFJCR) {
        newcardShape << Form("failBB%d_WLF rateParam * %s ((1.0-@0*@1)/(1.0-@1)) eff%dSFDoubleB_WLF,eff%dDoubleB_WLF\n",ao.year,plotBaseNames[kPlotWLF].Data(),ao.year,ao.year);
        newcardShape << Form("failBB%d_Wb  rateParam * %s ((1.0-@0*@1)/(1.0-@1)) eff%dSFDoubleB_Wb,eff%dDoubleB_Wb\n"  ,ao.year,plotBaseNames[kPlotWb].Data(),ao.year,ao.year);
      }
      // Shape uncertainty using BTV doubleB SF for Zbb, VZ(bb), ZH
      newcardShape << Form("CMS_doubleB shape ");
      for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
      if(ao.histo_Baseline[lep][ic] && ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0) {
        if(ic==kPlotWbb||ic==kPlotVZbb||ic==kPlotZH||ic==kPlotGGZH)
          newcardShape << "1.0 ";
        else
          newcardShape << "- ";
      }
      newcardShape<<std::endl;
      newcardShape<<Form("VjetsGluFrac shape ");
      for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
      if(ao.histo_Baseline[lep][ic] && ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0) {
        if(ic==kPlotWbb||ic==kPlotWb||ic==kPlotWLF)
          newcardShape << "1.0 ";
        else
          newcardShape << "- ";
      }
      newcardShape<<std::endl;
    } else {
      newcardShape << Form("SF%d_TT_Wln  rateParam * %s 1 [0.2,5]\n" ,ao.year,plotBaseNames[kPlotTT].Data());
      newcardShape << Form("SF%d_TT_Wln  rateParam * %s 1 [0.2,5]\n" ,ao.year,plotBaseNames[kPlotTop].Data());
      newcardShape << Form("SF%d_Wbb_Wln rateParam * %s 1 [0.2,5]\n" ,ao.year,plotBaseNames[kPlotWbb].Data());
      newcardShape << Form("SF%d_Wb_Wln  rateParam * %s 1 [0.2,5]\n" ,ao.year,plotBaseNames[kPlotWb].Data());
      newcardShape << Form("SF%d_WLF_Wln  rateParam * %s 1 [0.2,5]\n",ao.year,plotBaseNames[kPlotWLF].Data());
    }

    newcardShape << Form("* autoMCStats 1\n");
    newcardShape.close();
  }
}

void datacardsFromHistograms(
  TString dataCardDir,
  vhbbPlot::selectionType selection,
  bool useBoostedCategory=false,
  int MVAVarType=3,
  unsigned year=2016
) {
  struct analysisObjects ao;
  ao.useBoostedCategory=useBoostedCategory;
  ao.selection = selection;
  ao.year = year;
  GeneralTree gt;
  // Get shape histograms from file
  for(unsigned lep=0; lep<nLepSel; lep++) {

    char inFileDatacardsName[200];
    sprintf(
      inFileDatacardsName,
      "MitVHBBAnalysis/datacards/%s/datacard_W%sH%s.root",
      dataCardDir.Data(),
      leptonStrings[lep].Data(),
      selectionNames[selection].Data()
    );
    TFile* infile = new TFile(inFileDatacardsName,"read");

    for(unsigned ic=kPlotData; ic!=nPlotCategories; ic++) {
      ao.histo_Baseline    [lep][ic] = (TH1F*)infile->Get(Form("histo_%s"                , plotBaseNames[ic].Data()));
      if(!ao.histo_Baseline[lep][ic]) continue;
      ao.histo_Baseline    [lep][ic]->SetDirectory(0);
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() <= 0 && ic!=kPlotData) continue;
      if(ic<kPlotVZbb) continue;
      ao.histo_pileupUp    [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_pileupUp"       , plotBaseNames[ic].Data()));
      ao.histo_pileupDown  [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_pileupDown"     , plotBaseNames[ic].Data()));
      ao.histo_VHCorrUp    [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_VHCorrUp"       , plotBaseNames[ic].Data()));
      ao.histo_VHCorrDown  [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_VHCorrDown"     , plotBaseNames[ic].Data()));
      ao.histo_QCDScaleUp  [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_QCDScale%sUp"   , plotBaseNames[ic].Data(),plotBaseNames[ic].Data()));
      ao.histo_QCDScaleDown[lep][ic] = (TH1F*)infile->Get(Form("histo_%s_QCDScale%sDown" , plotBaseNames[ic].Data(),plotBaseNames[ic].Data()));
      ao.histo_eleSFUp     [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_eleSFUp"        , plotBaseNames[ic].Data()));                         
      ao.histo_eleSFDown   [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_eleSFDown"      , plotBaseNames[ic].Data()));                         
      ao.histo_muSFUp      [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_muSFUp"         , plotBaseNames[ic].Data()));                         
      ao.histo_muSFDown    [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_muSFDown"       , plotBaseNames[ic].Data()));                         
      ao.histo_pileupUp    [lep][ic]->SetDirectory(0);
      ao.histo_pileupDown  [lep][ic]->SetDirectory(0);
      ao.histo_VHCorrUp    [lep][ic]->SetDirectory(0);
      ao.histo_VHCorrDown  [lep][ic]->SetDirectory(0);
      ao.histo_QCDScaleUp  [lep][ic]->SetDirectory(0);
      ao.histo_QCDScaleDown[lep][ic]->SetDirectory(0);
      ao.histo_eleSFUp     [lep][ic]->SetDirectory(0);
      ao.histo_eleSFDown   [lep][ic]->SetDirectory(0);
      ao.histo_muSFUp      [lep][ic]->SetDirectory(0);
      ao.histo_muSFDown    [lep][ic]->SetDirectory(0);
      for(unsigned iJES=0; iJES<NJES; iJES++) { 
        if(iJES==(unsigned)shiftjes::kJESTotalUp || iJES==(unsigned)shiftjes::kJESTotalDown) continue;
        ao.histo_jes[iJES][lep][ic] = (TH1F*)infile->Get(Form("histo_%s_%s", plotBaseNames[ic].Data(), jesName(static_cast<shiftjes>(iJES)).Data()));
        ao.histo_jes[iJES][lep][ic]->SetDirectory(0);
      }
      for(unsigned iPt=0; iPt<5; iPt++)
      for(unsigned iEta=0; iEta<3; iEta++)
      for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
        GeneralTree::csvShift shift = gt.csvShifts[iShift];
        if (shift==GeneralTree::csvCent) continue;
        ao.histo_btag[iShift][iPt][iEta][lep][ic] = (TH1F*)infile->Get(Form("histo_%s_CMS_VH_btag%d_pt%d_eta%d_%s",plotBaseNames[ic].Data(),year,iPt,iEta,btagShiftName(shift)));
        ao.histo_btag[iShift][iPt][iEta][lep][ic]->SetDirectory(0);
      }
      if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) {
        ao.histo_VGluUp     [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_VjetsGluFracUp"    , plotBaseNames[ic].Data()));
        ao.histo_VGluDown   [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_VjetsGluFracDown"  , plotBaseNames[ic].Data()));
        ao.histo_doubleBUp  [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_CMS_doubleBUp"     , plotBaseNames[ic].Data()));
        ao.histo_doubleBDown[lep][ic] = (TH1F*)infile->Get(Form("histo_%s_CMS_doubleBDown"   , plotBaseNames[ic].Data()));
        ao.histo_VGluUp     [lep][ic]->SetDirectory(0);
        ao.histo_VGluDown   [lep][ic]->SetDirectory(0);
        ao.histo_doubleBUp  [lep][ic]->SetDirectory(0);
        ao.histo_doubleBDown[lep][ic]->SetDirectory(0);
      }
    }
    infile->Close();
  }
  
  writeDatacards(ao, dataCardDir);
}
