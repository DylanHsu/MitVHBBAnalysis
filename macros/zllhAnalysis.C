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
#include "PandaAnalysis/Utilities/src/CSVHelper.cc"

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
#include "MitAnalysisRunII/panda/macros/80x/auxiliar.h"
#include "MitAnalysisRunII/panda/macros/80x/common.h"
#include "vhbbPlot.h"
#include "PandaAnalysis/Flat/interface/Common.h"

TString ntupleDir2016 = "/mnt/hadoop/scratch/dhsu/dylansVHSkims/2016/v_009_vhbb1/dilepton";
TString ntupleDir2017 = "/mnt/hadoop/scratch/dhsu/dylansVHSkims/2017/v_010_vhbb2/dilepton";
const bool useHtBinnedVJetsKFactor=false;
const int NJES = (int)shiftjes::N; // Number of JES variations
const int nLepSel=3; // Number of lepton selections
const int nPlots=50; // Max number of plots
const unsigned char nBinsZpt = 2;
const int nThreads=10;
vector<float> binsZpt = {50,150,3000};
vector<TString> leptonStrings={"mm","ee","em"};
float sf_training=1.4286;
using namespace vhbbPlot;
std::mutex mvaTreeMutex;

struct analysisObjects {
  // Physics
  selectionType selection;
  std::map<selectionType,vector<TString>> cuts;
  float isojetBtagCut;
  // Configuration parameters
  bool useBoostedCategory;
  int MVAVarType; 
  unsigned year;
  unsigned debug;
  double lumi;
  unsigned long whichTriggers;
  char binZpt=-1;
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
  TH2D *kfactors_ZJets=0;
  BTagCalibrationReader *deepcsvSFs=0; 
  BTagCalibration *deepcsvCalib=0; 
  CSVHelper *cmvaReweighter=0;
  vector<double> ZjetsEWKCorr;
  // MVA training tree
  TTree *mvaTree=0; TFile *mvaFile=0;
  // Common vars
  float mva_ZBosonPt, mva_ZBosonM, mva_CosThetaCS, mva_CosThetaStar;
  float mva_weight;
  unsigned char mva_category;
  ULong64_t mva_eventNumber;
  // Resolved only vars
  float mva_sumEtSoft1, mva_nSoft2, mva_nSoft5, mva_nSoft10;
  float mva_bjet1Pt, mva_bjet2Pt, mva_bjet1btag, mva_bjet2btag;
  float mva_hbbpt, mva_hbbm, mva_dPhiZH, mva_ptBalanceZH;
  float mva_dRBjets, mva_dEtaBjets, mva_dRZH, mva_dEtaZH;
  float mva_nAddJet;
  // Boosted only vars
  float mva_nIsojet, mva_MSD, mva_Tau21SD, mva_fjPt;
  float mva_ptBalanceZHFJ, mva_dEtaZHFJ, mva_dPhiZHFJ, mva_mTZHFJ;
  float mva_ptBalanceL1L2, mva_dRL1L2;
  float mva_lepton1Pt, mva_lepton2Pt, mva_lepton1Eta, mva_lepton2Eta, mva_deltaM;
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
// binZpt:
//   must be less than nBinsZpt, negative turns off Zpt binning
// batchSampleName, type:

void zllhAnalysis(
  TString dataCardDir,
  vhbbPlot::selectionType selection,
  bool useBoostedCategory=false,
  int MVAVarType=3,
  char binZpt=-1, 
  unsigned year=2016,
  unsigned debug=0,
  bool multithread=false,
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
  ao.binZpt = binZpt;
  TString ntupleDir = (year==2016)? ntupleDir2016:ntupleDir2017;

  // Analysis Cuts
  ao.isojetBtagCut = (ao.year==2017)? deepcsvLoose : cmvaLoose;
  ao.cuts[kZllHLightFlavorCR  ] = {"ZpT","bveto","Zmass"                               , "boostedVeto"             };
  ao.cuts[kZllHHeavyFlavorCR  ] = {"ZpT","btag" ,"ZmassTight","lowMET","dPhiZH","mjjSB", "boostedVeto"             };
  ao.cuts[kZllH2TopCR         ] = {"ZpT","btag" ,"ZmassSB"                             , "boostedVeto"             };
  ao.cuts[kZllHSR             ] = {"ZpT","btag" ,"Zmass"              ,"dPhiZH","mjj"  , "boostedVeto", "vetoTrain"};
  ao.cuts[kZllHPresel         ] = {"ZpT"                                               , "boostedVeto"             };
  ao.cuts[kZllHLightFlavorFJCR] = {"boostedCat","ZpTFJ","pTFJ","dPhiZHFJ","mSD"   ,         "Zmass"     , "bvetoFJ"             };
  ao.cuts[kZllHHeavyFlavorFJCR] = {"boostedCat","ZpTFJ","pTFJ","dPhiZHFJ","mSD_SB",         "Zmass"     , "btagFJ"              };
  ao.cuts[kZllHTT1bFJCR       ] = {"boostedCat","ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "1ijb", "Zmass"     , "bvetoFJ"             };
  ao.cuts[kZllHTT2bFJCR       ] = {"boostedCat","ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "1ijb", "Zmass"     , "btagFJ"              };
  ao.cuts[kZllHFJSR           ] = {"boostedCat","ZpTFJ","pTFJ","dPhiZHFJ","mSD_SR",         "Zmass"     , "btagFJ" , "vetoTrain"};
  ao.cuts[kZllHFJPresel       ] = {"boostedCat","ZpTFJ","pTFJ","dPhiZHFJ"                                                       };
  /////////////////////////////
  // List of Samples
  vector<pair<TString,vhbbPlot::sampleType>> samples;
  
  bool isBatchMode= (batchSampleName!="");
  TString batchSuffix="";
  if(isBatchMode) {
    // Handle batch mode for Condor
    //multithread=false; // force single threading
    dataCardDir = dataCardDir + "/split/"; // write output in a subdirectory
    batchSuffix = "_"+batchSampleName; // add a suffix to the output with this sample's name
    samples.emplace_back(batchSampleName, batchSampleType);
  } else if(year==2016) {
    samples.emplace_back("LeptonPDSalad"                  , vhbbPlot::kData   );
    samples.emplace_back("SingleTop_tW"                   , vhbbPlot::kTop    );       
    samples.emplace_back("SingleTop_tbarW"                , vhbbPlot::kTop    );       
    samples.emplace_back("TTTo2L2Nu"                      , vhbbPlot::kTT     );       
    samples.emplace_back("WWTo2L2Nu"                      , vhbbPlot::kWW     );       
    samples.emplace_back("WZTo2L2Q"                       , vhbbPlot::kVZ     );       
    samples.emplace_back("ZZTo2L2Q"                       , vhbbPlot::kVZ     );       
    samples.emplace_back("ZJets_ht100to200"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_ht200to400"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_ht400to600"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_ht600to800"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_ht800to1200"              , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_ht1200to2500"             , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_ht2500toinf"              , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_pt50to100"                , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_pt100to250"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_pt250to400"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_pt400to650"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_pt650toinf"               , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_bHadrons"                 , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_bHadrons_pt100to200"      , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_bHadrons_pt200toinf"      , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_bQuarks"                  , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_bQuarks_pt100to200"       , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_bQuarks_pt200toinf"       , vhbbPlot::kZjets  );       
    samples.emplace_back("ZJets_m10"                      , vhbbPlot::kZjets  );       
    samples.emplace_back("ZllHbb_mH125"                   , vhbbPlot::kZH     );       
    samples.emplace_back("ggZllHbb_mH125"                 , vhbbPlot::kZH     );       
  } else if(year==2017) {
    //samples.emplace_back("DoubleEG"                       , vhbbPlot::kData   );
    //samples.emplace_back("DoubleMuon"                     , vhbbPlot::kData   );
    samples.emplace_back("LeptonPDSalad"                  , vhbbPlot::kData   );
    samples.emplace_back("Diboson_wz_CP5"                 , vhbbPlot::kVZ     );
    samples.emplace_back("Diboson_zz_CP5"                 , vhbbPlot::kVZ     );
    samples.emplace_back("Diboson_ww_CP5"                 , vhbbPlot::kWW     );
    samples.emplace_back("TTTo2L2Nu_CP5"                  , vhbbPlot::kTT     );
    samples.emplace_back("SingleTop_tW_CP5"               , vhbbPlot::kTop    );
    samples.emplace_back("SingleTop_tbarW_CP5"            , vhbbPlot::kTop    );
    // 0 files produced by CMS for DY1JetsToLL_M-50_LHEZpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8, work around it
    samples.emplace_back("Z1Jets_ZpT150to250_CP5"         , vhbbPlot::kZjets  );
    samples.emplace_back("Z1Jets_ZpT250to400_CP5"         , vhbbPlot::kZjets  );
    samples.emplace_back("Z1Jets_ZpT400toinf_CP5"         , vhbbPlot::kZjets  );
    samples.emplace_back("Z2Jets_ZpT50to150_CP5"          , vhbbPlot::kZjets  );
    samples.emplace_back("Z2Jets_ZpT150to250_CP5"         , vhbbPlot::kZjets  );
    samples.emplace_back("Z2Jets_ZpT250to400_CP5"         , vhbbPlot::kZjets  );
    samples.emplace_back("Z2Jets_ZpT400toinf_CP5"         , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_inclNLO_CP5"              , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_m4_ht70to100_CP5"         , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_m4_ht100to200_CP5"        , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_m4_ht200to400_CP5"        , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_m4_ht400to600_CP5"        , vhbbPlot::kZjets  );
    samples.emplace_back("ZJets_m4_ht600toinf_CP5"        , vhbbPlot::kZjets  );
    samples.emplace_back("ZllHbb_mH125"                   , vhbbPlot::kZH     );
    samples.emplace_back("ggZllHbb_mH125"                 , vhbbPlot::kZH     );
    // Ntuples we produce but do not currently use
    //samples.emplace_back("ZJets_ht100to200_CP5"           , vhbbPlot::kZjets  );
    //samples.emplace_back("ZJets_ht200to400_CP5"           , vhbbPlot::kZjets  );
    //samples.emplace_back("ZJets_ht400to600_CP5"           , vhbbPlot::kZjets  );
    //samples.emplace_back("ZJets_ht600to800_CP5"           , vhbbPlot::kZjets  );
    //samples.emplace_back("ZJets_ht800to1200_CP5"          , vhbbPlot::kZjets  );
    //samples.emplace_back("ZJets_ht1200to2500_CP5"         , vhbbPlot::kZjets  );
    //samples.emplace_back("ZJets_ht2500toinf_CP5"          , vhbbPlot::kZjets  );
  }
  if(multithread) std::random_shuffle(samples.begin(),samples.end());
  // End List of Samples
  /////////////////////////////
  
  // Load Shared Objects for ACLIC
  gSystem->Load("libPandaAnalysisFlat.so");
  
  assert(dataCardDir!="");
  if(dataCardDir!="") system(Form("mkdir -p MitVHBBAnalysis/datacards/%s",dataCardDir.Data()));
  ao.whichTriggers = (1<<pa::kSingleEleTrig) | (1<<pa::kDoubleEleTrig) | (1<<pa::kSingleMuTrig) | (1<<pa::kDoubleMuTrig) | (1<<pa::kEMuTrig);
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

  // Choice of the MVA variable type, binning, and the name
  // This can be different for each control region
  TString binZptSuffix="";
  if(binZpt>=0 && binZpt<nBinsZpt && 
    !(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel))
    binZptSuffix = Form("_ZptBin%d",binZpt);
  // Define the shape variable
  if(ao.MVAVarType==1) {
    // 1 - simple pT variable
    ao.MVAVarName="Higgs p_{T} classifier [GeV]";
    if(selection==kZllHSR) {
      ao.MVAbins={100,120,140,160,180,200,250,300,350};
      ao.MVAVarName="H(bb) pT";
      ao.shapeType="ptShape";
    } else if(selection==kZllHFJSR) {
      ao.MVAbins={250,300,350,400,450,500,550,600};
      ao.MVAVarName="H(bb) pT";
      ao.shapeType="ptShape";
    } else if(selection==kZllHLightFlavorCR) {
      ao.MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667};
      ao.MVAVarName="Subleading H(bb) CMVA";
      ao.shapeType="lesserCMVAShape";
    } else if(selection==kZllHHeavyFlavorCR || selection==kZllH2TopCR || selection==kZllHPresel) {
      ao.MVAbins={-0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000, 0.7333, 0.8667, 1.0000};
      ao.MVAVarName="Subleading H(bb) CMVA";
      ao.shapeType="lesserCMVAShape";
    } else if((selection>=kZllHLightFlavorFJCR && selection<kZllHFJSR) || selection==kZllHFJPresel) {
      if(selection==kZllHHeavyFlavorFJCR)
        ao.MVAbins={40,45,50,55,60,65,70,75,80};
      else
        ao.MVAbins={40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,160,180,200};
      ao.MVAVarName="Fatjet soft drop mass [GeV]";
      ao.shapeType="softDropMassShape";
    }
    // 2 - multiclass BDT in SR, subleading CMVA in CR (not implemented)
  } else if(ao.MVAVarType==3) {
    // 3 - normal BDT in SR, subleading CMVA in CR
    if(selection==kZllHLightFlavorCR) {
      ao.MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667};
      ao.MVAVarName="Subleading H(bb) CMVA";
      ao.shapeType="lesserCMVAShape";
    } else if(selection==kZllHHeavyFlavorCR || selection==kZllH2TopCR || selection==kZllHPresel) {
      ao.MVAbins={-0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000, 0.7333, 0.8667, 1.0000};
      ao.MVAVarName="Subleading H(bb) CMVA";
      ao.shapeType="lesserCMVAShape";
    } else if(selection==kZllHSR) {
      if(binZpt==0)
        ao.MVAbins={-1,0,0.4,0.5,0.6,0.7,0.8,1};
      else
        ao.MVAbins={-1,0,0.4,0.5,0.6,0.7,0.8,1};
      ao.MVAVarName="BDT Output";
      ao.shapeType="singleClassBDTShape"; 
    } else if((selection>=kZllHLightFlavorFJCR && selection<kZllHFJSR) || selection==kZllHFJPresel) {
      if(selection==kZllHHeavyFlavorFJCR)
        ao.MVAbins={40,45,50,55,60,65,70,75,80};
      else
        ao.MVAbins={40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,160,180,200};
      ao.MVAVarName="Fatjet soft drop mass [GeV]";
      ao.shapeType="softDropMassShape";
    } else if(selection==kZllHFJSR) {
      ao.MVAbins={-0.8,-0.4,0,0.4,0.8};
      ao.MVAVarName="BDT Output";
      ao.shapeType="singleClassBDTShape"; 
    }
  } else throw std::runtime_error("bad ao.MVAVarType");
  
  // Declare histograms for plotting
  printf("Building plotting histograms, please wait...\n");
  ao.xmin.resize(nPlots);
  ao.xmax.resize(nPlots);
  ao.nbins.resize(nPlots);
  ao.histoNames.resize(nPlots); 
  ao.histoTitles.resize(nPlots);
  { int p=0;
    ao.histoNames[p]="MVAVar"                  ; ao.histoTitles[p]=""                      ;                                                p++;
    ao.histoNames[p]="lepton1Pt"               ; ao.histoTitles[p]="Lepton 1 p_{T} [GeV]"  ; ao.nbins[p]=  23; ao.xmin[p]=    20; ao.xmax[p]=   250; p++; 
    ao.histoNames[p]="lepton2Pt"               ; ao.histoTitles[p]="Lepton 2 p_{T} [GeV]"  ; ao.nbins[p]=  23; ao.xmin[p]=    20; ao.xmax[p]=   250; p++; 
    ao.histoNames[p]="lepton1Eta"              ; ao.histoTitles[p]="Lepton 1 #eta"         ; ao.nbins[p]=  25; ao.xmin[p]=  -2.5; ao.xmax[p]=   2.5; p++; 
    ao.histoNames[p]="lepton2Eta"              ; ao.histoTitles[p]="Lepton 2 #eta"         ; ao.nbins[p]=  25; ao.xmin[p]=  -2.5; ao.xmax[p]=   2.5; p++; 
    ao.histoNames[p]="ZBosonPt"                ; ao.histoTitles[p]="Z boson pT [GeV]"      ; ao.nbins[p]=  50; ao.xmin[p]=     0; ao.xmax[p]=   500; p++; 
    ao.histoNames[p]="ZBosonEta"               ; ao.histoTitles[p]="Z boson #eta"          ; ao.nbins[p]=  25; ao.xmin[p]=  -2.5; ao.xmax[p]=   2.5; p++; 
    ao.histoNames[p]="ZBosonPhi"               ; ao.histoTitles[p]="Z boson phi"           ; ao.nbins[p]=  32; ao.xmin[p]=     0; ao.xmax[p]= 3.142; p++; 
    ao.histoNames[p]="ZBosonM"                 ; ao.histoTitles[p]="Z boson mass"          ; ao.nbins[p]=  60; ao.xmin[p]=     0; ao.xmax[p]=   120; p++; 
    ao.histoNames[p]="ZBosonLep1CosThetaCS"    ; ao.histoTitles[p]="Z(ll) cos#theta^{CS}"  ; ao.nbins[p]=  20; ao.xmin[p]=    -1; ao.xmax[p]=    1.; p++; 
    ao.histoNames[p]="bdtValue"                ; ao.histoTitles[p]="BDT Output"            ; ao.nbins[p]=  40; ao.xmin[p]=    -1; ao.xmax[p]=    1.; p++; 
    if(selection>=kZllHLightFlavorFJCR && selection<=kZllHFJPresel) {
      // fatjet only plots
      ao.histoNames[p]="mSD"                     ; ao.histoTitles[p]="Fatjet mSD [GeV]"      ; ao.nbins[p]=  32; ao.xmin[p]=    40; ao.xmax[p]=   200; p++; 
      ao.histoNames[p]="pTFJ"                    ; ao.histoTitles[p]="Fatjet pT [GeV]"       ; ao.nbins[p]=  25; ao.xmin[p]=   250; ao.xmax[p]=   600; p++; 
      ao.histoNames[p]="Tau21SD"                 ; ao.histoTitles[p]="#tau_{2}/#tau_{1} SD"  ; ao.nbins[p]=  20; ao.xmin[p]=     0; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="doubleB"                 ; ao.histoTitles[p]="Fatjet double b-tag"   ; ao.nbins[p]=  40; ao.xmin[p]=   -1.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="ZBosonLep1CosThetaStarFJ"; ao.histoTitles[p]="cos#theta* Z(ll)+FJ"   ; ao.nbins[p]=  20; ao.xmin[p]=    -1; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="deltaEtaZHFJ"            ; ao.histoTitles[p]="|#Delta#eta(Z,FJ)|"    ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="deltaPhiZHFJ"            ; ao.histoTitles[p]="#Delta#phi(Z,FJ) [Rad]"; ao.nbins[p]=  20; ao.xmin[p]= 1.571; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="mTZHFJ"                  ; ao.histoTitles[p]="m_{T} (Z+FJ) [GeV]"    ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=  200.; p++; 
      ao.histoNames[p]="dRL1L2"                  ; ao.histoTitles[p]="#DeltaR(l1,l2)"        ; ao.nbins[p]=  20; ao.xmin[p]=    0.; ao.xmax[p]=    5.; p++; 
      ao.histoNames[p]="ptBalanceZHFJ"           ; ao.histoTitles[p]="|FJ pT / Z pT|"        ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=    3.; p++; 
      ao.histoNames[p]="ptBalanceL1L2"           ; ao.histoTitles[p]="p_{T}^{l1}/p_{T}^{l2}" ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=    3.; p++; 
      ao.histoNames[p]="nIsojet"                 ; ao.histoTitles[p]="N isojets"             ; ao.nbins[p]=   8; ao.xmin[p]=    0.; ao.xmax[p]=    8.; p++; 
      ao.histoNames[p]="isojetNBtags"            ; ao.histoTitles[p]="N isojet b-tags"       ; ao.nbins[p]=   4; ao.xmin[p]=    0.; ao.xmax[p]=    4.; p++; 
    } else {
      ao.histoNames[p]="Mjj"                     ; ao.histoTitles[p]="Dijet mass [GeV]"      ; ao.nbins[p]=  25; ao.xmin[p]=     0; ao.xmax[p]=   250; p++; 
      ao.histoNames[p]="pTjj"                    ; ao.histoTitles[p]="Dijet pT [GeV]"        ; ao.nbins[p]=  18; ao.xmin[p]=    50; ao.xmax[p]=   350; p++; 
      ao.histoNames[p]="bjet1Pt"                 ; ao.histoTitles[p]="B-jet 1 pT [GeV]"      ; ao.nbins[p]=  38; ao.xmin[p]=    20; ao.xmax[p]=   400; p++; 
      ao.histoNames[p]="bjet2Pt"                 ; ao.histoTitles[p]="B-jet 2 pT [GeV]"      ; ao.nbins[p]=  38; ao.xmin[p]=    20; ao.xmax[p]=   400; p++; 
      ao.histoNames[p]="bjet1btag"               ; ao.histoTitles[p]="B-jet 1 btag"          ; ao.nbins[p]=  40; ao.xmin[p]=   -1.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="bjet2btag"               ; ao.histoTitles[p]="B-jet 2 btag"          ; ao.nbins[p]=  40; ao.xmin[p]=   -1.; ao.xmax[p]=    1.; p++; 
      ao.histoNames[p]="nJet"                    ; ao.histoTitles[p]="N central AK4CHS jets" ; ao.nbins[p]=   8; ao.xmin[p]=    0.; ao.xmax[p]=    8.; p++; 
      ao.histoNames[p]="deltaPhiZH"              ; ao.histoTitles[p]="#Delta#phi(Z,H) [Rad]" ; ao.nbins[p]=  20; ao.xmin[p]= 1.571; ao.xmax[p]= 3.142; p++; 
      ao.histoNames[p]="ptBalanceZH"             ; ao.histoTitles[p]="|H pT / Z pT|"         ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=    3.; p++; 
      ao.histoNames[p]="sumEtSoft1"              ; ao.histoTitles[p]="#sum E_{T}(soft 1)"    ; ao.nbins[p]=  30; ao.xmin[p]=    0.; ao.xmax[p]=  300.; p++; 
      ao.histoNames[p]="nSoft2"                  ; ao.histoTitles[p]="N^{soft}_{2}"          ; ao.nbins[p]=  25; ao.xmin[p]=    0.; ao.xmax[p]=   25.; p++; 
      ao.histoNames[p]="nSoft5"                  ; ao.histoTitles[p]="N^{soft}_{5}"          ; ao.nbins[p]=  12; ao.xmin[p]=    0.; ao.xmax[p]=   12.; p++; 
      ao.histoNames[p]="nSoft10"                 ; ao.histoTitles[p]="N^{soft}_{10}"         ; ao.nbins[p]=   8; ao.xmin[p]=    0.; ao.xmax[p]=    8.; p++; 
      ao.histoNames[p]="ZBosonLep1CosThetaStar"  ; ao.histoTitles[p]="cos#theta* Z(ll)+jj"   ; ao.nbins[p]=  20; ao.xmin[p]=    -1; ao.xmax[p]=    1.; p++; 
    }
  }
  
  for(unsigned lep=0; lep<nLepSel; lep++) 
  for(int p=0; p<nPlots; p++) {
    if(ao.histoNames[p]=="") continue;
    for(unsigned char ic=0; ic<nPlotCategories; ic++) {
      if(ao.histoNames[p]=="MVAVar") {
        ao.histos[lep][p][ic] = new TH1F(
          Form("histo%d_Z%sH%s_%s",
            ic,
            leptonStrings[lep].Data(),
            selectionNames[(int)selection].Data(),
            ao.histoNames[p].Data()
          ), ao.MVAVarName, ao.MVAbins.size()-1, ao.MVAbins.data());
      } else {
        ao.histos[lep][p][ic] = new TH1F(
          Form("histo%d_Z%sH%s_%s",
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
    ao.histo_Baseline     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s"                , plotBaseNames[ic].Data()));
    if(ic<kPlotVZbb) continue;
    ao.histo_pileupUp     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_pileupUp"       , plotBaseNames[ic].Data()));
    ao.histo_pileupDown   [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_pileupDown"     , plotBaseNames[ic].Data()));
    ao.histo_VHCorrUp     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VHCorrUp"       , plotBaseNames[ic].Data()));
    ao.histo_VHCorrDown   [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VHCorrDown"     , plotBaseNames[ic].Data()));
    ao.histo_QCDr1f2      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr1f2"        , plotBaseNames[ic].Data()));
    ao.histo_QCDr1f5      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr1f5"        , plotBaseNames[ic].Data()));
    ao.histo_QCDr2f1      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr2f1"        , plotBaseNames[ic].Data()));
    ao.histo_QCDr2f2      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr2f2"        , plotBaseNames[ic].Data()));
    ao.histo_QCDr5f1      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr5f1"        , plotBaseNames[ic].Data()));
    ao.histo_QCDr5f5      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDr5f5"        , plotBaseNames[ic].Data()));
    ao.histo_QCDScaleUp   [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDScale%sUp"   , plotBaseNames[ic].Data(),plotBaseNames[ic].Data()));
    ao.histo_QCDScaleDown [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_QCDScale%sDown" , plotBaseNames[ic].Data(),plotBaseNames[ic].Data()));
    ao.histo_eleSFUp      [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_eleSFUp"        , plotBaseNames[ic].Data()));
    ao.histo_eleSFDown    [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_eleSFDown"      , plotBaseNames[ic].Data()));
    ao.histo_muSFUp       [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_muSFUp"         , plotBaseNames[ic].Data()));
    ao.histo_muSFDown     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_muSFDown"       , plotBaseNames[ic].Data()));
    ao.histo_VGluUp       [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VGluUp"         , plotBaseNames[ic].Data()));
    ao.histo_VGluDown     [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_VGluDown"       , plotBaseNames[ic].Data()));
    ao.histo_doubleBUp    [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_doubleBUp"      , plotBaseNames[ic].Data()));
    ao.histo_doubleBDown  [lep][ic] = (TH1F*)ao.histos[lep][0][ic]->Clone(Form("histo_%s_doubleBDown"    , plotBaseNames[ic].Data()));
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
          Form("histo_%s_CMS_VH_btag_pt%d_eta%d_%s",
          plotBaseNames[ic].Data(),iPt,iEta,btagShiftName(shift))
        );
    }
  }
  printf("Done building uncertainty histograms\n");


  ////////////////////////////////////////////////////////////////////////
  // Load corrections to apply offline
  // CMVA reweighting for 2016, DeepCSV reweighting for 2017
  std::vector<std::string> btagSystNames;
  if(ao.year==2016) {
    ao.cmvaReweighter = new CSVHelper(
      "PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root", 
      "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root", 
       5);
  } else if(ao.year==2017) {
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
      "PandaAnalysis/data/csv/DeepCSV_94XSF_V2_B_F.csv");
    ao.deepcsvSFs->load(*(ao.deepcsvCalib), BTagEntry::FLAV_B, "iterativeFit");
    ao.deepcsvSFs->load(*(ao.deepcsvCalib), BTagEntry::FLAV_C, "iterativeFit");
    ao.deepcsvSFs->load(*(ao.deepcsvCalib), BTagEntry::FLAV_UDSG, "iterativeFit");
    
  }
  ao.ZjetsEWKCorr = EWKCorrPars(kZjets);
  if(!useHtBinnedVJetsKFactor) {
    TFile *kfactorsFile = TFile::Open("PandaAnalysis/data/higgs/hbb_kfactors.root","read");
    ao.kfactors_ZJets = (TH2D*) kfactorsFile->Get("h_ZJets");
    ao.kfactors_ZJets->SetDirectory(0);
    kfactorsFile->Close();
  }

  // Done loading offline corrections
  ////////////////////////////////////////////////////////////////////////
  // Setup MVA training tree if applicable (Not implemented yet for boosted)
  if(selection==kZllHSR) {
    system(Form("mkdir -p MitVHBBAnalysis/mva/%s",dataCardDir.Data()));
    ao.mvaFile = new TFile(Form("MitVHBBAnalysis/mva/%s/ZllHSR_mvaTree%s%s.root",
      dataCardDir.Data(),binZptSuffix.Data(),batchSuffix.Data()),"recreate");
    ao.mvaTree = new TTree("mvaTree","mvaTree");
    ao.mvaTree->Branch("sumEtSoft1"  , &ao.mva_sumEtSoft1  ); 
    ao.mvaTree->Branch("nSoft2"      , &ao.mva_nSoft2      ); 
    ao.mvaTree->Branch("nSoft5"      , &ao.mva_nSoft5      ); 
    ao.mvaTree->Branch("nSoft10"     , &ao.mva_nSoft10     ); 
    ao.mvaTree->Branch("bjet1Pt"     , &ao.mva_bjet1Pt     ); 
    ao.mvaTree->Branch("bjet2Pt"     , &ao.mva_bjet2Pt     ); 
    ao.mvaTree->Branch("bjet1btag"   , &ao.mva_bjet1btag   ); 
    ao.mvaTree->Branch("bjet2btag"   , &ao.mva_bjet2btag   ); 
    ao.mvaTree->Branch("ZBosonPt"    , &ao.mva_ZBosonPt    ); 
    ao.mvaTree->Branch("ZBosonM"     , &ao.mva_ZBosonM     ); 
    ao.mvaTree->Branch("CosThetaCS"  , &ao.mva_CosThetaCS  ); 
    ao.mvaTree->Branch("CosThetaStar", &ao.mva_CosThetaStar); 
    ao.mvaTree->Branch("hbbpt"       , &ao.mva_hbbpt       ); 
    ao.mvaTree->Branch("hbbm"        , &ao.mva_hbbm        ); 
    ao.mvaTree->Branch("dPhiZH"      , &ao.mva_dPhiZH      ); 
    ao.mvaTree->Branch("ptBalanceZH" , &ao.mva_ptBalanceZH ); 
    ao.mvaTree->Branch("dRBjets"     , &ao.mva_dRBjets     ); 
    ao.mvaTree->Branch("dEtaBjets"   , &ao.mva_dEtaBjets   ); 
    ao.mvaTree->Branch("dRZH"        , &ao.mva_dRZH        ); 
    ao.mvaTree->Branch("dEtaZH"      , &ao.mva_dEtaZH      ); 
    ao.mvaTree->Branch("nAddJet"     , &ao.mva_nAddJet     ); 
    ao.mvaTree->Branch("weight"      , &ao.mva_weight      ); 
    ao.mvaTree->Branch("category"    , &ao.mva_category    ); 
    ao.mvaTree->Branch("eventNumber" , &ao.mva_eventNumber ); 
  } else if(selection==kZllHFJSR) {
    system(Form("mkdir -p MitVHBBAnalysis/mva/%s",dataCardDir.Data()));
    ao.mvaFile = new TFile(Form("MitVHBBAnalysis/mva/%s/ZllHFJSR_mvaTree%s.root",
      dataCardDir.Data(),batchSuffix.Data()),"recreate");
    ao.mvaTree = new TTree("mvaTree","mvaTree");
    ao.mvaTree->Branch("nIsojet"       , &ao.mva_nIsojet        ); 
    ao.mvaTree->Branch("MSD"           , &ao.mva_MSD            ); 
    ao.mvaTree->Branch("fjPt"          , &ao.mva_fjPt           ); 
    ao.mvaTree->Branch("Tau21SD"       , &ao.mva_Tau21SD        ); 
    ao.mvaTree->Branch("ptBalanceZHFJ" , &ao.mva_ptBalanceZHFJ  ); 
    ao.mvaTree->Branch("dEtaZHFJ"      , &ao.mva_dEtaZHFJ       ); 
    ao.mvaTree->Branch("dPhiZHFJ"      , &ao.mva_dPhiZHFJ       ); 
    ao.mvaTree->Branch("mTZHFJ"        , &ao.mva_mTZHFJ         ); 
    ao.mvaTree->Branch("ptBalanceL1L2" , &ao.mva_ptBalanceL1L2  ); 
    ao.mvaTree->Branch("dRL1L2"        , &ao.mva_dRL1L2         ); 
    ao.mvaTree->Branch("ZBosonPt"      , &ao.mva_ZBosonPt       ); 
    ao.mvaTree->Branch("CosThetaCS"    , &ao.mva_CosThetaCS     ); 
    ao.mvaTree->Branch("lepton1Pt"     , &ao.mva_lepton1Pt      ); 
    ao.mvaTree->Branch("lepton2Pt"     , &ao.mva_lepton2Pt      ); 
    ao.mvaTree->Branch("lepton1Eta"    , &ao.mva_lepton1Eta     ); 
    ao.mvaTree->Branch("lepton2Eta"    , &ao.mva_lepton2Eta     ); 
    ao.mvaTree->Branch("deltaM"        , &ao.mva_deltaM         ); 
    ao.mvaTree->Branch("weight"        , &ao.mva_weight         ); 
    ao.mvaTree->Branch("category"      , &ao.mva_category       ); 
    ao.mvaTree->Branch("eventNumber"   , &ao.mva_eventNumber    ); 
  }
  
  // Instantiate TMVA reader
  if(MVAVarType>1) {
    TString bdtWeights="";
    if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) 
      bdtWeights = "MitVHBBAnalysis/weights/bdt_BDT_singleClass_boosted_may24.weights.xml";
    else if(ao.binZpt==0)
      bdtWeights = (ao.year==2016)?
        "MitVHBBAnalysis/weights/bdt_BDT_singleClass_resolved_ZptBin0.weights.xml":
        "MitVHBBAnalysis/weights/bdt_BDT_singleClass_resolved_ZptBin0_2017.weights.xml";
    else if(ao.binZpt==1) 
      bdtWeights = (ao.year==2016)?
        "MitVHBBAnalysis/weights/bdt_BDT_singleClass_resolved_ZptBin1.weights.xml":
        "MitVHBBAnalysis/weights/bdt_BDT_singleClass_resolved_ZptBin1_2017.weights.xml";
    if(bdtWeights!="") for(unsigned nThread=0; nThread < (multithread? nThreads:1); nThread++) {
      TMVA::Reader *theReader = new TMVA::Reader("Silent");
      // This object is never deleted, which is a small memory leak,
      // but the TMVA Reader destructor has issues
      
      // Vars are hardcoded for now, could make it more general if we care
      if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) {
        theReader->AddVariable("nIsojet"       , &ao.mvaInputs[nThread][ 0]); 
        theReader->AddVariable("fjPt"          , &ao.mvaInputs[nThread][ 1]); 
        theReader->AddVariable("MSD"           , &ao.mvaInputs[nThread][ 2]); 
        theReader->AddVariable("Tau21SD"       , &ao.mvaInputs[nThread][ 3]); 
        theReader->AddVariable("ptBalanceZHFJ" , &ao.mvaInputs[nThread][ 4]); 
        theReader->AddVariable("dEtaZHFJ"      , &ao.mvaInputs[nThread][ 5]); 
        theReader->AddVariable("dPhiZHFJ"      , &ao.mvaInputs[nThread][ 6]); 
        theReader->AddVariable("mTZHFJ"        , &ao.mvaInputs[nThread][ 7]); 
        theReader->AddVariable("ptBalanceL1L2" , &ao.mvaInputs[nThread][ 8]); 
        theReader->AddVariable("dRL1L2"        , &ao.mvaInputs[nThread][ 9]); 
        theReader->AddVariable("ZBosonPt"      , &ao.mvaInputs[nThread][10]); 
        theReader->AddVariable("lepton1Pt"     , &ao.mvaInputs[nThread][11]); 
        theReader->AddVariable("lepton2Pt"     , &ao.mvaInputs[nThread][12]); 
        theReader->AddVariable("deltaM"        , &ao.mvaInputs[nThread][13]); 
        theReader->AddVariable("CosThetaCS"    , &ao.mvaInputs[nThread][14]);
      } else {
        theReader->AddVariable("sumEtSoft1"  , &ao.mvaInputs[nThread][ 0]);
        theReader->AddVariable("bjet1Pt"     , &ao.mvaInputs[nThread][ 1]);
        theReader->AddVariable("bjet2Pt"     , &ao.mvaInputs[nThread][ 2]);
        theReader->AddVariable("bjet1btag"   , &ao.mvaInputs[nThread][ 3]);
        theReader->AddVariable("bjet2btag"   , &ao.mvaInputs[nThread][ 4]);
        theReader->AddVariable("ZBosonPt"    , &ao.mvaInputs[nThread][ 5]);
        theReader->AddVariable("ZBosonM"     , &ao.mvaInputs[nThread][ 6]);
        theReader->AddVariable("CosThetaCS"  , &ao.mvaInputs[nThread][ 7]);
        theReader->AddVariable("CosThetaStar", &ao.mvaInputs[nThread][ 8]);
        theReader->AddVariable("hbbpt"       , &ao.mvaInputs[nThread][ 9]);
        theReader->AddVariable("hbbm"        , &ao.mvaInputs[nThread][10]);
        theReader->AddVariable("dPhiZH"      , &ao.mvaInputs[nThread][11]);
        theReader->AddVariable("ptBalanceZH" , &ao.mvaInputs[nThread][12]);
        theReader->AddVariable("dRBjets"     , &ao.mvaInputs[nThread][13]);
        theReader->AddVariable("dEtaBjets"   , &ao.mvaInputs[nThread][14]);
        theReader->AddVariable("nAddJet"     , &ao.mvaInputs[nThread][15]);
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
  if(selection==kZllHSR || selection==kZllHFJSR) {
    ao.mvaFile->cd();
    ao.mvaTree->Write();
    ao.mvaFile->Close();
  }
  if(ao.deepcsvSFs) delete ao.deepcsvSFs;
  if(ao.deepcsvCalib) delete ao.deepcsvCalib;
  if(ao.cmvaReweighter) delete ao.cmvaReweighter;

  // renormalize V+Glu shapes to the nominal norm
  if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel)
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
    if(leptonStrings[lep]=="em" && !(selection==kZllHSR||selection==kZllHFJSR)) continue; // avoiding unnecessary histograms
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
      if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) {
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
  }
  
  // Write shape histograms to file
  for(unsigned lep=0; lep<nLepSel; lep++) {
    // avoiding unnecessary histograms
    if(leptonStrings[lep]=="em" && !(selection==kZllHSR||selection==kZllHFJSR)) continue; 

    char outFileDatacardsName[200];
    sprintf(
      outFileDatacardsName,
      "MitVHBBAnalysis/datacards/%s/datacard_Z%sH%s%s%s.root",
      dataCardDir.Data(),
      leptonStrings[lep].Data(),
      selectionNames[selection].Data(),
      binZptSuffix.Data(),
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
        // Special symmetrization procedure
        TString shiftName(jesName(static_cast<shiftjes>(iJES)).Data());
        if(shiftName.EndsWith("Down")){
          for(int nb=1; nb<=ao.histo_Baseline[lep][ic]->GetNbinsX(); nb++){
            double diff = ao.histo_jes[iJES-1][lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb);
            double mean = ao.histo_Baseline[lep][ic]->GetBinContent(nb);
            ao.histo_jes[iJES][lep][ic]->SetBinContent(nb,mean-diff);
          }
        }
        ao.histo_jes[iJES][lep][ic]->Write();
      }
      for(unsigned iPt=0; iPt<5; iPt++)
      for(unsigned iEta=0; iEta<3; iEta++)
      for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
        GeneralTree::csvShift shift = gt.csvShifts[iShift];
        if (shift==GeneralTree::csvCent) continue;
        // Special symmetrization procedure
        TString shiftName(btagShiftName(shift));
        if(shiftName.EndsWith("Down")){
          for(int nb=1; nb<=ao.histo_Baseline[lep][ic]->GetNbinsX(); nb++){
            double diff = ao.histo_btag[iShift-1][iPt][iEta][lep][ic]->GetBinContent(nb)-ao.histo_Baseline[lep][ic]->GetBinContent(nb);
            double mean = ao.histo_Baseline[lep][ic]->GetBinContent(nb);
            ao.histo_btag[iShift][iPt][iEta][lep][ic]->SetBinContent(nb,mean-diff);
          }
        }
        ao.histo_btag[iShift][iPt][iEta][lep][ic]->Write();
      }
      if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) {
        ao.histo_VGluUp     [lep][ic]->Write();
        ao.histo_VGluDown   [lep][ic]->Write();
        ao.histo_doubleBUp  [lep][ic]->Write();
        ao.histo_doubleBDown[lep][ic]->Write();
      }
    }
    outFileDatacards->Close();
  }

  // Writing datacards - need to move this to separate function
  if(!isBatchMode)
    writeDatacards(ao, dataCardDir);

  // Write plots
  char regionName[128];
  for(unsigned lep=0; lep<nLepSel; lep++) {
    sprintf(regionName, "Z%sH%s",leptonStrings[lep].Data(),selectionNames[selection].Data());
    TString plotFileName = Form("MitVHBBAnalysis/datacards/%s/plots_%s%s%s.root",dataCardDir.Data(),regionName,binZptSuffix.Data(),batchSuffix.Data());
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
  
  // Sample properties
  // Only use events with HT<100 for NLO pt binned samples in 2016
  bool isLowMassZjets = sampleName.Contains("ZJets_m10") || sampleName.Contains("ZJets_m4");
  bool isBQuarkEnriched = sampleName.Contains("bQuarks");
  bool isBHadronEnriched = sampleName.Contains("bHadrons");
  bool isInclusiveZjets = type==vhbbPlot::kZjets && !sampleName.Contains("_ht") && !sampleName.Contains("_pt");
  bool useNPNLOLookup = (ao.year==2017 && sampleName=="ZJets_inclNLO_CP5");
  bool isV12jets = sampleName.Contains("Z1Jets") || sampleName.Contains("Z2Jets");
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
  gt.useCMVA        = true;
  gt.is_breg        = false;
  // Branches not in GeneralTree;
  std::map<TString, void*> extraAddresses;
  float normalizedWeight; unsigned char npnlo;
  extraAddresses["normalizedWeight"] = &normalizedWeight;
  if(useNPNLOLookup)
    extraAddresses["npnlo"] = &npnlo;
    
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
      if(isNLOZjets && !isLowMassZjets) { 
        // for low pT M>50 NLO Z+jets in 2016
        bLoad(b["lheHT"],ientry);
        if(gt.lheHT>=100) continue;
      }
      // B-enriched sample stitching
      // Let the b-enriched samples eat 90% of the XS where there are b quarks or status 2 b hadrons
      if(type==kZjets) {
        bLoad(b["trueGenBosonPt"],ientry); // LHE Z boson pT
        bLoad(b["nStatus2BHadrons"],ientry); // number of B hadrons at matrix element level
        bLoad(b["nB"],ientry); // number of B quarks
        if(isInclusiveZjets && (isBQuarkEnriched||isBHadronEnriched) && gt.trueGenBosonPt>=100) continue;
        bool hasBQuarks  = gt.nB > 0;
        bool hasBHadrons = gt.nStatus2BHadrons>0 && gt.nB==0;
        // Orthogonalize        
        if(isBQuarkEnriched && !hasBQuarks) continue;
        if(isBHadronEnriched && !hasBHadrons) continue;
        // Downweight
        if(hasBQuarks) {
          if(isBQuarkEnriched) stitchWeight = 0.9;
          else                 stitchWeight = 0.1;
        } else if(hasBHadrons) {
          if(isBHadronEnriched) stitchWeight = 0.9;
          else                  stitchWeight = 0.1;
        }
      }
    } else if(ao.year==2017) {
      if(useNPNLOLookup) {
        // Here, we perform the lookup of the NPNLO (number of NLO partons) for inclusive Z+jets
        // This is not a quantity in Panda, we must use a lookup table. Can potentially add this to the ntuples as an extra tree later
        bLoad(b["npnlo"],ientry);
        bLoad(b["eventNumber"],ientry);
        bLoad(b["trueGenBosonPt"],ientry); 
        if(npnlo==255) {
          printf("WARNING: NPNLO=255 for eventtNumber %llu\n", gt.eventNumber);
          continue;
        }
        if(npnlo==1 && gt.trueGenBosonPt<150) stitchWeight=1; // HACK because of missing files!
        else if(npnlo==1 || npnlo==2) stitchWeight=0.2;
        else stitchWeight=1;
      } else if(isV12jets) {
        stitchWeight=0.8;
      }
    }

    //////////////////////
    // Clear variables
    typeLepSel=99; // 0: Z(mm), 1: Z(ee), 2: e-mu

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
    if(gt.nLooseLep!=2) continue; 
    if(ao.debug) printf("  Passed lepton multiplicity\n");

    // Trigger
    bLoad(b["trigger"],ientry);
    bool passTrigger = (gt.trigger & ao.whichTriggers) !=0;
    if(!passTrigger) continue;
    if(ao.debug) printf("  Passed trigger\n");

    // Lepton ID and isolation
    bLoad(b["nLooseElectron"],ientry);
    bLoad(b["nTightElectron"],ientry);
    bLoad(b["nLooseMuon"],ientry);
    bLoad(b["nTightMuon"],ientry);
    bLoad(b["muonSelBit"],ientry);
    bLoad(b["muonPt"],ientry);
    bLoad(b["muonPdgId"],ientry);
    bLoad(b["electronSelBit"],ientry);
    bLoad(b["electronPt"],ientry);
    bLoad(b["electronPdgId"],ientry);

    if(gt.nLooseMuon==2 && gt.nLooseElectron==0 &&
      gt.muonPt[0]>25 && gt.muonPt[1]>10 && 
      gt.muonPdgId[0]+gt.muonPdgId[1]==0 &&
      (gt.muonSelBit[0] & pa::kTight)!=0 &&
      (gt.muonSelBit[1] & pa::kTight)!=0
    ) typeLepSel=0;
    else if(gt.nLooseElectron==2 && gt.nLooseMuon==0 &&
      gt.electronPt[0]>25 && gt.electronPt[1]>15 && 
      gt.electronPdgId[0]+gt.electronPdgId[1]==0 &&
      (gt.electronSelBit[0] & pa::kEleMvaWP90)!=0 &&
      (gt.electronSelBit[1] & pa::kEleMvaWP90)!=0
    ) typeLepSel=1;
    else if(gt.nLooseMuon==1 && gt.nLooseElectron==1 &&
      ((gt.muonPt[0]>25 && gt.electronPt[0]>15) || (gt.muonPt[0]>15 && gt.electronPt[0]>25)) && 
      gt.muonPdgId[0]*gt.electronPdgId[0]<0 &&
      (gt.muonSelBit[0]     & pa::kTight     )!=0 &&
      (gt.electronSelBit[0] & pa::kEleMvaWP90)!=0
    ) typeLepSel=2;
    if(typeLepSel!=0 && typeLepSel!=1 && typeLepSel!=2) continue;
    if(ao.debug) printf("  Passed lepton ID/iso multiplicity\n");
    
    // Z boson basics
    bLoad(b["ZBosonPt"],ientry);
    bLoad(b["ZBosonEta"],ientry);
    bLoad(b["ZBosonPhi"],ientry);
    bLoad(b["ZBosonM"],ientry);
    // Apply the ZpT window for this pt bin if doing resolved, provided
    // we have chosen to do so by choosing a non-negative value of binZpt
    float preselMinZpt=30, preselMaxZpt=9999;
    if(
      ao.binZpt>=0 && ao.binZpt<nBinsZpt && 
      !(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel)
    ) {
      preselMinZpt = binsZpt[ao.binZpt];
      preselMaxZpt = binsZpt[ao.binZpt+1];
    }
    if(gt.ZBosonPt<preselMinZpt || gt.ZBosonPt>=preselMaxZpt) continue;
    if(gt.ZBosonM<10) continue;
    if(ao.debug) printf("  Passed Z boson reconstruction\n");

    // Lepton kinematics
    float lepton1Pt,lepton1Eta,lepton1Phi,lepton1RelIso,lepton1D0,lepton1DZ,
          lepton2Pt,lepton2Eta,lepton2Phi,lepton2RelIso,lepton2D0,lepton2DZ;
    if (typeLepSel==0) {
      bLoad(b["muonEta"],ientry);
      bLoad(b["muonPhi"],ientry);
      bLoad(b["muonD0"],ientry);
      bLoad(b["muonDZ"],ientry);
      bLoad(b["muonCombIso"],ientry);
      bLoad(b["muonPdgId"],ientry);
      lepton1Pt     = gt.muonPt[0];
      lepton1Eta    = gt.muonEta[0];
      lepton1Phi    = gt.muonPhi[0]; 
      lepton1D0     = gt.muonD0[0];
      lepton1DZ     = gt.muonDZ[0];
      lepton1RelIso = gt.muonCombIso[0]/gt.muonPt[0];
      lepton2Pt     = gt.muonPt[1];
      lepton2Eta    = gt.muonEta[1];
      lepton2Phi    = gt.muonPhi[1]; 
      lepton2D0     = gt.muonD0[1];
      lepton2DZ     = gt.muonDZ[1];
      lepton2RelIso = gt.muonCombIso[1]/gt.muonPt[1];
    } else if(typeLepSel==1) {
      bLoad(b["electronEta"],ientry);
      bLoad(b["electronPhi"],ientry);
      bLoad(b["electronD0"],ientry);
      bLoad(b["electronDZ"],ientry);
      bLoad(b["electronCombIso"],ientry);
      bLoad(b["electronPdgId"],ientry);
      lepton1Pt     = gt.electronPt[0]; 
      lepton1Eta    = gt.electronEta[0];
      lepton1Phi    = gt.electronPhi[0];
      lepton1RelIso = gt.electronCombIso[0]/gt.electronPt[0];
      lepton1D0     = gt.electronD0[0];
      lepton1DZ     = gt.electronDZ[0];
      lepton2Pt     = gt.electronPt[1]; 
      lepton2Eta    = gt.electronEta[1];
      lepton2Phi    = gt.electronPhi[1];
      lepton2RelIso = gt.electronCombIso[1]/gt.electronPt[1];
      lepton2D0     = gt.electronD0[1];
      lepton2DZ     = gt.electronDZ[1];
    } else if (typeLepSel==2) {
      bLoad(b["muonEta"],ientry);
      bLoad(b["muonPhi"],ientry);
      bLoad(b["muonD0"],ientry);
      bLoad(b["muonDZ"],ientry);
      bLoad(b["muonCombIso"],ientry);
      bLoad(b["muonPdgId"],ientry);
      bLoad(b["electronEta"],ientry);
      bLoad(b["electronPhi"],ientry);
      bLoad(b["electronD0"],ientry);
      bLoad(b["electronDZ"],ientry);
      bLoad(b["electronCombIso"],ientry);
      bLoad(b["electronPdgId"],ientry);
      lepton1Pt     = gt.muonPt[0];
      lepton1Eta    = gt.muonEta[0];
      lepton1Phi    = gt.muonPhi[0]; 
      lepton1D0     = gt.muonD0[0];
      lepton1DZ     = gt.muonDZ[0];
      lepton1RelIso = gt.muonCombIso[0]/gt.muonPt[0];
      lepton2Pt     = gt.electronPt[0]; 
      lepton2Eta    = gt.electronEta[0];
      lepton2Phi    = gt.electronPhi[0];
      lepton2D0     = gt.electronD0[0];
      lepton2DZ     = gt.electronDZ[0];
      lepton2RelIso = gt.electronCombIso[0]/gt.electronPt[0];
    }
    if(ao.debug) printf("  Passed lepton kinematics\n");

    // Jet multiplicity
    bool isBoostedCategory=false;
    float deltaPhiZHFJ=-1;
    if(ao.useBoostedCategory) { 
      //bLoad(b["fjMSD_corr"],ientry);
      bLoad(b["fjMSD"],ientry); // TEMPORARY DGH
      bLoad(b["fjPt"],ientry);
      bLoad(b["fjEta"],ientry);
      bLoad(b["fjPhi"],ientry);
      if(gt.fjPt[0]>0) // protection against NaN values of fjPhi
        deltaPhiZHFJ = fabs(TVector2::Phi_mpi_pi(gt.ZBosonPhi-gt.fjPhi));
      if(
        gt.fjPt[0] >= 250 && 
        //gt.fjMSD_corr[0] >= 40 &&
        gt.fjMSD[0] >= 40 && // TEMPORARY DGH
        fabs(gt.fjEta) < 2.4 &&
        gt.ZBosonPt >= 250 &&
        deltaPhiZHFJ >= 2.5
      ) isBoostedCategory=true;
      // If we consider a boosted category splitting,
      // only put boosted (resolved) events in boosted (resolved) regions
      if(isBoostedCategory ^ (ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel))
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
      bLoad(b["nJet"],ientry);
      bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
      bLoad(b["hbbpt"],ientry);
      bLoad(b["hbbm_reg"],ientry);
      if(
        gt.nJet[0]<2 || 
        gt.hbbpt[0]<50 || 
        gt.hbbm_reg[0]<0 ||
        gt.hbbm_reg[0]>250) 
        continue;
    }
    if(ao.debug) printf("  Passed jet kinematics\n");
    if(ao.debug) printf("Passed preselection!\n");
    
    // Met
    bLoad(b["pfmet"],ientry);
    //bLoad(b["pfmetphi"],ientry);
    
    // Jet B-tagging
    bool bjet1IsLoose=false, bjet2IsLoose=false;
    float bjet1btag=-2, bjet2btag=-2;
    bLoad(b["jotCMVA"],ientry);
    bLoad(b["jotCSV"],ientry);
    if(ao.year==2016) {
      bjet1btag = gt.jotCMVA[gt.hbbjtidx[0][0]];
      bjet2btag = gt.jotCMVA[gt.hbbjtidx[0][1]];
      bjet1IsLoose = bjet1btag > cmvaLoose;
      bjet2IsLoose = bjet2btag > cmvaLoose;
    } else if(ao.year==2017) {
      bjet1btag = TMath::Max(gt.jotCSV[gt.hbbjtidx[0][0]],gt.jotCSV[gt.hbbjtidx[0][1]]);
      bjet2btag = TMath::Min(gt.jotCSV[gt.hbbjtidx[0][0]],gt.jotCSV[gt.hbbjtidx[0][1]]);
      bjet1IsLoose = bjet1btag > deepcsvLoose;
      bjet2IsLoose = bjet2btag > deepcsvLoose;
    }

    for(unsigned iJES=0; iJES<NJES; iJES++) {
      bLoad(b[Form("jotPt_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
      bLoad(b[Form("nJot_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
    }
    bLoad(b["nJot"],ientry);
    bLoad(b["jotPt"],ientry);
    bLoad(b["jotEta"],ientry);
    bLoad(b["jotPhi"],ientry);
    bLoad(b["jotFlav"],ientry);
    float bjet1Pt = gt.jotPt[0][gt.hbbjtidx[0][0]];
    float bjet2Pt = gt.jotPt[0][gt.hbbjtidx[0][1]];
    
    if(isBoostedCategory) {
      // Isojets for boosted category
      bLoad(b["fjEta"],ientry);
      bLoad(b["fjPhi"],ientry);
      for(unsigned iJES=0; iJES<NJES; iJES++) 
      for(unsigned char iJ=0; iJ<gt.nJot[iJES]; iJ++) {
        if(fabs(gt.jotEta[iJ])>2.4) continue;
        float dR2JetFatjet=pow(gt.jotEta[iJ]-gt.fjEta,2)+pow(TVector2::Phi_mpi_pi(gt.jotPhi[iJ]-gt.fjPhi),2);
        if(dR2JetFatjet<0.64) continue;
        
        float isojetBtag = (ao.year==2016)? gt.jotCMVA[iJ] : gt.jotCSV[iJ];
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
        if(ao.debug>=3) printf("jet with (pt,|eta|)=(%.2f,%.3f) => (iPt,iEta) = (%d,%d)\n",gt.jotPt[0][iJ],jetAbsEta,iPt,iEta);
        if(iPt>=0 && iEta>=0) {
          float btag = (ao.year==2016)? gt.jotCMVA[iJ] : gt.jotCSV[iJ];
          jetPts    [iPt][iEta].push_back(gt.jotPt[0][iJ]);
          jetEtas   [iPt][iEta].push_back(gt.jotEta[iJ]);
          jetFlavors[iPt][iEta].push_back(gt.jotFlav[iJ]);
          jetBtags  [iPt][iEta].push_back(btag);
        }
      }
    }
    
    // Load branches and calculate stuff for the cuts
    bLoad(b["eventNumber"],ientry);
    bLoad(b["ZBosonPt"],ientry);
    bLoad(b["ZBosonEta"],ientry);
    bLoad(b["ZBosonPhi"],ientry);
    if(isBoostedCategory) {
      bLoad(b["fjPt"],ientry);
      bLoad(b["fjPhi"],ientry);
      bLoad(b["fjEta"],ientry);
      //bLoad(b["fjMSD_corr"],ientry);
      bLoad(b["fjMSD"],ientry); // TEMPORARY DGH
      bLoad(b["fjDoubleCSV"],ientry);
      bLoad(b["fjTau21SD"],ientry);
    } else {
      bLoad(b["hbbpt_reg"],ientry);
      bLoad(b["hbbphi"],ientry);
      bLoad(b["hbbm_reg"],ientry);
    }
    float deltaPhiZH    = -1;
    float ptBalanceZH   = -1;
    float ptBalanceZHFJ = -1;
    float dEtaBjets     = -1;
    float dPhiBjets     = -1;
    float dRBjets       = -1;
    float dEtaZH        = -1;
    float dRZH          = -1;
    float ptBalanceL1L2 = -1;
    float dRL1L2        = -1;
    float mTZHFJ        = -1;
    float dEtaZHFJ      = -1;
    float deltaM        = -1;
    TLorentzVector ZBosonP4, fjP4, ZHFJP4;
    if(isBoostedCategory) {
      ZBosonP4.SetPtEtaPhiM(gt.ZBosonPt,gt.ZBosonEta,gt.ZBosonPhi,gt.ZBosonM);
      fjP4.SetPtEtaPhiM(gt.fjPt[0],gt.fjEta,gt.fjPhi,gt.fjMSD[0]);  // TEMPORARY DGH
      ZHFJP4 = ZBosonP4 + fjP4;
      ptBalanceZHFJ = gt.fjPt[0] / gt.ZBosonPt;
      mTZHFJ = ZHFJP4.Mt();
      ptBalanceL1L2 = lepton1Pt/lepton2Pt;
      dRL1L2 = sqrt(pow(lepton1Eta-lepton2Eta,2)+pow(TVector2::Phi_mpi_pi(lepton1Phi-lepton2Phi),2));
      dEtaZHFJ = fabs(gt.fjEta - gt.ZBosonEta);
      deltaM = fabs(gt.ZBosonM-91.1876);
    } else {
      deltaPhiZH    = fabs(TVector2::Phi_mpi_pi(gt.hbbphi[0] - gt.ZBosonPhi));
      ptBalanceZH   = gt.hbbpt_reg[0] /  gt.ZBosonPt;
      dEtaBjets     = fabs(gt.jotEta[gt.hbbjtidx[0][0]]-gt.jotEta[gt.hbbjtidx[0][1]]);
      dPhiBjets     = fabs(TVector2::Phi_mpi_pi(gt.jotPhi[gt.hbbjtidx[0][0]]-gt.jotPhi[gt.hbbjtidx[0][1]]));
      dRBjets       = sqrt(dEtaBjets*dEtaBjets + dPhiBjets*dPhiBjets);
      dEtaZH        = fabs(gt.ZBosonEta - gt.hbbeta[0]);     
      dRZH          = sqrt(dEtaZH*dEtaZH + deltaPhiZH*deltaPhiZH);
    }
    // deltaPhiZHFJ computed already in Jet multiplicity section
    bool vetoTrainEvts = (
      ao.MVAVarType>1 &&
      ao.selection==kZllHSR && 
      category!=kPlotData &&
      (gt.eventNumber%10)<3
    );

    // Set Selection Bits
    std::map<TString, bool> cut;
    for(unsigned iJES=0; iJES<NJES; iJES++) {
      if(iJES==0) {
        if(isBoostedCategory) {
          cut["ZpTFJ"      ] = gt.ZBosonPt > 250;
          cut["dPhiZHFJ"   ] = fabs(gt.fjPhi - gt.ZBosonPhi) > 2.5;
          cut["btagFJ"     ] = gt.fjDoubleCSV > doubleBCut;
          cut["bvetoFJ"    ] = gt.fjDoubleCSV < doubleBCut;
        } else {
          cut["ZpT"        ] = gt.ZBosonPt > 50;
          cut["btag"       ] = bjet1IsLoose && bjet2IsLoose;
          cut["bveto"      ] = !bjet1IsLoose && !bjet2IsLoose;
        }
        cut["Zmass"      ] = gt.ZBosonM >= 75 && gt.ZBosonM < 105;
        cut["ZmassTight" ] = gt.ZBosonM >= 85 && gt.ZBosonM < 97;
        cut["ZmassSB"    ] = (gt.ZBosonM >= 10 && gt.ZBosonM < 75) || (gt.ZBosonM >= 105 && gt.ZBosonM < 120);
        cut["lowMET"     ] = gt.pfmet[0] < 60;
        cut["boostedCat" ] = isBoostedCategory;
        cut["boostedVeto"] = !isBoostedCategory;
        cut["vetoTrain"  ] = (ao.MVAVarType==1 || (gt.eventNumber%10)>=3 || category==kPlotData);
      } else {
        if(isBoostedCategory) {
          bLoad(b[Form("fjPt_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
          bLoad(b[Form("fjPhi_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
          //bLoad(b[Form("fjMSD_corr_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
          bLoad(b[Form("fjMSD_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry); // TEMPORARY DGH
        } else {
          bLoad(b[Form("hbbpt_reg_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
          bLoad(b[Form("hbbphi_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
          bLoad(b[Form("hbbm_reg_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
        }
      }
      if(iJES==(int)shiftjes::kJESTotalUp ||
         iJES==(int)shiftjes::kJESTotalDown) {
        bLoad(b[Form("pfmet_%s",jesName(static_cast<shiftjes>(iJES)).Data())],ientry);
        cut["lowMET"] = gt.pfmet[iJES] < 60;
      }
      if(isBoostedCategory) {
        //cut["mSD"     ] = gt.fjMSD_corr[iJES] >= 40;
        //cut["mSD_SR"  ] = gt.fjMSD_corr[iJES] >= 80 && gt.fjMSD_corr[iJES]<150;
        //cut["mSD_SB"  ] = cut["mSD"] && gt.fjMSD_corr[iJES]<80;
        // TEMPORARY DGH
        cut["mSD"     ] = gt.fjMSD[iJES] >= 40;
        cut["mSD_SR"  ] = gt.fjMSD[iJES] >= 80 && gt.fjMSD[iJES]<150;
        cut["mSD_SB"  ] = cut["mSD"] && gt.fjMSD[iJES]<80;
        cut["pTFJ"    ] = gt.fjPt[iJES] > 250;
        cut["0ijb"    ] = isojetNBtags[iJES]==0; // not used in SR
        cut["1ijb"    ] = !cut["0ijb"];
      } else {
        cut["dPhiZH"  ] = fabs(gt.hbbphi[iJES] - gt.ZBosonPhi) > 2.5;
        cut["pTjj"    ] = gt.hbbpt_reg[iJES] > 100; // not used
        cut["mjj"     ] = gt.hbbm_reg[iJES] >= 90 && gt.hbbm_reg[iJES] < 150;
        cut["mjjSB"   ] = !cut["mjj"] && gt.hbbm_reg[iJES]<250;
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
        gt.sf_ewkV=ao.ZjetsEWKCorr[0]+ao.ZjetsEWKCorr[1]*
          (TMath::Power((gt.trueGenBosonPt+ao.ZjetsEWKCorr[2]),ao.ZjetsEWKCorr[3]));
        if(isNLOZjets || isLowMassZjets) {
          gt.sf_qcdV=1;
        } else if(useHtBinnedVJetsKFactor) {
          bLoad(b["lheHT"],ientry);
          gt.sf_qcdV=qcdKFactor(type, gt.lheHT);
        } else {
          bLoad(b["lheHT"],ientry);
          bLoad(b["trueGenBosonPt"],ientry);
          gt.sf_qcdV=ao.kfactors_ZJets->GetBinContent(ao.kfactors_ZJets->FindBin(
            TMath::Min(1.39,gt.trueGenBosonPt/1000.),
            TMath::Min(1.39,gt.lheHT         /1000.)
          ));
        }
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
        weight *= 1.0; // gt.sf_muTrig; 
        weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
        weight *= gt.muonSfReco[1] * gt.muonSfTight[1];
        weight_muSF = (1+gt.muonSfUnc[0]) * (1+gt.muonSfUnc[1]);
      } else if(typeLepSel==1) {
        bLoad(b["electronSfReco"],ientry);
        bLoad(b["electronSfMvaWP90"],ientry);
        bLoad(b["electronSfUnc"],ientry);
        bLoad(b["sf_eleTrig"],ientry);
        weight *= 1.0; // gt.sf_eleTrig; 
        weight *= gt.electronSfReco[0] * gt.electronSfMvaWP90[0];
        weight *= gt.electronSfReco[1] * gt.electronSfMvaWP90[1];
        weight_elSF = (1+gt.electronSfUnc[0]) * (1+gt.electronSfUnc[1]);
      } else if (typeLepSel==2) {
        bLoad(b["muonSfReco"],ientry);
        bLoad(b["muonSfTight"],ientry);
        bLoad(b["muonSfUnc"],ientry);
        bLoad(b["electronSfReco"],ientry);
        bLoad(b["electronSfMvaWP90"],ientry);
        bLoad(b["electronSfUnc"],ientry);
        bLoad(b["sf_eleTrig"],ientry);
        weight *= 1.0; // no trigger data/MC SF
        weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
        weight *= gt.electronSfReco[0] * gt.electronSfMvaWP90[0];
        weight_muSF = (1+gt.muonSfUnc[0]);
        weight_elSF = (1+gt.electronSfUnc[0]);
      }

      float recorrect_vhEWK=1, recorrect_vhEWKUp=1, recorrect_vhEWKDown=1;
      if(type==kVZ) {
        bLoad(b["sf_wz"],ientry);
        bLoad(b["sf_zz"],ientry);
        weight *= gt.sf_wz * gt.sf_zz;
      } else if(type==kZH) {
        bLoad(b["sf_vh"],ientry);
        bLoad(b["sf_vhUp"],ientry);
        bLoad(b["sf_vhDown"],ientry);
        //TString vhChannel=lepton1Charge>0?"WplusH":"WminusH";
        recorrect_vhEWK     = vhEWKCorr("ZllH",gt.sf_vh    );
        recorrect_vhEWKUp   = vhEWKCorr("ZllH",gt.sf_vhUp  );
        recorrect_vhEWKDown = vhEWKCorr("ZllH",gt.sf_vhDown);
        weight *= recorrect_vhEWK;
        weight_VHCorrUp = recorrect_vhEWKUp/recorrect_vhEWK;
        weight_VHCorrDown = recorrect_vhEWKDown/recorrect_vhEWK;
      }
      
      //bLoad(b["sf_cmvaWeight_Cent"],ientry);
      //weight *= gt.sf_csvWeights[GeneralTree::csvCent];
      // Reweight the events here if they were not used for training
      
      // Hack for the central Btag weights 
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        if(ao.year==2016) {
          // in 2016 we can use the CSVHelper to calculate the total shift for each jet kinematic bin
          double cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
          weight *= ao.cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetBtags[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCent, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
        } else if(ao.year==2017) {
          // in 2017, we have to calculate the reshape factor for each jet in each kinematic bin
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
      }
      if((ao.selection==kZllHSR || ao.selection==kZllHFJSR) && ao.MVAVarType>1)
        weight *= sf_training;
    
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        if(ao.year==2016) {
          // in 2016 we can use the CSVHelper to calculate the total shift for each jet kinematic bin
          double cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
          double centralWeight = ao.cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetBtags[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCent, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
          for(unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
            GeneralTree::csvShift theShift = gt.csvShifts[iShift];
            weight_btag[iShift][iPt][iEta] = ao.cmvaReweighter->getCSVWeight(
              jetPts[iPt][iEta], jetEtas[iPt][iEta], jetBtags[iPt][iEta], jetFlavors[iPt][iEta],
              theShift,
              cmvaWgtHF, cmvaWgtLF, cmvaWgtCF
            )/centralWeight; 
            if(ao.debug>=3) printf("iPt=%d, iEta=%d, %zu jets, iShift=%d, weight_btag[iShift][iPt][iEta] = %.3f\n", iPt,iEta,jetPts[iPt][iEta].size(),iShift, weight_btag[iShift][iPt][iEta]);
          }
        } else if(ao.year==2017) {
          // in 2017, we have to calculate the reshape factor for each jet in each kinematic bin
          for(unsigned iShift=1; iShift<GeneralTree::nCsvShifts; iShift++) {
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
              if(reshapeFactor>0.001) weight_btag[iShift][iPt][iEta] = reshapeFactor/weight_btag[0][iPt][iEta];
              else                    weight_btag[iShift][iPt][iEta] = 1;
            }
          }
        }
      }
      if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) { // Boosted only weighting
      bLoad(b["fjHighestPtGen"],ientry);
        if(gt.fjHighestPtGen==21 && (
          category==kPlotWbb || category==kPlotWb || category==kPlotWLF ||
          category==kPlotZbb || category==kPlotZb || category==kPlotZLF)
        ) { weight_VGluUp *= 1.2; weight_VGluDown *= 0.8; }
        
        bLoad(b["fjPt"],ientry);
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Boosted_event_topologies
        if(category==kPlotZH || category==kPlotGGZH || category==kPlotZbb) {
          bool passBB = ao.selection==kZllHFJSR || ao.selection==kZllHHeavyFlavorFJCR || ao.selection==kZllHTT2bFJCR;
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
    if(ao.MVAVarType>1) for(unsigned iJES=0; iJES<NJES; iJES++) {
      if((iJES==0 && ao.selection>=kZllHLightFlavorCR && ao.selection<=kZllHPresel) || ao.selection==kZllHSR) {
        ao.mvaInputs[nThread][ 0] = gt.sumEtSoft1                                           ; //"sumEtSoft1"  
        ao.mvaInputs[nThread][ 1] = gt.jotPt[iJES][gt.hbbjtidx[0][0]]                       ; //"bjet1Pt"     
        ao.mvaInputs[nThread][ 2] = gt.jotPt[iJES][gt.hbbjtidx[0][1]]                       ; //"bjet2Pt"     
        ao.mvaInputs[nThread][ 3] = bjet1btag                                               ; //"bjet1btag"   
        ao.mvaInputs[nThread][ 4] = bjet2btag                                               ; //"bjet2btag"   
        ao.mvaInputs[nThread][ 5] = gt.ZBosonPt                                             ; //"ZBosonPt"    
        ao.mvaInputs[nThread][ 6] = gt.ZBosonM                                              ; //"ZBosonM"     
        ao.mvaInputs[nThread][ 7] = gt.ZBosonLep1CosThetaCS                                 ; //"CosThetaCS"  
        ao.mvaInputs[nThread][ 8] = gt.ZBosonLep1CosThetaStar                               ; //"CosThetaStar"
        ao.mvaInputs[nThread][ 9] = gt.hbbpt_reg[iJES]                                      ; //"hbbpt"       
        ao.mvaInputs[nThread][10] = gt.hbbm_reg[iJES]                                       ; //"hbbm"        
        ao.mvaInputs[nThread][11] = fabs(TVector2::Phi_mpi_pi(gt.hbbphi[iJES]-gt.ZBosonPhi)); //"dPhiZH"      
        ao.mvaInputs[nThread][12] = gt.hbbpt_reg[iJES]/gt.ZBosonPt                          ; //"ptBalanceZH" 
        ao.mvaInputs[nThread][13] = dRBjets                                                 ; //"dRBjets"     
        ao.mvaInputs[nThread][14] = dEtaBjets                                               ; //"dEtaBjets"   
        ao.mvaInputs[nThread][15] = gt.nJet[iJES]-2                                         ; //"nAddJet"     
        bdtValue[iJES] = ao.reader[nThread]->EvaluateMVA("BDT");
      } else if((iJES==0 && ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) || ao.selection==kZllHFJSR) { 
        TLorentzVector fjP4_jes, ZHFJP4_jes;
        fjP4_jes.SetPtEtaPhiM(gt.fjPt[iJES],gt.fjEta,gt.fjPhi,gt.fjMSD[iJES]);  // TEMPORARY DGH
        ZHFJP4_jes = ZBosonP4 + fjP4_jes;
        ao.mvaInputs[nThread][ 0] = nIsojet[iJES]                   ; //"nIsojet"       
        ao.mvaInputs[nThread][ 1] = gt.fjPt[iJES]                   ; //"fjPt"          
        ao.mvaInputs[nThread][ 2] = gt.fjMSD[iJES]                  ; //"MSD"           
        ao.mvaInputs[nThread][ 3] = gt.fjTau21SD                    ; //"Tau21SD"       
        ao.mvaInputs[nThread][ 4] = gt.fjPt[iJES]/gt.ZBosonPt;      ; //"ptBalanceZHFJ" 
        ao.mvaInputs[nThread][ 5] = dEtaZHFJ                        ; //"dEtaZHFJ"      
        ao.mvaInputs[nThread][ 6] = deltaPhiZHFJ                    ; //"dPhiZHFJ"      
        ao.mvaInputs[nThread][ 7] = ZHFJP4_jes.Mt()                 ; //"mTZHFJ"        
        ao.mvaInputs[nThread][ 8] = ptBalanceL1L2                   ; //"ptBalanceL1L2" 
        ao.mvaInputs[nThread][ 9] = dRL1L2                          ; //"dRL1L2"        
        ao.mvaInputs[nThread][10] = gt.ZBosonPt                     ; //"ZBosonPt"      
        ao.mvaInputs[nThread][11] = lepton1Pt                       ; //"lepton1Pt"     
        ao.mvaInputs[nThread][12] = lepton2Pt                       ; //"lepton2Pt"     
        ao.mvaInputs[nThread][13] = deltaM                          ; //"deltaM"        
        ao.mvaInputs[nThread][14] = gt.ZBosonLep1CosThetaCS         ; //"CosThetaCS"        
        bdtValue[iJES] = ao.reader[nThread]->EvaluateMVA("BDT");
      }
      switch(ao.MVAVarType) {
        case 1:
        default:
          if(ao.selection==kZllHSR)
            MVAVar[iJES]=gt.hbbpt_reg[iJES];
          else if(ao.selection==kZllHFJSR)
            MVAVar[iJES]=gt.fjPt[iJES];
          else if(ao.selection==kZllHLightFlavorCR ||
            ao.selection==kZllHHeavyFlavorCR       || 
            ao.selection==kZllH2TopCR)
            MVAVar[iJES]=bjet2btag;
          else if(ao.selection==kZllHLightFlavorFJCR ||
            ao.selection==kZllHHeavyFlavorFJCR       ||
            ao.selection==kZllHTT2bFJCR              ||
            ao.selection==kZllHTT1bFJCR)
            //MVAVar[iJES]=gt.fjMSD_corr[iJES];
            MVAVar[iJES]=gt.fjMSD[iJES]; // TEMPORARY DGH
          break;
        case 3:
          if(ao.selection==kZllHSR || ao.selection==kZllHFJSR)
            MVAVar[iJES]=bdtValue[iJES];
          else if(ao.selection==kZllHLightFlavorCR ||
            ao.selection==kZllHHeavyFlavorCR       || 
            ao.selection==kZllH2TopCR)
            MVAVar[iJES]=bjet2btag;
          else if(ao.selection==kZllHLightFlavorFJCR ||
            ao.selection==kZllHHeavyFlavorFJCR       ||
            ao.selection==kZllHTT2bFJCR              ||
            ao.selection==kZllHTT1bFJCR)
            //MVAVar[iJES]=gt.fjMSD_corr[iJES];
            MVAVar[iJES]=gt.fjMSD[iJES]; // TEMPORARY DGH
          break;
      }
    }

    // Fill the uncertainty histograms
    if(category==kData)  {
      bool passFullSel = (selectionBits[0] & ao.selection) != 0;
      if(passFullSel) {
        ao.histo_Baseline    [typeLepSel][category]->Fill(MVAVar[0], weight);
      }
    } else if(category>=kPlotVZbb)  {
      bool passFullSel = (selectionBits[0] & ao.selection) != 0;
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
        if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) { // Boosted only weighting
          ao.histo_VGluUp      [typeLepSel][category]->Fill(MVAVar[0], weight*weight_VGluUp);
          ao.histo_VGluDown    [typeLepSel][category]->Fill(MVAVar[0], weight*weight_VGluDown);
          ao.histo_doubleBUp   [typeLepSel][category]->Fill(MVAVar[0], weight*weight_doubleBUp);
          ao.histo_doubleBDown [typeLepSel][category]->Fill(MVAVar[0], weight*weight_doubleBDown);
        }
      }
      for(unsigned iJES=0; iJES<NJES; iJES++) {
        bool passFullSelJES = (selectionBits[iJES] & ao.selection) != 0;
        if(!passFullSelJES) continue;
        ao.histo_jes[iJES][typeLepSel][category]->Fill(MVAVar[iJES], weight);
      }
    }

    // Fill the plotting histograms and MVA tree (if applicable)
    bool passFullSel = (selectionBits[0] & ao.selection) != 0;
    if(passFullSel && ao.debug>=3) printf("\tPassed this sel\n");
    bLoad(b["ZBosonLep1CosThetaCS"],ientry);
    bLoad(b["ZBosonLep1CosThetaStar"],ientry);
    bLoad(b["ZBosonLep1CosThetaStarFJ"],ientry);
    bLoad(b["sumEtSoft1"],ientry);
    bLoad(b["nSoft2"],ientry);
    bLoad(b["nSoft5"],ientry);
    bLoad(b["nSoft10"],ientry);
    //bLoad(b["fjMSD_corr"],ientry);
    bLoad(b["fjMSD"],ientry); // TEMPORARY DGH
    bLoad(b["fjPt"],ientry);
    bLoad(b["fjDoubleCSV"],ientry);
    // Lock the mutex and fill the MVA tree
    if(ao.selection==kZllHSR && passFullSel && category!=kPlotData) {
      mvaTreeMutex.lock();
      ao.mva_sumEtSoft1   = gt.sumEtSoft1            ; 
      ao.mva_nSoft2       = gt.nSoft2                ; 
      ao.mva_nSoft5       = gt.nSoft5                ; 
      ao.mva_nSoft10      = gt.nSoft10               ; 
      ao.mva_bjet1Pt      = bjet1Pt                  ; 
      ao.mva_bjet2Pt      = bjet2Pt                  ; 
      ao.mva_bjet1btag    = bjet1btag                ; 
      ao.mva_bjet2btag    = bjet2btag                ; 
      ao.mva_ZBosonPt     = gt.ZBosonPt              ; 
      ao.mva_ZBosonM      = gt.ZBosonM               ; 
      ao.mva_CosThetaCS   = gt.ZBosonLep1CosThetaCS  ; 
      ao.mva_CosThetaStar = gt.ZBosonLep1CosThetaStar; 
      ao.mva_hbbpt        = gt.hbbpt_reg[0]          ; 
      ao.mva_hbbm         = gt.hbbm_reg[0]           ; 
      ao.mva_dPhiZH       = deltaPhiZH               ; 
      ao.mva_ptBalanceZH  = ptBalanceZH              ; 
      ao.mva_dRBjets      = dRBjets                  ; 
      ao.mva_dEtaBjets    = dEtaBjets                ; 
      ao.mva_dRZH         = dRZH                     ; 
      ao.mva_dEtaZH       = dEtaZH                   ; 
      ao.mva_nAddJet      = gt.nJet[0]-2             ; 
      ao.mva_weight       = weight                   ; 
      ao.mva_category     = category                 ;
      ao.mva_eventNumber  = gt.eventNumber           ;
      ao.mvaTree->Fill();
      mvaTreeMutex.unlock();
    } else if(ao.selection==kZllHFJSR && passFullSel && category!=kPlotData) {
      mvaTreeMutex.lock();
      ao.mva_nIsojet       = nIsojet[0]      ; 
      //ao.mva_MSD           = fjMSD_corr[0]  ;
      ao.mva_fjPt          = gt.fjPt[0]      ; 
      ao.mva_MSD           = gt.fjMSD[0]     ; // TEMPORARY DGH 
      ao.mva_Tau21SD       = gt.fjTau21SD    ; 
      ao.mva_ptBalanceZHFJ = ptBalanceZHFJ   ; 
      ao.mva_dEtaZHFJ      = dEtaZHFJ        ; 
      ao.mva_dPhiZHFJ      = deltaPhiZHFJ    ; 
      ao.mva_mTZHFJ        = mTZHFJ          ; 
      ao.mva_ptBalanceL1L2 = ptBalanceL1L2   ; 
      ao.mva_dRL1L2        = dRL1L2          ; 
      ao.mva_ZBosonPt      = gt.ZBosonPt     ;
      ao.mva_CosThetaCS    = gt.ZBosonLep1CosThetaCS;
      ao.mva_lepton1Pt     = lepton1Pt       ; 
      ao.mva_lepton2Pt     = lepton2Pt       ; 
      ao.mva_lepton1Eta    = fabs(lepton1Eta);
      ao.mva_lepton2Eta    = fabs(lepton2Eta);
      ao.mva_deltaM        = deltaM          ;
      ao.mva_weight        = weight          ; 
      ao.mva_category      = category        ;
      ao.mva_eventNumber   = gt.eventNumber  ;
      ao.mvaTree->Fill();
      mvaTreeMutex.unlock();
    }
    float theVar;
    for(int p=0; p<nPlots; p++) { 
      bool makePlot=false;
      // Variables -- change the makePlot for n-1 later
      if      (ao.histoNames[p]=="MVAVar"                  ) { theVar = MVAVar[0]                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="lepton1Pt"               ) { theVar = lepton1Pt                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="lepton2Pt"               ) { theVar = lepton2Pt                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="lepton1Eta"              ) { theVar = lepton1Eta                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="lepton2Eta"              ) { theVar = lepton2Eta                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ZBosonPt"                ) { theVar = gt.ZBosonPt                ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ZBosonEta"               ) { theVar = gt.ZBosonEta               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ZBosonPhi"               ) { theVar = gt.ZBosonPhi               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ZBosonM"                 ) { theVar = gt.ZBosonM                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ZBosonLep1CosThetaCS"    ) { theVar = gt.ZBosonLep1CosThetaCS    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ZBosonLep1CosThetaStar"  ) { theVar = gt.ZBosonLep1CosThetaStar  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="Mjj"                     ) { theVar = gt.hbbm_reg[0]             ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="pTjj"                    ) { theVar = gt.hbbpt_reg[0]            ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet1Pt"                 ) { theVar = bjet1Pt                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet2Pt"                 ) { theVar = bjet2Pt                    ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet1btag"               ) { theVar = bjet1btag                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bjet2btag"               ) { theVar = bjet2btag                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="deltaPhiZH"              ) { theVar = deltaPhiZH                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ptBalanceZH"             ) { theVar = ptBalanceZH                ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nJet"                    ) { theVar = gt.nJet[0]                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="sumEtSoft1"              ) { theVar = gt.sumEtSoft1              ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nSoft2"                  ) { theVar = gt.nSoft2                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nSoft5"                  ) { theVar = gt.nSoft5                  ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nSoft10"                 ) { theVar = gt.nSoft10                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="bdtValue"                ) { theVar = bdtValue[0]                ; makePlot = passFullSel; }
      // fatjet
      //else if (ao.histoNames[p]=="mSD"                     ) { theVar = gt.fjMSD_corr[0]           ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="mSD"                     ) { theVar = gt.fjMSD[0]                ; makePlot = passFullSel; } // TEMPORARY DGH
      else if (ao.histoNames[p]=="pTFJ"                    ) { theVar = gt.fjPt[0]                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="Tau21SD"                 ) { theVar = gt.fjTau21SD               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="doubleB"                 ) { theVar = gt.fjDoubleCSV             ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ZBosonLep1CosThetaStarFJ") { theVar = gt.ZBosonLep1CosThetaStarFJ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="deltaEtaZHFJ"            ) { theVar = dEtaZHFJ                   ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="deltaPhiZHFJ"            ) { theVar = deltaPhiZHFJ               ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="mTZHFJ"                  ) { theVar = mTZHFJ                     ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="dRL1L2"                  ) { theVar = dRL1L2                     ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ptBalanceZHFJ"           ) { theVar = ptBalanceZHFJ              ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="ptBalanceL1L2"           ) { theVar = ptBalanceL1L2              ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="nIsojet"                 ) { theVar = nIsojet[0]                 ; makePlot = passFullSel; }
      else if (ao.histoNames[p]=="isojetNBtags"            ) { theVar = isojetNBtags[0]            ; makePlot = passFullSel; }
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
  TString binZptSuffix="";
  if(ao.binZpt>=0 && ao.binZpt<nBinsZpt && 
    !(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel))
    binZptSuffix = Form("_ZptBin%d",ao.binZpt);

  for(unsigned lep=0; lep<nLepSel; lep++) {
    // avoiding unnecessary histograms
    if(leptonStrings[lep]=="em" && !(ao.selection==kZllHSR||ao.selection==kZllHFJSR)) continue; 
    ofstream newcardShape;
    newcardShape.open(Form("MitVHBBAnalysis/datacards/%s/datacard_Z%sH%s%s.txt",
                      dataCardDir.Data(),leptonStrings[lep].Data(),selectionNames[ao.selection].Data(),binZptSuffix.Data()));
    newcardShape << Form("imax * number of channels\n");
    newcardShape << Form("jmax * number of background minus 1\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");

    newcardShape << Form("shapes      *   * datacard_Z%sH%s%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",leptonStrings[lep].Data(),selectionNames[ao.selection].Data(),binZptSuffix.Data());
    newcardShape << Form("shapes data_obs * datacard_Z%sH%s%s.root  histo_%s\n",leptonStrings[lep].Data(),selectionNames[ao.selection].Data(),binZptSuffix.Data(),plotBaseNames[kPlotData].Data());

    newcardShape << Form("Observation %f\n",ao.histo_Baseline[lep][kPlotData]->GetSumOfWeights());

    newcardShape << Form("bin   ");
    for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
    if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0)
      newcardShape << Form("chZ%sH%s  ",leptonStrings[lep].Data(),selectionNames[ao.selection].Data());
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
      if(ic==kPlotZH)
        newcardShape << Form("%d  ",0);
      else if(ic==kPlotGGZH)
        newcardShape << Form("%d  ",-1);
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
      newcardShape << Form("CMS_VH_btag_pt%d_eta%d_%s    shape   ",iPt,iEta,nuisanceName.Data());
      for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0) {
        if((ic==kPlotZbb||ic==kPlotZb||ic==kPlotZLF) && ao.selection==kZllH2TopCR) 
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
      if(ic==kPlotVZbb||ic==kPlotVVLF||ic==kPlotZbb||ic==kPlotZb||ic==kPlotZLF||ic==kPlotZH)
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

    newcardShape<<Form("CMS_VH_TopNorm lnN ");
    for(int ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights()<=0)
        continue;
      newcardShape<< (ic==kPlotTop? "1.15 ":"- ");
    }
    newcardShape<<std::endl;
    
    newcardShape<<Form("CMS_VH_VVNorm lnN ");
    for(int ic=kPlotVZbb; ic!=nPlotCategories; ic++) {
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights()<=0)
        continue;
      newcardShape<< ((ic==kPlotVZbb||ic==kPlotVVLF)? "1.15 ":"- ");
    }
    newcardShape<<std::endl;

    // Normalization and double B scale factors
    if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) {
      // Boosted norms
      newcardShape << Form("SF_ZHFFJ_Zll rateParam * %s 1 [0.2,5]\n",plotBaseNames[kPlotZbb].Data());
      newcardShape << Form("SF_ZHFFJ_Zll rateParam * %s 1 [0.2,5]\n",plotBaseNames[kPlotZb].Data());
      newcardShape << Form("SF_ZLFFJ_Zll rateParam * %s 1 [0.2,5]\n",plotBaseNames[kPlotZLF].Data());
      // In situ measurement of doubleB SF for ZLF, Zb, eff checked in kZllHFJPresel
      newcardShape << "effDoubleB_ZLF  param 0.024 0.0024" << std::endl; 
      newcardShape << "effDoubleB_Zb   param 0.20  0.02  " << std::endl;
      newcardShape << "effSFDoubleB_ZLF extArg 1.0 [0.1,10]" << std::endl;
      newcardShape << "effSFDoubleB_Zb  extArg 1.0 [0.1,10]" << std::endl;
      if(ao.selection==kZllHHeavyFlavorFJCR || ao.selection==kZllHTT2bFJCR || ao.selection==kZllHFJSR) {
        newcardShape << Form("passBB_ZLF rateParam * %s (@0*1.0) effSFDoubleB_ZLF\n" ,plotBaseNames[kPlotZLF].Data());
        newcardShape << Form("passBB_Zb  rateParam * %s (@0*1.0) effSFDoubleB_Zb \n" ,plotBaseNames[kPlotZb].Data());
      } else if(ao.selection==kZllHLightFlavorFJCR || ao.selection==kZllHTT1bFJCR) {
        newcardShape << Form("failBB_ZLF rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_ZLF,effDoubleB_ZLF\n" ,plotBaseNames[kPlotZLF].Data());
        newcardShape << Form("failBB_Zb  rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_Zb,effDoubleB_Zb\n" ,plotBaseNames[kPlotZb].Data());
      }
      // Shape uncertainty using BTV doubleB SF for Zbb, VZ(bb), ZH
      newcardShape << Form("CMS_doubleB shape ");
      for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0) {
        if(ic==kPlotZbb||ic==kPlotVZbb||ic==kPlotZH||ic==kPlotGGZH)
          newcardShape << "1.0 ";
        else
          newcardShape << "- ";
      }
      newcardShape<<Form("VjetsGluFrac shape ");
      for(unsigned ic=kPlotVZbb; ic!=nPlotCategories; ic++)
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() > 0) {
        if(ic==kPlotZbb||ic==kPlotZb||ic==kPlotZLF)
          newcardShape << "1.0 ";
        else
          newcardShape << "- ";
      }
    } else {
      newcardShape << Form("SF_TT%s_Zll  rateParam * %s 1 [0.2,5]\n" ,binZptSuffix.Data(),plotBaseNames[kPlotTT].Data());
      newcardShape << Form("SF_Zbb%s_Zll rateParam * %s 1 [0.2,5]\n" ,binZptSuffix.Data(),plotBaseNames[kPlotZbb].Data());
      newcardShape << Form("SF_Zb%s_Zll  rateParam * %s 1 [0.2,5]\n" ,binZptSuffix.Data(),plotBaseNames[kPlotZb].Data());
      newcardShape << Form("SF_ZLF%s_Zll  rateParam * %s 1 [0.2,5]\n",binZptSuffix.Data(),plotBaseNames[kPlotZLF].Data());
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
  char binZpt=-1, 
  unsigned year=2016
) {
  struct analysisObjects ao;
  ao.selection = selection;
  GeneralTree gt;
  TString binZptSuffix="";
  if(binZpt>=0 && binZpt<nBinsZpt && 
    !(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel))
    binZptSuffix = Form("_ZptBin%d",binZpt);
  // Get shape histograms from file
  for(unsigned lep=0; lep<nLepSel; lep++) {
    // avoiding unnecessary histograms
    if(leptonStrings[lep]=="em" && !(selection==kZllHSR||selection==kZllHFJSR)) continue; 

    char inFileDatacardsName[200];
    sprintf(
      inFileDatacardsName,
      "MitVHBBAnalysis/datacards/%s/datacard_Z%sH%s%s.root",
      dataCardDir.Data(),
      leptonStrings[lep].Data(),
      selectionNames[selection].Data(),
      binZptSuffix.Data()
    );
    TFile* infile = new TFile(inFileDatacardsName,"read");

    for(unsigned ic=kPlotData; ic!=nPlotCategories; ic++) {
      if(ao.histo_Baseline[lep][ic]->GetSumOfWeights() <= 0 && ic!=kPlotData) continue;
      ao.histo_Baseline    [lep][ic] = (TH1F*)infile->Get(Form("histo_%s"                , plotBaseNames[ic].Data()));
      ao.histo_Baseline    [lep][ic]->SetDirectory(0);
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
        ao.histo_btag[iShift][iPt][iEta][lep][ic] = (TH1F*)infile->Get(Form("histo_%s_CMS_VH_btag_pt%d_eta%d_%s",plotBaseNames[ic].Data(),iPt,iEta,btagShiftName(shift)));
        ao.histo_btag[iShift][iPt][iEta][lep][ic]->SetDirectory(0);
      }
      if(ao.selection>=kZllHLightFlavorFJCR && ao.selection<=kZllHFJPresel) {
        ao.histo_VGluUp     [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_VGluUp"         , plotBaseNames[ic].Data()));
        ao.histo_VGluDown   [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_VGluDown"       , plotBaseNames[ic].Data()));
        ao.histo_doubleBUp  [lep][ic] = (TH1F*)infile->Get(Form("histo_%s_doubleBUp"      , plotBaseNames[ic].Data()));
        ao.histo_doubleBDown[lep][ic] = (TH1F*)infile->Get(Form("histo_%s_doubleBDown"    , plotBaseNames[ic].Data()));
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
