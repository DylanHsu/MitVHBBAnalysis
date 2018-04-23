#include "PandaAnalysis/Flat/interface/GeneralTree.h"
#include "PandaAnalysis/Utilities/src/CSVHelper.cc"
#include <Compression.h>
#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <cassert>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unistd.h>
#include "vhbbPlot.h"
#include "formulas.h"

// Script for making light weight plotting/training trees for VH analysis

// useBoostedCategory - This controls whether to cut on the number of fatjets in the resolved regions
// False means you will use all possible events to do the resolved VH
// True means the events with a fatjet are reserved for boosted VH

const bool useHtBinnedVJetsKFactor=true;

using namespace vhbbPlot;
bool vhbbPlotSkim(
  TString inputFileName="", 
  TString outputFileName="",
  vhbbPlot::sampleType sample=nSampleTypes,
  vhbbPlot::selectionType selection = kWHLightFlavorCR,
  bool debug=false,
  Long64_t maxEntries=0,
  bool useBoostedCategory=false,
  int modSplitIndex=-1 // -1: no splitting. 0-9: split into 10 files
) {
  const double theLumi=35900.;
  bool useModSplit = (modSplitIndex>=0 && modSplitIndex<10);

  // Load Shared Objects for ACLIC
  bool loadPandaAnalysisFlat=(0==gSystem->Load("libPandaAnalysisFlat.so"));
  if(!loadPandaAnalysisFlat) { throw std::runtime_error("Error loading shared object libPandaAnalysisFlat.so"); return false; }

  // Handle File IO
  //bool batchMode=false;
  //if(inputFileName!="") batchMode=true;
  //if(batchMode) {
    if(inputFileName==outputFileName) { throw std::runtime_error("Input file and output file are the same file"); return false; }
    TFile *inputFile = TFile::Open(inputFileName, "READ");
    if(!inputFile) { throw std::runtime_error("Problem opening input file"); return false; }
  TTree *events = (TTree*)inputFile->Get("events");
  if(!events) { throw std::runtime_error("Problem loading tree"); return false; }

  TH1F *hDTotalMCWeight=0;
  float sfLHEweights_QCDr1f2;
  float sfLHEweights_QCDr1f5;
  float sfLHEweights_QCDr2f1;
  float sfLHEweights_QCDr2f2;
  float sfLHEweights_QCDr5f1;
  float sfLHEweights_QCDr5f5;
  if(sample!=kData) {
    hDTotalMCWeight = (TH1F*)inputFile->Get("hDTotalMCWeight");
    assert(hDTotalMCWeight); assert(hDTotalMCWeight->GetNbinsX()>=7);
    hDTotalMCWeight->SetDirectory(0);
    float sumLHEweights_nominal = hDTotalMCWeight->GetBinContent(1);
    float sumLHEweights_QCDr1f2 = hDTotalMCWeight->GetBinContent(2);
    float sumLHEweights_QCDr1f5 = hDTotalMCWeight->GetBinContent(3);
    float sumLHEweights_QCDr2f1 = hDTotalMCWeight->GetBinContent(4);
    float sumLHEweights_QCDr2f2 = hDTotalMCWeight->GetBinContent(5);
    float sumLHEweights_QCDr5f1 = hDTotalMCWeight->GetBinContent(6);
    float sumLHEweights_QCDr5f5 = hDTotalMCWeight->GetBinContent(7);
    sfLHEweights_QCDr1f2 = sumLHEweights_nominal / sumLHEweights_QCDr1f2; 
    sfLHEweights_QCDr1f5 = sumLHEweights_nominal / sumLHEweights_QCDr1f5; 
    sfLHEweights_QCDr2f1 = sumLHEweights_nominal / sumLHEweights_QCDr2f1; 
    sfLHEweights_QCDr2f2 = sumLHEweights_nominal / sumLHEweights_QCDr2f2; 
    sfLHEweights_QCDr5f1 = sumLHEweights_nominal / sumLHEweights_QCDr5f1; 
    sfLHEweights_QCDr5f5 = sumLHEweights_nominal / sumLHEweights_QCDr5f5; 
  }
 
  GeneralTree gt;
  gt.monohiggs      = true;
  gt.vbf            = false;
  gt.fatjet         = true;
  gt.leptonic       = true;
  gt.hfCounting     = true;
  gt.btagWeights    = true;
  gt.useCMVA        = true;
  // Branches not in GeneralTree;
  std::map<TString, void*> extraAddresses;
  float normalizedWeight; extraAddresses["normalizedWeight"] = &normalizedWeight;
  // Map of the branches
  std::map<TString, TBranch*> b;
  TTree *dummyTree = new TTree("dummyTree","dummyTree");
  gt.WriteTree(dummyTree);
  TObjArray *listOfBranches = events->GetListOfBranches();
  for(unsigned iB=0; iB<(unsigned)listOfBranches->GetEntries(); iB++) {
    TBranch *branch = (TBranch*)listOfBranches->At(iB);
    TString branchName = branch->GetName();
    //if(branchName.Contains("fj1GenPt") || branchName.Contains("fj1Gen")) continue;
    //if(branchName.Contains("genFj1")) continue;
    bool isExtraBranch=false;
    auto x = extraAddresses.find(branchName);
    isExtraBranch= (x!=extraAddresses.end());
    if(isExtraBranch) {
      b[x->first] = events->GetBranch(x->first);
      if(!b[x->first]) { throw std::runtime_error(Form("Extra branch %s could not be found in events tree\n", x->first.Data())); return false; }
      b[x->first]->SetAddress(x->second);
      if(debug) printf("Booking extra branch \"%s\"\n", x->first.Data());
      continue;
    }
    TBranch *dummyBranch = dummyTree->GetBranch(branchName);
    if(!dummyBranch) { throw std::runtime_error(Form("WARNING: Couldn't find branch \"%s\" in the dummyTree, skipping", branchName.Data())); continue; }
    void* address=dummyBranch->GetAddress();
    b[branchName] = events->GetBranch(branchName);
    if(!b[branchName]) { throw std::runtime_error(Form("Branch \"%s\" could not be found in events tree", x->first.Data())); return false; }
    b[branchName]->SetAddress(address);
    if(debug) printf("Booking GeneralTree branch \"%s\"\n", branchName.Data());
  }
   
  // Load corrections to apply offline
  CSVHelper *cmvaReweighter = new CSVHelper("PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root"   , "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root"   , 5);
  
  TFile *kfactorsFile = TFile::Open("PandaAnalysis/data/kfactors.root","read"); assert(kfactorsFile);
  TH1F *hWjets_relErr_QCDfScaleUp   =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/fact_up"  )->Clone("hWjets_relErr_QCDfScaleUp"  ); hWjets_relErr_QCDfScaleUp  ->SetDirectory(0);
  TH1F *hWjets_relErr_QCDfScaleDown =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/fact_down")->Clone("hWjets_relErr_QCDfScaleDown"); hWjets_relErr_QCDfScaleDown->SetDirectory(0);
  TH1F *hWjets_relErr_QCDrScaleUp   =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/ren_up"   )->Clone("hWjets_relErr_QCDrScaleUp"  ); hWjets_relErr_QCDrScaleUp  ->SetDirectory(0);
  TH1F *hWjets_relErr_QCDrScaleDown =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/ren_down" )->Clone("hWjets_relErr_QCDrScaleDown"); hWjets_relErr_QCDrScaleDown->SetDirectory(0);
  TH1F *hZjets_relErr_QCDfScaleUp   =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/fact_up"  )->Clone("hWjets_relErr_QCDfScaleUp"  ); hWjets_relErr_QCDfScaleUp  ->SetDirectory(0);
  TH1F *hZjets_relErr_QCDfScaleDown =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/fact_down")->Clone("hWjets_relErr_QCDfScaleDown"); hWjets_relErr_QCDfScaleDown->SetDirectory(0);
  TH1F *hZjets_relErr_QCDrScaleUp   =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/ren_up"   )->Clone("hWjets_relErr_QCDrScaleUp"  ); hWjets_relErr_QCDrScaleUp  ->SetDirectory(0);
  TH1F *hZjets_relErr_QCDrScaleDown =(TH1F*)kfactorsFile->Get("WJets_012j_NLO/ren_down" )->Clone("hWjets_relErr_QCDrScaleDown"); hWjets_relErr_QCDrScaleDown->SetDirectory(0);
  hWjets_relErr_QCDfScaleUp   ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  hWjets_relErr_QCDfScaleDown ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  hWjets_relErr_QCDrScaleUp   ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  hWjets_relErr_QCDrScaleDown ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  hZjets_relErr_QCDfScaleUp   ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  hZjets_relErr_QCDfScaleDown ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  hZjets_relErr_QCDrScaleUp   ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  hZjets_relErr_QCDrScaleDown ->Divide((TH1F*)kfactorsFile->Get("WJets_012j_NLO/nominal"));
  // Done loading corrections

  // Define the output tree
  TFile *outputFile = TFile::Open(outputFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  TTree *plotTree = new TTree("plotTree","Tree for making plots");
  unsigned char typeLepSel, theCategory;
  unsigned selectionBits, selectionBits_jesUp, selectionBits_jesDown, nMinusOneBits;
  float lepton1Pt, lepton1Eta, lepton1Phi, lepton1RelIso, lepton1D0, lepton1DZ;
  int lepton1Flav, lepton1Charge;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1Eta, hbbJet1Phi;
  float hbbJet2Pt, hbbJet2Eta, hbbJet2Phi;
  float hbbJet1PtUp, hbbJet1PtDown, hbbJet2PtUp, hbbJet2PtDown;
  float hbbDijetPt, hbbDijetPtUp, hbbDijetPtDown;
  float hbbDijetMass, hbbDijetMassUp, hbbDijetMassDown;
  float bDiscrMin, bDiscrMax;
  float deltaPhiLep1Met;
  float deltaPhiVH, deltaPhiVH_jesUp, deltaPhiVH_jesDown;
  float topWBosonPt, topWBosonPt_jesUp, topWBosonPt_jesDown, topWBosonPhi_jesUp, topWBosonPhi_jesDown;
  float mT_jesUp, mT_jesDown; // need to be in PandaExpress ntuples
  float fj1MSD_corr,fj1MSD_corr_jesUp,fj1MSD_corr_jesDown; 
  
  // For boosted categories, the number of 30 GeV AK4 jets not in the fat jet
  int nIsojet=0, nIsojet_jesUp=0, nIsojet_jesDown=0; 
  vector<unsigned char> isojets, isojets_jesUp, isojets_jesDown;
  unsigned char isojetNBtags, isojetNBtags_jesUp, isojetNBtags_jesDown;   
  
  // CMVA jet kinematic decorrelated weight nuisances
  vector<double> jetPts[5][3], jetPtsUp[5][3], jetPtsDown[5][3], jetEtas[5][3], jetCMVAs[5][3];
  vector<int> jetFlavors[5][3];
  for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
    jetPts    [iPt][iEta].reserve(20);
    jetPtsUp  [iPt][iEta].reserve(20);
    jetPtsDown[iPt][iEta].reserve(20);
    jetEtas   [iPt][iEta].reserve(20);
    jetCMVAs  [iPt][iEta].reserve(20);
    jetFlavors[iPt][iEta].reserve(20);
  }

  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_QCDr1f2, weight_QCDr1f5, weight_QCDr2f1, weight_QCDr2f2, weight_QCDr5f1, weight_QCDr5f5;
  float weight_NLOQCDrUp, weight_NLOQCDrDown, weight_NLOQCDfUp, weight_NLOQCDfDown; // Renorm/factorization scales of the NLO QCD V+jets k-factors *deep breath*
  float weight_lepSFUp;
  float weight_cmvaLFUp      [5][3] , weight_cmvaLFDown      [5][3]; 
  float weight_cmvaHFUp      [5][3] , weight_cmvaHFDown      [5][3]; 
  float weight_cmvaHFStats1Up[5][3] , weight_cmvaHFStats1Down[5][3]; 
  float weight_cmvaHFStats2Up[5][3] , weight_cmvaHFStats2Down[5][3]; 
  float weight_cmvaLFStats1Up[5][3] , weight_cmvaLFStats1Down[5][3]; 
  float weight_cmvaLFStats2Up[5][3] , weight_cmvaLFStats2Down[5][3]; 
  float weight_cmvaCErr1Up   [5][3] , weight_cmvaCErr1Down   [5][3]; 
  float weight_cmvaCErr2Up   [5][3] , weight_cmvaCErr2Down   [5][3]; 
  float weight_cmvaJESUp     [5][3] , weight_cmvaJESDown     [5][3]; 
  float weight_VHCorrUp, weight_VHCorrDown; 

  float weight_btagBUp, weight_btagBDown;
  float weight_btagMUp, weight_btagMDown;
  float *sf_sjbtag0      = (float*)dummyTree->GetBranch("sf_sjbtag0"     )->GetAddress();
  float *sf_sjbtag0BUp   = (float*)dummyTree->GetBranch("sf_sjbtag0BUp"  )->GetAddress();
  float *sf_sjbtag0BDown = (float*)dummyTree->GetBranch("sf_sjbtag0BDown")->GetAddress();
  float *sf_sjbtag0MUp   = (float*)dummyTree->GetBranch("sf_sjbtag0MUp"  )->GetAddress();
  float *sf_sjbtag0MDown = (float*)dummyTree->GetBranch("sf_sjbtag0MDown")->GetAddress();
  float *sf_sjbtagGT0      = (float*)dummyTree->GetBranch("sf_sjbtagGT0"     )->GetAddress();
  float *sf_sjbtagGT0BUp   = (float*)dummyTree->GetBranch("sf_sjbtagGT0BUp"  )->GetAddress();
  float *sf_sjbtagGT0BDown = (float*)dummyTree->GetBranch("sf_sjbtagGT0BDown")->GetAddress();
  float *sf_sjbtagGT0MUp   = (float*)dummyTree->GetBranch("sf_sjbtagGT0MUp"  )->GetAddress();
  float *sf_sjbtagGT0MDown = (float*)dummyTree->GetBranch("sf_sjbtagGT0MDown")->GetAddress();
  //float *sf_sjbtag2      = (float*)dummyTree->GetBranch("sf_sjbtag2"     )->GetAddress();
  //float *sf_sjbtag2BUp   = (float*)dummyTree->GetBranch("sf_sjbtag2BUp"  )->GetAddress();
  //float *sf_sjbtag2BDown = (float*)dummyTree->GetBranch("sf_sjbtag2BDown")->GetAddress();
  //float *sf_sjbtag2MUp   = (float*)dummyTree->GetBranch("sf_sjbtag2MUp"  )->GetAddress();
  //float *sf_sjbtag2MDown = (float*)dummyTree->GetBranch("sf_sjbtag2MDown")->GetAddress();
  
  plotTree->Branch("runNumber"       , &gt.runNumber    );
  plotTree->Branch("lumiNumber"      , &gt.lumiNumber   );
  plotTree->Branch("eventNumber"     , &gt.eventNumber  );
  plotTree->Branch("selectionBits"   , &selectionBits   );
  plotTree->Branch("selectionBits_jesUp"   , &selectionBits_jesUp   );
  plotTree->Branch("selectionBits_jesDown" , &selectionBits_jesDown );
  plotTree->Branch("nMinusOneBits"   , &nMinusOneBits   );
  plotTree->Branch("theCategory"     , &theCategory     );
  plotTree->Branch("pfmet"           , &gt.pfmet        );
  plotTree->Branch("pfmetUp"         , &gt.pfmetUp      );
  plotTree->Branch("pfmetDown"       , &gt.pfmetDown    );
  plotTree->Branch("pfmetsig"        , &gt.pfmetsig     );
  plotTree->Branch("puppimetsig"     , &gt.puppimetsig  );
  plotTree->Branch("pfmetphi"        , &gt.pfmetphi     );
  plotTree->Branch("topWBosonPt"             , &gt.topWBosonPt           );
  plotTree->Branch("topWBosonPt_jesUp"       , &topWBosonPt_jesUp        );
  plotTree->Branch("topWBosonPt_jesDown"     , &topWBosonPt_jesDown      );
  plotTree->Branch("topWBosonPhi"            , &gt.topWBosonPhi          );
  plotTree->Branch("topWBosonPhi_jesUp"      , &topWBosonPhi_jesUp       );
  plotTree->Branch("topWBosonPhi_jesDown"    , &topWBosonPhi_jesDown     );
  plotTree->Branch("mT"                      , &gt.mT                    );
  plotTree->Branch("mT_jesUp"                , &mT_jesUp                 );
  plotTree->Branch("mT_jesDown"              , &mT_jesDown               );
  plotTree->Branch("typeLepSel"              , &typeLepSel               );
  plotTree->Branch("lepton1Pt"               , &lepton1Pt                );
  plotTree->Branch("lepton1Eta"              , &lepton1Eta               );
  plotTree->Branch("lepton1Phi"              , &lepton1Phi               );
  plotTree->Branch("lepton1RelIso"           , &lepton1RelIso            );
  plotTree->Branch("lepton1D0"               , &lepton1D0                );
  plotTree->Branch("lepton1DZ"               , &lepton1DZ                );
  plotTree->Branch("lepton1Flav"             , &lepton1Flav              );
  plotTree->Branch("lepton1Charge"           , &lepton1Charge            );
  plotTree->Branch("deltaPhiLep1Met"         , &deltaPhiLep1Met          );
  plotTree->Branch("deltaPhiVH"              , &deltaPhiVH               );
  plotTree->Branch("deltaPhiVH_jesUp"        , &deltaPhiVH_jesUp         );
  plotTree->Branch("deltaPhiVH_jesDown"      , &deltaPhiVH_jesDown       );
  plotTree->Branch("nFatjet"                 , &gt.nFatjet               );

  plotTree->Branch("weight"                  , &weight                   );
  plotTree->Branch("weight_VHCorrUp"         , &weight_VHCorrUp          );
  plotTree->Branch("weight_VHCorrDown"       , &weight_VHCorrDown        );
  plotTree->Branch("weight_pdfUp"            , &weight_pdfUp             );
  plotTree->Branch("weight_pdfDown"          , &weight_pdfDown           );
  plotTree->Branch("weight_QCDr1f2"          , &weight_QCDr1f2           );
  plotTree->Branch("weight_QCDr1f5"          , &weight_QCDr1f5           );
  plotTree->Branch("weight_QCDr2f1"          , &weight_QCDr2f1           );
  plotTree->Branch("weight_QCDr2f2"          , &weight_QCDr2f2           );
  plotTree->Branch("weight_QCDr5f1"          , &weight_QCDr5f1           );
  plotTree->Branch("weight_QCDr5f5"          , &weight_QCDr5f5           );
  plotTree->Branch("weight_NLOQCDrUp"        , &weight_NLOQCDrUp         );
  plotTree->Branch("weight_NLOQCDrDown"      , &weight_NLOQCDrDown       );
  plotTree->Branch("weight_NLOQCDfUp"        , &weight_NLOQCDfUp         );
  plotTree->Branch("weight_NLOQCDfDown"      , &weight_NLOQCDfDown       );
  plotTree->Branch("weight_lepSFUp"          , &weight_lepSFUp           );
  plotTree->Branch("weight_btagBUp"          , &weight_btagBUp           );
  plotTree->Branch("weight_btagBDown"        , &weight_btagBDown         );
  plotTree->Branch("weight_btagMUp"          , &weight_btagMUp           );
  plotTree->Branch("weight_btagMDown"        , &weight_btagMDown         );
  plotTree->Branch("weight_cmvaLFUp"         , weight_cmvaLFUp         , "weight_cmvaLFUp[5][3]/F"        );
  plotTree->Branch("weight_cmvaHFUp"         , weight_cmvaHFUp         , "weight_cmvaHFUp[5][3]/F"        );
  plotTree->Branch("weight_cmvaHFStats1Up"   , weight_cmvaHFStats1Up   , "weight_cmvaHFStats1Up[5][3]/F"  );
  plotTree->Branch("weight_cmvaHFStats2Up"   , weight_cmvaHFStats2Up   , "weight_cmvaHFStats2Up[5][3]/F"  );
  plotTree->Branch("weight_cmvaLFStats1Up"   , weight_cmvaLFStats1Up   , "weight_cmvaLFStats1Up[5][3]/F"  );
  plotTree->Branch("weight_cmvaLFStats2Up"   , weight_cmvaLFStats2Up   , "weight_cmvaLFStats2Up[5][3]/F"  );
  plotTree->Branch("weight_cmvaCErr1Up"      , weight_cmvaCErr1Up      , "weight_cmvaCErr1Up[5][3]/F"     );
  plotTree->Branch("weight_cmvaCErr2Up"      , weight_cmvaCErr2Up      , "weight_cmvaCErr2Up[5][3]/F"     );
  plotTree->Branch("weight_cmvaJESUp"        , weight_cmvaJESUp        , "weight_cmvaJESUp[5][3]/F"       );
  plotTree->Branch("weight_cmvaLFDown"       , weight_cmvaLFDown       , "weight_cmvaLFDown[5][3]/F"      );
  plotTree->Branch("weight_cmvaHFDown"       , weight_cmvaHFDown       , "weight_cmvaHFDown[5][3]/F"      );
  plotTree->Branch("weight_cmvaHFStats1Down" , weight_cmvaHFStats1Down , "weight_cmvaHFStats1Down[5][3]/F");
  plotTree->Branch("weight_cmvaHFStats2Down" , weight_cmvaHFStats2Down , "weight_cmvaHFStats2Down[5][3]/F");
  plotTree->Branch("weight_cmvaLFStats1Down" , weight_cmvaLFStats1Down , "weight_cmvaLFStats1Down[5][3]/F");
  plotTree->Branch("weight_cmvaLFStats2Down" , weight_cmvaLFStats2Down , "weight_cmvaLFStats2Down[5][3]/F");
  plotTree->Branch("weight_cmvaCErr1Down"    , weight_cmvaCErr1Down    , "weight_cmvaCErr1Down[5][3]/F"   );
  plotTree->Branch("weight_cmvaCErr2Down"    , weight_cmvaCErr2Down    , "weight_cmvaCErr2Down[5][3]/F"   );
  plotTree->Branch("weight_cmvaJESDown"      , weight_cmvaJESDown      , "weight_cmvaJESDown[5][3]/F"     );
  if(selection>=kWHLightFlavorCR && selection<=kWHPresel) {
    plotTree->Branch("nJet"                    , &gt.nJet                  );
    plotTree->Branch("nJet_jesUp"              , &gt.nJet_jesUp            );
    plotTree->Branch("nJet_jesDown"            , &gt.nJet_jesDown          );
    plotTree->Branch("hbbDijetPt"              , &hbbDijetPt               );
    plotTree->Branch("hbbDijetPtUp"            , &hbbDijetPtUp             );
    plotTree->Branch("hbbDijetPtDown"          , &hbbDijetPtDown           );
    plotTree->Branch("hbbDijetPt"              , &hbbDijetPt               );
    plotTree->Branch("hbbDijetMass"            , &hbbDijetMass             );
    plotTree->Branch("hbbDijetMassUp"          , &hbbDijetMassUp           );
    plotTree->Branch("hbbDijetMassDown"        , &hbbDijetMassDown         );
    plotTree->Branch("bDiscrMin"               , &bDiscrMin                );
    plotTree->Branch("bDiscrMax"               , &bDiscrMax                );
    plotTree->Branch("hbbJet1Pt"               , &hbbJet1Pt                );
    plotTree->Branch("hbbJet1Eta"              , &hbbJet1Eta               );
    plotTree->Branch("hbbJet1Phi"              , &hbbJet1Phi               );
    plotTree->Branch("hbbJet2Pt"               , &hbbJet2Pt                );
    plotTree->Branch("hbbJet2Eta"              , &hbbJet2Eta               );
    plotTree->Branch("hbbJet2Phi"              , &hbbJet2Phi               );
    plotTree->Branch("hbbJet1PtUp"             , &hbbJet1PtUp              );
    plotTree->Branch("hbbJet1PtDown"           , &hbbJet1PtDown            );
    plotTree->Branch("hbbJet2PtUp"             , &hbbJet2PtUp              );
    plotTree->Branch("hbbJet2PtDown"           , &hbbJet2PtDown            );
    plotTree->Branch("topMassLep1Met"          , &gt.topMassLep1Met        );
    plotTree->Branch("topMassLep1Met_jesUp"    , &gt.topMassLep1Met_jesUp  );
    plotTree->Branch("topMassLep1Met_jesDown"  , &gt.topMassLep1Met_jesDown);
    plotTree->Branch("topWBosonEta"            , &gt.topWBosonEta          );
    plotTree->Branch("topWBosonCosThetaCS"     , &gt.topWBosonCosThetaCS   );
    plotTree->Branch("sumEtSoft1"              , &gt.sumEtSoft1            );
    plotTree->Branch("nSoft2"                  , &gt.nSoft2                );
    plotTree->Branch("nSoft5"                  , &gt.nSoft5                );
    plotTree->Branch("nSoft10"                 , &gt.nSoft10               );
    plotTree->Branch("hbbCosThetaJJ"           , &gt.hbbCosThetaJJ         );
    plotTree->Branch("hbbCosThetaCSJ1"         , &gt.hbbCosThetaCSJ1       );
  }
  // Save fatjet variables if we are in a boosted category, or if we aren't using a boosted category
  if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
    plotTree->Branch("nIsojet"               , &nIsojet                  );
    plotTree->Branch("isojetNBtags"          , &isojetNBtags             );
    plotTree->Branch("fj1Tau32"              , &gt.fj1Tau32              );
    plotTree->Branch("fj1Tau21"              , &gt.fj1Tau21              );
    plotTree->Branch("fj1Tau32SD"            , &gt.fj1Tau32SD            ); 
    plotTree->Branch("fj1Tau21SD"            , &gt.fj1Tau21SD            ); 
    plotTree->Branch("fj1MSD"                , &gt.fj1MSD                );
    plotTree->Branch("fj1MSDScaleUp"         , &gt.fj1MSDScaleUp         );    
    plotTree->Branch("fj1MSDScaleDown"       , &gt.fj1MSDScaleDown       );      
    plotTree->Branch("fj1MSD_corr"           , &fj1MSD_corr              );      
    plotTree->Branch("fj1MSD_corr_jesUp"     , &fj1MSD_corr_jesUp        );
    plotTree->Branch("fj1MSD_corr_jesDown"   , &fj1MSD_corr_jesDown      );    
    plotTree->Branch("fj1MSDSmeared"         , &gt.fj1MSDSmeared         );    
    plotTree->Branch("fj1MSDSmearedUp"       , &gt.fj1MSDSmearedUp       );      
    plotTree->Branch("fj1MSDSmearedDown"     , &gt.fj1MSDSmearedDown     );        
    plotTree->Branch("fj1Pt"                 , &gt.fj1Pt                 );
    plotTree->Branch("fj1PtScaleUp"          , &gt.fj1PtScaleUp          );   
    plotTree->Branch("fj1PtScaleDown"        , &gt.fj1PtScaleDown        );     
    plotTree->Branch("fj1PtSmeared"          , &gt.fj1PtSmeared          );   
    plotTree->Branch("fj1PtSmearedUp"        , &gt.fj1PtSmearedUp        );     
    plotTree->Branch("fj1PtSmearedDown"      , &gt.fj1PtSmearedDown      );       
    plotTree->Branch("fj1HighestPtGen"       , &gt.fj1HighestPtGen       );       
    plotTree->Branch("fj1Phi"                , &gt.fj1Phi                );
    plotTree->Branch("fj1Eta"                , &gt.fj1Eta                );
    plotTree->Branch("fj1M"                  , &gt.fj1M                  );
    plotTree->Branch("fj1MaxCSV"             , &gt.fj1MaxCSV             );
    plotTree->Branch("fj1MinCSV"             , &gt.fj1MinCSV             );
    plotTree->Branch("fj1DoubleCSV"          , &gt.fj1DoubleCSV          );   
    plotTree->Branch("fj1HTTMass"            , &gt.fj1HTTMass            ); 
    plotTree->Branch("fj1HTTFRec"            , &gt.fj1HTTFRec            ); 
    plotTree->Branch("fj1SDEFrac100"         , &gt.fj1SDEFrac100         );    
    plotTree->Branch("fj1ECFN_1_1_05", (float*)dummyTree->GetBranch("fj1ECFN_1_1_05")->GetAddress());
    plotTree->Branch("fj1ECFN_2_1_05", (float*)dummyTree->GetBranch("fj1ECFN_2_1_05")->GetAddress());
    plotTree->Branch("fj1ECFN_3_1_05", (float*)dummyTree->GetBranch("fj1ECFN_3_1_05")->GetAddress());
    plotTree->Branch("fj1ECFN_1_2_05", (float*)dummyTree->GetBranch("fj1ECFN_1_2_05")->GetAddress());
    plotTree->Branch("fj1ECFN_2_2_05", (float*)dummyTree->GetBranch("fj1ECFN_2_2_05")->GetAddress());
    plotTree->Branch("fj1ECFN_3_2_05", (float*)dummyTree->GetBranch("fj1ECFN_3_2_05")->GetAddress());
    plotTree->Branch("fj1ECFN_1_3_05", (float*)dummyTree->GetBranch("fj1ECFN_1_3_05")->GetAddress());
    plotTree->Branch("fj1ECFN_2_3_05", (float*)dummyTree->GetBranch("fj1ECFN_2_3_05")->GetAddress());
    plotTree->Branch("fj1ECFN_3_3_05", (float*)dummyTree->GetBranch("fj1ECFN_3_3_05")->GetAddress());
    plotTree->Branch("fj1ECFN_1_4_05", (float*)dummyTree->GetBranch("fj1ECFN_1_4_05")->GetAddress());
    plotTree->Branch("fj1ECFN_2_4_05", (float*)dummyTree->GetBranch("fj1ECFN_2_4_05")->GetAddress());
    plotTree->Branch("fj1ECFN_3_4_05", (float*)dummyTree->GetBranch("fj1ECFN_3_4_05")->GetAddress());
    plotTree->Branch("fj1ECFN_1_1_10", (float*)dummyTree->GetBranch("fj1ECFN_1_1_10")->GetAddress());
    plotTree->Branch("fj1ECFN_2_1_10", (float*)dummyTree->GetBranch("fj1ECFN_2_1_10")->GetAddress());
    plotTree->Branch("fj1ECFN_3_1_10", (float*)dummyTree->GetBranch("fj1ECFN_3_1_10")->GetAddress());
    plotTree->Branch("fj1ECFN_1_2_10", (float*)dummyTree->GetBranch("fj1ECFN_1_2_10")->GetAddress());
    plotTree->Branch("fj1ECFN_2_2_10", (float*)dummyTree->GetBranch("fj1ECFN_2_2_10")->GetAddress());
    plotTree->Branch("fj1ECFN_3_2_10", (float*)dummyTree->GetBranch("fj1ECFN_3_2_10")->GetAddress());
    plotTree->Branch("fj1ECFN_1_3_10", (float*)dummyTree->GetBranch("fj1ECFN_1_3_10")->GetAddress());
    plotTree->Branch("fj1ECFN_2_3_10", (float*)dummyTree->GetBranch("fj1ECFN_2_3_10")->GetAddress());
    plotTree->Branch("fj1ECFN_3_3_10", (float*)dummyTree->GetBranch("fj1ECFN_3_3_10")->GetAddress());
    plotTree->Branch("fj1ECFN_1_4_10", (float*)dummyTree->GetBranch("fj1ECFN_1_4_10")->GetAddress());
    plotTree->Branch("fj1ECFN_2_4_10", (float*)dummyTree->GetBranch("fj1ECFN_2_4_10")->GetAddress());
    plotTree->Branch("fj1ECFN_3_4_10", (float*)dummyTree->GetBranch("fj1ECFN_3_4_10")->GetAddress());
    plotTree->Branch("fj1ECFN_1_1_20", (float*)dummyTree->GetBranch("fj1ECFN_1_1_20")->GetAddress());
    plotTree->Branch("fj1ECFN_2_1_20", (float*)dummyTree->GetBranch("fj1ECFN_2_1_20")->GetAddress());
    plotTree->Branch("fj1ECFN_3_1_20", (float*)dummyTree->GetBranch("fj1ECFN_3_1_20")->GetAddress());
    plotTree->Branch("fj1ECFN_1_2_20", (float*)dummyTree->GetBranch("fj1ECFN_1_2_20")->GetAddress());
    plotTree->Branch("fj1ECFN_2_2_20", (float*)dummyTree->GetBranch("fj1ECFN_2_2_20")->GetAddress());
    plotTree->Branch("fj1ECFN_3_2_20", (float*)dummyTree->GetBranch("fj1ECFN_3_2_20")->GetAddress());
    plotTree->Branch("fj1ECFN_1_3_20", (float*)dummyTree->GetBranch("fj1ECFN_1_3_20")->GetAddress());
    plotTree->Branch("fj1ECFN_2_3_20", (float*)dummyTree->GetBranch("fj1ECFN_2_3_20")->GetAddress());
    plotTree->Branch("fj1ECFN_3_3_20", (float*)dummyTree->GetBranch("fj1ECFN_3_3_20")->GetAddress());
    plotTree->Branch("fj1ECFN_1_4_20", (float*)dummyTree->GetBranch("fj1ECFN_1_4_20")->GetAddress());
    plotTree->Branch("fj1ECFN_2_4_20", (float*)dummyTree->GetBranch("fj1ECFN_2_4_20")->GetAddress());
    plotTree->Branch("fj1ECFN_3_4_20", (float*)dummyTree->GetBranch("fj1ECFN_3_4_20")->GetAddress());
    plotTree->Branch("fj1ECFN_1_1_40", (float*)dummyTree->GetBranch("fj1ECFN_1_1_40")->GetAddress());
    plotTree->Branch("fj1ECFN_2_1_40", (float*)dummyTree->GetBranch("fj1ECFN_2_1_40")->GetAddress());
    plotTree->Branch("fj1ECFN_3_1_40", (float*)dummyTree->GetBranch("fj1ECFN_3_1_40")->GetAddress());
    plotTree->Branch("fj1ECFN_1_2_40", (float*)dummyTree->GetBranch("fj1ECFN_1_2_40")->GetAddress());
    plotTree->Branch("fj1ECFN_2_2_40", (float*)dummyTree->GetBranch("fj1ECFN_2_2_40")->GetAddress());
    plotTree->Branch("fj1ECFN_3_2_40", (float*)dummyTree->GetBranch("fj1ECFN_3_2_40")->GetAddress());
    plotTree->Branch("fj1ECFN_1_3_40", (float*)dummyTree->GetBranch("fj1ECFN_1_3_40")->GetAddress());
    plotTree->Branch("fj1ECFN_2_3_40", (float*)dummyTree->GetBranch("fj1ECFN_2_3_40")->GetAddress());
    plotTree->Branch("fj1ECFN_3_3_40", (float*)dummyTree->GetBranch("fj1ECFN_3_3_40")->GetAddress());
    plotTree->Branch("fj1ECFN_1_4_40", (float*)dummyTree->GetBranch("fj1ECFN_1_4_40")->GetAddress());
    plotTree->Branch("fj1ECFN_2_4_40", (float*)dummyTree->GetBranch("fj1ECFN_2_4_40")->GetAddress());
    plotTree->Branch("fj1ECFN_3_4_40", (float*)dummyTree->GetBranch("fj1ECFN_3_4_40")->GetAddress());
  }
  
  vector<double> theEWKCorrPars; 
  if(useHtBinnedVJetsKFactor) {
    unsigned htLow=0, htHigh=0;
    if(sample==kWjets || sample==kZjets) {
      // Determine the HT bin
      // Assume the sample names are of the form "/path/to/WJets_ht100to200.root"
      //string inputFileNameStr(inputFileName);
      //size_t htWord    = inputFileNameStr.find("ht");
      //size_t toWord    = inputFileNameStr.find("to");
      //size_t lastDot   = inputFileNameStr.find_last_of(".");
      //string htLowStr   = inputFileNameStr.substr(htWord+2, toWord-htWord-2);
      //string htHighStr  = inputFileNameStr.substr(toWord+2, lastDot-toWord-2);
      //htLow   = atoi(htLowStr.c_str());
      //if(htHighStr=="inf"||htHighStr=="Inf")
      //  htHigh=99999;
      //else
      //  htHigh  = atoi(htHighStr.c_str());
      //if(htLow==0 || htHigh==0) {
      //  throw std::runtime_error(Form("Warning: Error parsing the filename \"%s\", probably it is not of the form \"/path/to/WJets_ht100to200.root\" (go fix that)", inputFileName.Data()));
      //  return false;
      //}
      theEWKCorrPars = EWKCorrPars(sample);
      if(theEWKCorrPars.size()<4) {
        throw std::runtime_error("Warning: theEWKCorrPars does not have 4 array elements, make sure function EWKCorrPars is being called correctly!");
        return false;
      }
    } 
  }

  delete dummyTree; // no longer needed
  Long64_t nentries = events->GetEntries();
  if(maxEntries!=0) nentries = TMath::Min((ULong64_t)maxEntries, (ULong64_t)nentries);
  for (Long64_t ientry=0; ientry<nentries; ientry++) {
    if(debug) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries);
    else if(ientry%100000==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries);
    if(useModSplit && ientry%10!=modSplitIndex) continue;

    int nBytesRead=0;
    theCategory=-1; // plot category 
    typeLepSel=99; // 0: mixed e-mu, 1: all mu, 2: all e, 99: undefined
    for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
      jetPts    [iPt][iEta].clear();
      jetPtsUp  [iPt][iEta].clear();
      jetPtsDown[iPt][iEta].clear();
      jetEtas   [iPt][iEta].clear();
      jetCMVAs  [iPt][iEta].clear();
      jetFlavors[iPt][iEta].clear();
    }
 
    // Vectors to keep around in memory
    TVector2 vectorBosonV2;
    TLorentzVector hbbJet1P4, hbbJet2P4, hbbDijetP4;
    
    bLoad(b["runNumber"],ientry);
    bLoad(b["lumiNumber"],ientry);
    bLoad(b["eventNumber"],ientry);
    // Analysis Preselection
    int nLooseLep=0;
    int iTightLep;
    if(selection>=kWHLightFlavorCR && selection<=kWHPresel) { // WH Resolved Category
      {
        
        // Lepton multiplicity
        bLoad(b["nLooseLep"],ientry);
        if(gt.nLooseLep<1) continue; //N_al = 0
        if(debug) printf("Passed lepton multiplicity\n");
        
        // Trigger
        if(sample==kData) {
          bLoad(b["trigger"],ientry);
          if( (gt.trigger & (1<<3 | 1<<1))==0 ) continue;
        }

        bLoad(b["metFilter"],ientry);
        if(gt.metFilter!=1) continue;
        if(debug) printf("Passed MET filters\n");
     
        // Jet multiplicity
        bLoad(b["nJet"],ientry);
        if(debug) printf("hello from line %d (njet=%d)\n",__LINE__,gt.nJet);
        if     (gt.nJet<2) continue;
        if(debug) printf("hello from line %d\n",__LINE__);
        if(useBoostedCategory) { 
          bLoad(b["nFatjet"],ientry);
          bLoad(b["fj1MSD_corr"],ientry);
          bLoad(b["fj1Pt"],ientry);
          bLoad(b["fj1Eta"],ientry);
          bLoad(b["topWBosonPt"],ientry);
          bLoad(b["topWBosonPhi"],ientry);
          bLoad(b["fj1Phi"],ientry);
          float temp_deltaPhiVH = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.fj1Phi));
          if(
            gt.nFatjet != 0 && 
            gt.fj1Pt >= 250 && 
            gt.fj1MSD_corr >= 40 &&
            fabs(gt.fj1Eta) < 2.4 &&
            gt.topWBosonPt >= 250 &&
            temp_deltaPhiVH >= 2.5
          ) continue;
        }
        // Jet kinematics
        bLoad(b["nJot"],ientry);
        bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
        if(debug) printf("hello from line %d\n",__LINE__);
        if(gt.hbbjtidx[0]>gt.nJot || gt.hbbjtidx[1]>gt.nJot) continue; // Bug fix
        if(debug) printf("hello from line %d\n",__LINE__);
        bLoad(b["jetPt"],ientry);
        bLoad(b["jetRegFac"],ientry);
        bLoad(b["jetEta"],ientry);
        bLoad(b["jetPhi"],ientry);
        hbbJet1Pt=gt.jetRegFac[0]*gt.jetPt[gt.hbbjtidx[0]]; hbbJet1Eta=gt.jetEta[gt.hbbjtidx[0]]; hbbJet1Phi=gt.jetPhi[gt.hbbjtidx[0]];
        hbbJet2Pt=gt.jetRegFac[1]*gt.jetPt[gt.hbbjtidx[1]]; hbbJet2Eta=gt.jetEta[gt.hbbjtidx[1]]; hbbJet2Phi=gt.jetPhi[gt.hbbjtidx[1]];
        if(debug) printf("hbb jet1 pt %.2f, jet2 pt %.2f\n", hbbJet1Pt, hbbJet2Pt);
        if(hbbJet1Pt<25 || hbbJet2Pt<25) continue;
        
        bLoad(b["hbbpt_reg"],ientry);
        bLoad(b["hbbeta"],ientry);
        bLoad(b["hbbphi"],ientry);
        bLoad(b["hbbm_reg"],ientry);
        hbbDijetPt = gt.hbbpt_reg;
        hbbDijetMass = gt.hbbm_reg;
        if(debug) printf("hbbDijetPt %.2f, hbbDijetMass %.2f\n", hbbDijetPt, hbbDijetMass); 
        if(hbbDijetPt<50.) continue;
        if(hbbDijetMass<0 || hbbDijetMass>250.) continue;
        if(debug) printf("passed jet kinematics\n");

        // Lepton ID and isolation
        bLoad(b["nLooseElectron"],ientry);
        bLoad(b["nTightElectron"],ientry);
        bLoad(b["nLooseMuon"],ientry);
        bLoad(b["nTightMuon"],ientry);
        bLoad(b["muonSelBit"],ientry);
        bLoad(b["muonPt"],ientry);
        bLoad(b["electronSelBit"],ientry);
        bLoad(b["electronPt"],ientry);
        for(int i=0; i<gt.nLooseMuon; i++) {
          if(
            typeLepSel==99 && 
            (sample!=kData || (gt.trigger & 1<<3)!=0) &&
            (gt.muonSelBit[i]& 1<<3)!=0 && gt.muonPt[i]>25
          ) {
            typeLepSel=1;
            iTightLep=i;
          }
          if((gt.muonSelBit[i]&1<<0)!=0 && gt.muonPt[i]>5) nLooseLep++;
        }
        for(int i=0; i<gt.nLooseElectron; i++) {
          if(
            typeLepSel==99 &&
            (sample!=kData || (gt.trigger & 1<<1)!=0) && 
            (gt.electronSelBit[i]& 1<<6)!=0 && gt.electronPt[i]>30 /*kEleMvaWP80*/
          ) {
            typeLepSel=2;
            iTightLep=i;
          }
          if((gt.electronSelBit[i]&1<<5)!=0 && gt.electronPt[i]>15) nLooseLep++;
        }
        if(typeLepSel!=1 && typeLepSel!=2) continue;
        if(debug) printf("Passed lepton ID/iso multiplicity\n");

        // Lepton kinematics
        if     (typeLepSel==1) {
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
          lepton1Flav   = 13;
          lepton1Charge = gt.muonPdgId[0]>0? -1:1; 
        } else if(typeLepSel==2) {
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
          lepton1Flav   = 11;
          lepton1Charge = gt.electronPdgId[0]>0? -1:1;
        } else continue;
        if(debug) printf("Passed lepton kinematics\n");
        
        // W reconstruction
        bLoad(b["topWBosonPt"],ientry);
        if(gt.topWBosonPt<50.) continue;
        if(debug) printf("passed W reco\n");
      } // end preselection
      
      // Met
      bLoad(b["pfmet"],ientry);
      bLoad(b["pfmetUp"],ientry);
      bLoad(b["pfmetDown"],ientry);
      bLoad(b["pfmetphi"],ientry);
      bLoad(b["pfmetsig"],ientry);
      bLoad(b["puppimetsig"],ientry);
      
      // Jet B-tagging
      bLoad(b["jetCMVA"],ientry);
      bDiscrMax=gt.jetCMVA[gt.hbbjtidx[0]];
      bDiscrMin=gt.jetCMVA[gt.hbbjtidx[1]];
      
      // Top reconstruction
      bLoad(b["topMassLep1Met"],ientry);
      bLoad(b["topMassLep1Met_jesUp"],ientry);
      bLoad(b["topMassLep1Met_jesDown"],ientry);
      bLoad(b["topWBosonPt"],ientry);
      bLoad(b["topWBosonEta"],ientry);
      bLoad(b["topWBosonPhi"],ientry);
      bLoad(b["mT"],ientry);
      
      // mT, W pT up down -- temporary
      { 
        TVector2 metV2Up, metV2Down, lepton1V2, WBosonV2Up, WBosonV2Down;
        metV2Up.SetMagPhi(gt.pfmetUp, gt.pfmetphi);
        metV2Down.SetMagPhi(gt.pfmetDown, gt.pfmetphi);
        lepton1V2.SetMagPhi(lepton1Pt, lepton1Phi);
        mT_jesUp   = TMath::Sqrt(2.*gt.pfmetUp*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        mT_jesDown = TMath::Sqrt(2.*gt.pfmetDown*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        WBosonV2Up=metV2Up+lepton1V2;
        WBosonV2Down=metV2Down+lepton1V2;
        topWBosonPt_jesUp   = WBosonV2Up.Mod();
        topWBosonPt_jesDown = WBosonV2Down.Mod();
        topWBosonPhi_jesUp   = WBosonV2Up.Phi();
        topWBosonPhi_jesDown = WBosonV2Down.Phi();
        //if(debug) printf("topWBosonPt %.1f, topWBosonPt_jesUp %.1f, topWBosonPt_jesDown %.1f\n", topWBosonPt, topWBosonPt_jesUp, topWBosonPt_jesDown);
      }
      // Other stuff
      
      bLoad(b["hbbpt_reg_jesUp"],ientry);
      bLoad(b["hbbpt_reg_jesDown"],ientry);
      bLoad(b["hbbm_reg_jesUp"],ientry);
      bLoad(b["hbbm_reg_jesDown"],ientry);
      hbbDijetPtUp = gt.hbbpt_reg_jesUp;
      hbbDijetPtDown = gt.hbbpt_reg_jesDown;
      hbbDijetMassUp = gt.hbbm_reg_jesUp;
      hbbDijetMassDown = gt.hbbm_reg_jesDown;
      bLoad(b["jetPtUp"],ientry);
      bLoad(b["jetPtDown"],ientry);
      hbbJet1PtUp = gt.jetRegFac[0]*gt.jetPtUp[gt.hbbjtidx[0]];
      hbbJet1PtDown = gt.jetRegFac[0]*gt.jetPtDown[gt.hbbjtidx[0]];
      hbbJet2PtUp = gt.jetRegFac[1]*gt.jetPtUp[gt.hbbjtidx[1]];
      hbbJet2PtDown = gt.jetRegFac[1]*gt.jetPtDown[gt.hbbjtidx[1]];
      bLoad(b["topWBosonCosThetaCS"],ientry);
      bLoad(b["hbbCosThetaJJ"],ientry);
      bLoad(b["hbbCosThetaCSJ1"],ientry);
      
      bLoad(b["hbbphi"],ientry);
      bLoad(b["hbbphi_jesUp"],ientry);
      bLoad(b["hbbphi_jesDown"],ientry);
      deltaPhiVH         = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.hbbphi        ));
      deltaPhiVH_jesUp   = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.hbbphi_jesUp  ));
      deltaPhiVH_jesDown = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.hbbphi_jesDown));
      deltaPhiLep1Met = fabs(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi));

      // Jet multiplicity for JES
      bLoad(b["nJet_jesUp"],ientry);
      bLoad(b["nJet_jesDown"],ientry);
      bLoad(b["sumEtSoft1"],ientry);
      bLoad(b["nSoft2"],ientry);
      bLoad(b["nSoft5"],ientry);
      bLoad(b["nSoft10"],ientry);
      
      bLoad(b["jetGenFlavor"],ientry);

      // CMVA jet kinematic decorrelated weight nuisances
      for(unsigned iJ=0; iJ<(unsigned)gt.nJot; iJ++) {
        int iPt=-1, iEta=-1;
        double jetAbsEta=fabs(gt.jetEta[iJ]);
        if      (gt.jetPt[iJ] >= 19.99 && gt.jetPt[iJ] < 30 ) iPt = 0;
        else if (gt.jetPt[iJ] >= 30    && gt.jetPt[iJ] < 40 ) iPt = 1;
        else if (gt.jetPt[iJ] >= 40    && gt.jetPt[iJ] < 60 ) iPt = 2;
        else if (gt.jetPt[iJ] >= 60    && gt.jetPt[iJ] < 100) iPt = 3;
        else if (gt.jetPt[iJ] >= 100                        ) iPt = 4;
        if      (jetAbsEta >= 0   && jetAbsEta < 0.8  ) iEta = 0;
        else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6  ) iEta = 1;
        else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41 ) iEta = 2;
        
        if(iPt>=0 && iEta>=0) {
          jetPts    [iPt][iEta].push_back(gt.jetPt        [iJ]);
          jetEtas   [iPt][iEta].push_back(gt.jetEta       [iJ]);
          jetCMVAs  [iPt][iEta].push_back(gt.jetCMVA      [iJ]);
          jetFlavors[iPt][iEta].push_back(gt.jetGenFlavor [iJ]);
          if(gt.jetPtUp  [iJ]>=20.) jetPtsUp  [iPt][iEta].push_back(gt.jetPtUp      [iJ]); // Choose iPt based on the varied jet Pt? not sure
          if(gt.jetPtDown[iJ]>=20.) jetPtsDown[iPt][iEta].push_back(gt.jetPtDown    [iJ]);
        }
      }
      
      // End WH Resolved Category 
    } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) { 
      // WH Boosted Category
      {
        // Lepton multiplicity
        bLoad(b["nLooseLep"],ientry);
        if(gt.nLooseLep<1) continue; //N_al = 0
        if(debug) printf("Passed lepton multiplicity\n");
        
        // Trigger
        if(sample==kData) {
          bLoad(b["trigger"],ientry);
          if( (gt.trigger & (1<<3 | 1<<1))==0 ) continue;
        }

        bLoad(b["metFilter"],ientry);
        if(gt.metFilter!=1) continue;
        if(debug) printf("Passed MET filters\n");
     
        // Jet multiplicity
        bLoad(b["nJet"],ientry);
        bLoad(b["nFatjet"],ientry);
        if     (gt.nFatjet==0) continue;
        // Jet kinematics
        bLoad(b["fj1Eta"],ientry);
        if(fabs(gt.fj1Eta)>2.4) continue;
        bLoad(b["fj1ECFN_2_4_20"],ientry);
        if(*((float*)b["fj1ECFN_2_4_20"]->GetAddress()) <= 0) continue;
        if(debug) printf("passed jet kinematics\n");

        // Lepton ID and isolation
        bLoad(b["nLooseElectron"],ientry);
        bLoad(b["nTightElectron"],ientry);
        bLoad(b["electronSelBit"],ientry);
        bLoad(b["electronPt"],ientry);
        bLoad(b["nLooseMuon"],ientry);
        bLoad(b["nTightMuon"],ientry);
        bLoad(b["muonSelBit"],ientry);
        bLoad(b["muonPt"],ientry);
        for(int i=0; i<gt.nLooseMuon; i++) {
          if(
            typeLepSel==99 && 
            (sample!=kData || (gt.trigger & 1<<3)!=0) &&
            (gt.muonSelBit[i]&1<<3)!=0 && gt.muonPt[i]>25
          ) {
            typeLepSel=1;
            iTightLep=i;
          }
          if((gt.muonSelBit[i]&1<<0)!=0 && gt.muonPt[i]>5) nLooseLep++;
          //printf("muon tight=%d, loose=%d, pT=%.1f\n", (gt.muonSelBit[i]&1<<3)!=0, (gt.muonSelBit[i]&1<<0)!=0, gt.muonPt[i]);
        }
        for(int i=0; i<gt.nLooseElectron; i++) {
          if(
            typeLepSel==99 &&
            (sample!=kData || (gt.trigger & 1<<1)!=0) && 
            (gt.electronSelBit[i]& 1<<6)!=0 && gt.electronPt[i]>30 /*kEleMvaWP80*/
          ) {
            typeLepSel=2;
            iTightLep=i;
          }
          if((gt.electronSelBit[i]&1<<5)!=0 && gt.electronPt[i]>15) nLooseLep++;
        }
        if(typeLepSel!=1 && typeLepSel!=2) continue;
        if(debug) printf("Passed lepton ID/iso multiplicity\n");

        // Lepton kinematics
        if     (typeLepSel==1) {
          bLoad(b["muonEta"],ientry);
          bLoad(b["muonPhi"],ientry);
          bLoad(b["muonD0"],ientry);
          bLoad(b["muonDZ"],ientry);
          bLoad(b["muonCombIso"],ientry);
          bLoad(b["muonPdgId"],ientry);
          lepton1Pt     = gt.muonPt[iTightLep];
          lepton1Eta    = gt.muonEta[iTightLep];
          lepton1Phi    = gt.muonPhi[iTightLep]; 
          lepton1D0     = gt.muonD0[iTightLep];
          lepton1DZ     = gt.muonDZ[iTightLep];
          lepton1RelIso = gt.muonCombIso[iTightLep]/gt.muonPt[iTightLep];
          lepton1Flav   = 13;
          lepton1Charge = gt.muonPdgId[iTightLep]>0? -1:1; 
        } else if(typeLepSel==2) {
          bLoad(b["electronEta"],ientry);
          bLoad(b["electronPhi"],ientry);
          bLoad(b["electronD0"],ientry);
          bLoad(b["electronDZ"],ientry);
          bLoad(b["electronCombIso"],ientry);
          bLoad(b["electronPdgId"],ientry);
          lepton1Pt     = gt.electronPt[iTightLep]; 
          lepton1Eta    = gt.electronEta[iTightLep];
          lepton1Phi    = gt.electronPhi[iTightLep];
          lepton1RelIso = gt.electronCombIso[iTightLep]/gt.electronPt[iTightLep];
          lepton1D0     = gt.electronD0[iTightLep];
          lepton1DZ     = gt.electronDZ[iTightLep];
          lepton1Flav   = 11;
          lepton1Charge = gt.electronPdgId[iTightLep]>0? -1:1;
        } else continue;
        if(debug) printf("Passed lepton kinematics\n");
        
        // W reconstruction
        bLoad(b["topWBosonPt"],ientry);
        if(gt.topWBosonPt<50.) continue;
        if(debug) printf("passed W reco\n");
      } // end preselection
      
      // Met
      bLoad(b["pfmet"],ientry);
      bLoad(b["pfmetUp"],ientry);
      bLoad(b["pfmetDown"],ientry);
      bLoad(b["pfmetphi"],ientry);
      bLoad(b["pfmetsig"],ientry);
      bLoad(b["puppimetsig"],ientry);
      
      // W reconstruction
      bLoad(b["topWBosonPt"],ientry);
      bLoad(b["topWBosonPhi"],ientry);
      bLoad(b["mT"],ientry);
      
      // mT, W pT up down
      { 
        TVector2 metV2Up, metV2Down, lepton1V2, WBosonV2Up, WBosonV2Down;
        metV2Up.SetMagPhi(gt.pfmetUp, gt.pfmetphi);
        metV2Down.SetMagPhi(gt.pfmetDown, gt.pfmetphi);
        lepton1V2.SetMagPhi(lepton1Pt, lepton1Phi);
        mT_jesUp   = TMath::Sqrt(2.*gt.pfmetUp*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        mT_jesDown = TMath::Sqrt(2.*gt.pfmetDown*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        WBosonV2Up=metV2Up+lepton1V2;
        WBosonV2Down=metV2Down+lepton1V2;
        topWBosonPt_jesUp   = WBosonV2Up.Mod();
        topWBosonPt_jesDown = WBosonV2Down.Mod();
        topWBosonPhi_jesUp   = WBosonV2Up.Phi();
        topWBosonPhi_jesDown = WBosonV2Down.Phi();
        //if(debug) printf("topWBosonPt %.1f, topWBosonPt_jesUp %.1f, topWBosonPt_jesDown %.1f\n", topWBosonPt, topWBosonPt_jesUp, topWBosonPt_jesDown);
      }
      
      // Other stuff
      bLoad(b["fj1Phi"],ientry);
      deltaPhiVH         = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.fj1Phi));
      deltaPhiVH_jesUp   = fabs(TVector2::Phi_mpi_pi(topWBosonPhi_jesUp  -gt.fj1Phi));
      deltaPhiVH_jesDown = fabs(TVector2::Phi_mpi_pi(topWBosonPhi_jesDown-gt.fj1Phi));
      deltaPhiLep1Met = fabs(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi));
      // End WH Boosted Category
    }
    if((selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel)) {
      bLoad(b["nFatjet"],ientry);
      bLoad(b["fj1Eta"],ientry);
      bLoad(b["fj1MSD"],ientry);
      bLoad(b["fj1MSD_corr"],ientry);
      // Fatjet Properties
      bLoad(b["fj1Tau32"],ientry);
      bLoad(b["fj1Tau21"],ientry);
      bLoad(b["fj1Tau32SD"],ientry); 
      bLoad(b["fj1Tau21SD"],ientry); 
      bLoad(b["fj1MSDScaleUp"],ientry);    
      bLoad(b["fj1MSDScaleDown"],ientry);      
      bLoad(b["fj1MSDSmeared"],ientry);    
      bLoad(b["fj1MSDSmearedUp"],ientry);      
      bLoad(b["fj1MSDSmearedDown"],ientry);        
      bLoad(b["fj1PtScaleUp"],ientry);   
      bLoad(b["fj1PtScaleDown"],ientry);     
      bLoad(b["fj1PtSmeared"],ientry);   
      bLoad(b["fj1PtSmearedUp"],ientry);     
      bLoad(b["fj1PtSmearedDown"],ientry);       
      bLoad(b["fj1HighestPtGen"],ientry);
      bLoad(b["fj1Pt"],ientry);
      bLoad(b["fj1Phi"],ientry);
      bLoad(b["fj1M"],ientry);
      bLoad(b["fj1MaxCSV"],ientry);
      bLoad(b["fj1MinCSV"],ientry);
      bLoad(b["fj1DoubleCSV"],ientry);   
      bLoad(b["fj1HTTMass"],ientry); 
      bLoad(b["fj1HTTFRec"],ientry); 
      bLoad(b["fj1SDEFrac100"],ientry);    
       
      fj1MSD_corr         = gt.fj1MSD_corr;
      fj1MSD_corr_jesUp   = fj1MSD_corr * gt.fj1MSDScaleUp  /gt.fj1MSD_corr;
      fj1MSD_corr_jesDown = fj1MSD_corr * gt.fj1MSDScaleDown/gt.fj1MSD_corr;
      // Ak4 jet properties
      bLoad(b["jetPt"],ientry);
      bLoad(b["jetEta"],ientry);
      bLoad(b["jetPhi"],ientry);
      bLoad(b["jetPtUp"],ientry);
      bLoad(b["jetPtDown"],ientry);
      bLoad(b["jetCMVA"],ientry);
      bLoad(b["jetIso"],ientry); 
      bLoad(b["isojetNBtags"],ientry);

      isojets.clear();
      isojets_jesUp.clear();
      isojets_jesDown.clear();
      isojetNBtags=0, isojetNBtags_jesUp=0, isojetNBtags_jesDown=0;
      for(unsigned char iJ=0; iJ<gt.nJot; iJ++) {
        //if(!gt.jetIso[iJ]) continue; 
        if(fabs(gt.jetEta[iJ])>2.4) continue;
        float dR2JetFatjet=pow(gt.jetEta[iJ]-gt.fj1Eta,2)+pow(TVector2::Phi_mpi_pi(gt.jetPhi[iJ]-gt.fj1Phi),2);
        if(dR2JetFatjet<0.64) continue;

        if(gt.jetPt[iJ]>30) {
          isojets.push_back(iJ);
          if(gt.jetCMVA[iJ]>bDiscrLoose) isojetNBtags++;
        } if(gt.jetPtUp[iJ]>30) {
          isojets_jesUp.push_back(iJ);
          if(gt.jetCMVA[iJ]>bDiscrLoose) isojetNBtags_jesUp++;
        } if(gt.jetPtDown[iJ]>30) {
          isojets_jesDown.push_back(iJ);
          if(gt.jetCMVA[iJ]>bDiscrLoose) isojetNBtags_jesDown++;
        }
        
        // CMVA jet kinematic decorrelated weight nuisances for the isojets only
        int iPt=-1, iEta=-1;
        double jetAbsEta=fabs(gt.jetEta[iJ]);
        //if      (gt.jetPt[iJ] >= 19.99 && gt.jetPt[iJ] < 30 ) iPt = 0;
        //else if (gt.jetPt[iJ] >= 30    && gt.jetPt[iJ] < 40 ) iPt = 1;
        if      (gt.jetPt[iJ] >= 30    && gt.jetPt[iJ] < 40 ) iPt = 1;
        else if (gt.jetPt[iJ] >= 40    && gt.jetPt[iJ] < 60 ) iPt = 2;
        else if (gt.jetPt[iJ] >= 60    && gt.jetPt[iJ] < 100) iPt = 3;
        else if (gt.jetPt[iJ] >= 100                        ) iPt = 4;
        if      (jetAbsEta >= 0   && jetAbsEta < 0.8  ) iEta = 0;
        else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6  ) iEta = 1;
        else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41 ) iEta = 2;
        
        if(iPt>=0 && iEta>=0) {
          jetPts    [iPt][iEta].push_back(gt.jetPt        [iJ]);
          jetEtas   [iPt][iEta].push_back(gt.jetEta       [iJ]);
          jetCMVAs  [iPt][iEta].push_back(gt.jetCMVA      [iJ]);
          jetFlavors[iPt][iEta].push_back(gt.jetGenFlavor [iJ]);
          if(gt.jetPtUp  [iJ]>=30.) jetPtsUp  [iPt][iEta].push_back(gt.jetPtUp      [iJ]); // Choose iPt based on the varied jet Pt? not sure
          if(gt.jetPtDown[iJ]>=30.) jetPtsDown[iPt][iEta].push_back(gt.jetPtDown    [iJ]);
        }
      }
      //isojetNBtags=gt.isojetNBtags;
      nIsojet=isojets.size();
      nIsojet_jesUp=isojets_jesUp.size();
      nIsojet_jesDown=isojets_jesDown.size();
 
      // Energy Correlation Functions
      bLoad(b["fj1ECFN_1_1_05"],ientry);
      bLoad(b["fj1ECFN_2_1_05"],ientry);
      bLoad(b["fj1ECFN_3_1_05"],ientry);
      bLoad(b["fj1ECFN_1_2_05"],ientry);
      bLoad(b["fj1ECFN_2_2_05"],ientry);
      bLoad(b["fj1ECFN_3_2_05"],ientry);
      bLoad(b["fj1ECFN_1_3_05"],ientry);
      bLoad(b["fj1ECFN_2_3_05"],ientry);
      bLoad(b["fj1ECFN_3_3_05"],ientry);
      bLoad(b["fj1ECFN_1_4_05"],ientry);
      bLoad(b["fj1ECFN_2_4_05"],ientry);
      bLoad(b["fj1ECFN_3_4_05"],ientry);
      bLoad(b["fj1ECFN_1_1_10"],ientry);
      bLoad(b["fj1ECFN_2_1_10"],ientry);
      bLoad(b["fj1ECFN_3_1_10"],ientry);
      bLoad(b["fj1ECFN_1_2_10"],ientry);
      bLoad(b["fj1ECFN_2_2_10"],ientry);
      bLoad(b["fj1ECFN_3_2_10"],ientry);
      bLoad(b["fj1ECFN_1_3_10"],ientry);
      bLoad(b["fj1ECFN_2_3_10"],ientry);
      bLoad(b["fj1ECFN_3_3_10"],ientry);
      bLoad(b["fj1ECFN_1_4_10"],ientry);
      bLoad(b["fj1ECFN_2_4_10"],ientry);
      bLoad(b["fj1ECFN_3_4_10"],ientry);
      bLoad(b["fj1ECFN_1_1_20"],ientry);
      bLoad(b["fj1ECFN_2_1_20"],ientry);
      bLoad(b["fj1ECFN_3_1_20"],ientry);
      bLoad(b["fj1ECFN_1_2_20"],ientry);
      bLoad(b["fj1ECFN_2_2_20"],ientry);
      bLoad(b["fj1ECFN_3_2_20"],ientry);
      bLoad(b["fj1ECFN_1_3_20"],ientry);
      bLoad(b["fj1ECFN_2_3_20"],ientry);
      bLoad(b["fj1ECFN_3_3_20"],ientry);
      bLoad(b["fj1ECFN_1_4_20"],ientry);
      bLoad(b["fj1ECFN_2_4_20"],ientry);
      bLoad(b["fj1ECFN_3_4_20"],ientry);
      bLoad(b["fj1ECFN_1_1_40"],ientry);
      bLoad(b["fj1ECFN_2_1_40"],ientry);
      bLoad(b["fj1ECFN_3_1_40"],ientry);
      bLoad(b["fj1ECFN_1_2_40"],ientry);
      bLoad(b["fj1ECFN_2_2_40"],ientry);
      bLoad(b["fj1ECFN_3_2_40"],ientry);
      bLoad(b["fj1ECFN_1_3_40"],ientry);
      bLoad(b["fj1ECFN_2_3_40"],ientry);
      bLoad(b["fj1ECFN_3_3_40"],ientry);
      bLoad(b["fj1ECFN_1_4_40"],ientry);
      bLoad(b["fj1ECFN_2_4_40"],ientry);
      bLoad(b["fj1ECFN_3_4_40"],ientry);
    }
    // Set Selection Bits
    selectionBits=0; nMinusOneBits=0; selectionBits_jesUp=0; selectionBits_jesDown=0;
    std::map<TString, bool> cut, cut_jesUp, cut_jesDown;
    
    if(selection>=kWHLightFlavorCR && selection<=kWHPresel) { // Begin WH Resolved Selection

      cut["2ndLepVeto" ] = nLooseLep==1;
      cut["ultraLepIso"] = lepton1RelIso<0.06;
      cut["WpT"        ] = gt.topWBosonPt>100;
      cut["pTjj"       ] = hbbDijetPt>100; 
      cut["lepton1IP"  ] = (typeLepSel==1 && lepton1D0<0.20 && lepton1DZ<0.50) || (typeLepSel==2 && ((fabs(lepton1Eta)<1.479 && lepton1D0<0.05 && lepton1DZ<0.10) || (fabs(lepton1Eta)>=1.479 && lepton1D0<0.10 && lepton1DZ<0.20)));
      cut["dPhiVH"     ] = deltaPhiVH > 2.5;
      cut["dPhiLep1Met"] = deltaPhiLep1Met < 2;
      cut["tightBTag"  ] = (bDiscrMax >= bDiscrTight);
      cut["mediumBVeto"] = (bDiscrMax < bDiscrMedium);
      cut["looseBTag"  ] = (bDiscrMax >= bDiscrLoose);
      cut["looseBTag2" ] = (bDiscrMin >= bDiscrLoose);
      cut["metSig"     ] = gt.pfmetsig>2;
      cut["nJet_WHTT"  ] = gt.nJet>=4;
      cut["nJet_WHHF"  ] = gt.nJet==2;
      cut["nJet_WHSR"  ] = gt.nJet<4;
      cut["mH_WHSR"    ] = ((hbbDijetMass >= 90) && (hbbDijetMass <  150));
      cut["mH_Flip"    ] = ((hbbDijetMass <  90) || (hbbDijetMass >= 150));

      vector<TString> cutsWHLightFlavorCR, cutsWHHeavyFlavorCR, cutsWH2TopCR, cutsWHSR;
      cutsWHLightFlavorCR ={"lepton1IP", "ultraLepIso", "2ndLepVeto","WpT","pTjj",         "dPhiLep1Met","looseBTag","mediumBVeto","metSig"};
      cutsWHHeavyFlavorCR ={"lepton1IP", "ultraLepIso", "2ndLepVeto","WpT","pTjj",         "dPhiLep1Met","nJet_WHHF","tightBTag","mH_Flip","metSig"};
      cutsWH2TopCR        ={"lepton1IP", "ultraLepIso", "2ndLepVeto","WpT","pTjj",         "dPhiLep1Met","nJet_WHTT","tightBTag"};
      cutsWHSR            ={"lepton1IP", "ultraLepIso", "2ndLepVeto","WpT","pTjj","dPhiVH","dPhiLep1Met","nJet_WHSR","tightBTag","looseBTag2","mH_WHSR"};

      if(passAllCuts( cut, cutsWHLightFlavorCR))   selectionBits |= kWHLightFlavorCR;
      if(passAllCuts( cut, cutsWHHeavyFlavorCR))   selectionBits |= kWHHeavyFlavorCR;
      if(passAllCuts( cut, cutsWH2TopCR       ))   selectionBits |= kWH2TopCR;
      if(passAllCuts( cut, cutsWHSR           ))   selectionBits |= kWHSR;
      
      if(passNMinusOne( cut, cutsWHLightFlavorCR)) nMinusOneBits |= kWHLightFlavorCR;
      if(passNMinusOne( cut, cutsWHHeavyFlavorCR)) nMinusOneBits |= kWHHeavyFlavorCR;
      if(passNMinusOne( cut, cutsWH2TopCR       )) nMinusOneBits |= kWH2TopCR;
      if(passNMinusOne( cut, cutsWHSR           )) nMinusOneBits |= kWHSR;

      cut_jesUp=cut;
      cut_jesUp["WpT"        ] = topWBosonPt_jesUp>100;
      cut_jesUp["pTjj"       ] = hbbDijetPtUp>100; 
      cut_jesUp["dPhiVH"     ] = deltaPhiVH_jesUp > 2.5;
      cut_jesUp["dPhiLep1Met"] = deltaPhiLep1Met < 2;
      cut_jesUp["nJet_WHTT"  ] = gt.nJet_jesUp>=4;
      cut_jesUp["nJet_WHHF"  ] = gt.nJet_jesUp==2;
      cut_jesUp["nJet_WHSR"  ] = gt.nJet_jesUp<4;
      cut_jesUp["mH_WHSR"    ] = ((hbbDijetMassUp >= 90) && (hbbDijetMassUp <  150));
      cut_jesUp["mH_Flip"    ] = ((hbbDijetMassUp < 90) || (hbbDijetMassUp >= 150 && hbbDijetMassUp < 250));
      if(passAllCuts( cut_jesUp, cutsWHLightFlavorCR, debug)) selectionBits_jesUp |= kWHLightFlavorCR;
      if(passAllCuts( cut_jesUp, cutsWHHeavyFlavorCR, debug)) selectionBits_jesUp |= kWHHeavyFlavorCR;
      if(passAllCuts( cut_jesUp, cutsWH2TopCR       , debug)) selectionBits_jesUp |= kWH2TopCR;
      if(passAllCuts( cut_jesUp, cutsWHSR           , debug)) selectionBits_jesUp |= kWHSR;
      
      cut_jesDown=cut;
      cut_jesDown["WpT"        ] = topWBosonPt_jesDown>100;
      cut_jesDown["pTjj"       ] = hbbDijetPtDown>100; 
      cut_jesDown["dPhiVH"     ] = deltaPhiVH_jesDown > 2.5;
      cut_jesDown["dPhiLep1Met"] = deltaPhiLep1Met < 2;
      cut_jesDown["nJet_WHTT"  ] = gt.nJet_jesDown>=4;
      cut_jesDown["nJet_WHHF"  ] = gt.nJet_jesDown==2;
      cut_jesDown["nJet_WHSR"  ] = gt.nJet_jesDown<4;
      cut_jesDown["mH_WHSR"    ] = ((hbbDijetMassDown >= 90) && (hbbDijetMassDown <  150));
      cut_jesDown["mH_Flip"    ] = ((hbbDijetMassDown < 90) || (hbbDijetMassDown >= 150 && hbbDijetMassDown < 250));
      if(passAllCuts( cut_jesDown, cutsWHLightFlavorCR, debug)) selectionBits_jesDown |= kWHLightFlavorCR;
      if(passAllCuts( cut_jesDown, cutsWHHeavyFlavorCR, debug)) selectionBits_jesDown |= kWHHeavyFlavorCR;
      if(passAllCuts( cut_jesDown, cutsWH2TopCR       , debug)) selectionBits_jesDown |= kWH2TopCR;
      if(passAllCuts( cut_jesDown, cutsWHSR           , debug)) selectionBits_jesDown |= kWHSR;
      // End WH Resolved Selection

    } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
      // Begin WH Boosted Selection
      float MSDcutoff=40;
      float MSDmin=80, MSDmax=150;
      cut["2ndLepVeto" ] = nLooseLep==1;
      cut["ultraLepIso"] = lepton1RelIso<0.06;
      cut["WpT"        ] = gt.topWBosonPt>250;
      cut["pTfj"       ] = gt.fj1Pt>250; 
      cut["lepton1Pt"  ] = ((typeLepSel==1 && lepton1Pt>25) || (typeLepSel==2 && lepton1Pt>30));
      cut["lepton1IP"  ] = (typeLepSel==1 && lepton1D0<0.20 && lepton1DZ<0.50) || (typeLepSel==2 && ((fabs(lepton1Eta)<1.479 && lepton1D0<0.05 && lepton1DZ<0.10) || (fabs(lepton1Eta)>=1.479 && lepton1D0<0.10 && lepton1DZ<0.20)));
      cut["dPhiVH"     ] = deltaPhiVH>2.5;
      cut["DoubleBTag" ] = gt.fj1DoubleCSV>=0.8;
      cut["DoubleBVeto"] = gt.fj1DoubleCSV<0.8;
      cut["mH_WHFJSR"  ] = (fj1MSD_corr>=MSDmin && fj1MSD_corr<MSDmax);
      cut["mH_WHHFCR"  ] = fj1MSD_corr<MSDmin;
      cut["mH_cutoff"  ] = fj1MSD_corr>MSDcutoff;
      cut["isojetBtag" ] = isojetNBtags>0;
      cut["isojetBVeto"] = isojetNBtags==0;
      //cut["metSig"     ] = gt.pfmetsig>2;
      
      vector<TString> cutsWHLightFlavorFJCR, cutsWHHeavyFlavorFJCR, cutsWHTT2bFJCR, cutsWHTT1bFJCR, cutsWHFJSR, cutsWHFJPresel;
      cutsWHLightFlavorFJCR ={"ecfSanity","ultraLepIso", "lepton1IP", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","DoubleBVeto","isojetBVeto"            };
      cutsWHHeavyFlavorFJCR ={"ecfSanity","ultraLepIso", "lepton1IP", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","DoubleBTag" ,"isojetBVeto","mH_WHHFCR"};
      cutsWHTT2bFJCR        ={"ecfSanity","ultraLepIso", "lepton1IP", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","DoubleBTag" ,"isojetBtag"             };
      cutsWHTT1bFJCR        ={"ecfSanity","ultraLepIso", "lepton1IP", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","DoubleBVeto","isojetBtag"             };
      cutsWHFJSR            ={"ecfSanity","ultraLepIso", "lepton1IP", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","DoubleBTag" ,"isojetBVeto","mH_WHFJSR"};
      cutsWHFJPresel        ={"ecfSanity","ultraLepIso", "lepton1IP", "2ndLepVeto","WpT","pTfj","lepton1Pt"                                                 };

      if(passAllCuts( cut, cutsWHLightFlavorFJCR, debug))   selectionBits |= kWHLightFlavorFJCR;
      if(passAllCuts( cut, cutsWHHeavyFlavorFJCR, debug))   selectionBits |= kWHHeavyFlavorFJCR;
      if(passAllCuts( cut, cutsWHTT2bFJCR       , debug))   selectionBits |= kWHTT2bFJCR;
      if(passAllCuts( cut, cutsWHTT1bFJCR       , debug))   selectionBits |= kWHTT1bFJCR;
      if(passAllCuts( cut, cutsWHFJSR           , debug))   selectionBits |= kWHFJSR;
      if(passAllCuts( cut, cutsWHFJPresel       , debug))   selectionBits |= kWHFJPresel;
      if(passNMinusOne( cut, cutsWHLightFlavorFJCR)) nMinusOneBits |= kWHLightFlavorFJCR;
      if(passNMinusOne( cut, cutsWHHeavyFlavorFJCR)) nMinusOneBits |= kWHHeavyFlavorFJCR;
      if(passNMinusOne( cut, cutsWHTT2bFJCR       )) nMinusOneBits |= kWHTT2bFJCR;
      if(passNMinusOne( cut, cutsWHTT1bFJCR       )) nMinusOneBits |= kWHTT1bFJCR;
      if(passNMinusOne( cut, cutsWHFJSR           )) nMinusOneBits |= kWHFJSR;
      if(passNMinusOne( cut, cutsWHFJPresel       )) nMinusOneBits |= kWHFJPresel;

      cut_jesUp=cut;
      cut_jesUp["WpT"        ] = topWBosonPt_jesUp>250;
      cut_jesUp["pTfj"       ] = gt.fj1PtScaleUp>250; 
      cut_jesUp["dPhiVH"     ] = deltaPhiVH_jesUp>2.5;
      cut_jesUp["DoubleBTag" ] = gt.fj1DoubleCSV>=0.8;
      cut_jesUp["DoubleBVeto"] = gt.fj1DoubleCSV<0.8;
      cut_jesUp["mH_WHFJSR"  ] = ((fj1MSD_corr_jesUp>=MSDmin));
      cut_jesUp["mH_WHHFCR"  ] = fj1MSD_corr_jesUp<MSDmin;
      cut_jesUp["mH_cutoff"  ] = fj1MSD_corr_jesUp>MSDcutoff;
      cut_jesUp["isojetBtag" ] = isojetNBtags_jesUp>0;
      cut_jesUp["isojetBVeto"] = isojetNBtags_jesUp==0;
      if(passAllCuts( cut_jesUp, cutsWHLightFlavorFJCR)) selectionBits_jesUp |= kWHLightFlavorFJCR;
      if(passAllCuts( cut_jesUp, cutsWHHeavyFlavorFJCR)) selectionBits_jesUp |= kWHHeavyFlavorFJCR;
      if(passAllCuts( cut_jesUp, cutsWHTT2bFJCR       )) selectionBits_jesUp |= kWHTT2bFJCR;
      if(passAllCuts( cut_jesUp, cutsWHTT1bFJCR       )) selectionBits_jesUp |= kWHTT1bFJCR;
      if(passAllCuts( cut_jesUp, cutsWHFJSR           )) selectionBits_jesUp |= kWHFJSR;
      
      cut_jesDown=cut;
      cut_jesDown["WpT"        ] = topWBosonPt_jesDown>250;
      cut_jesDown["pTfj"       ] = gt.fj1PtScaleDown>250; 
      cut_jesDown["dPhiVH"     ] = deltaPhiVH_jesDown>2.5;
      cut_jesDown["DoubleBTag" ] = gt.fj1DoubleCSV>=0.8;
      cut_jesDown["DoubleBVeto"] = gt.fj1DoubleCSV<0.8;
      cut_jesDown["mH_WHFJSR"  ] = ((fj1MSD_corr_jesDown>=MSDmin));
      cut_jesDown["mH_WHHFCR"  ] = fj1MSD_corr_jesDown<MSDmin;
      cut_jesDown["mH_cutoff"  ] = fj1MSD_corr_jesDown>MSDcutoff;
      cut_jesDown["isojetBtag" ] = isojetNBtags_jesDown>0;
      cut_jesDown["isojetBVeto"] = isojetNBtags_jesDown==0;
      if(passAllCuts( cut_jesDown, cutsWHLightFlavorFJCR)) selectionBits_jesDown |= kWHLightFlavorFJCR;
      if(passAllCuts( cut_jesDown, cutsWHHeavyFlavorFJCR)) selectionBits_jesDown |= kWHHeavyFlavorFJCR;
      if(passAllCuts( cut_jesDown, cutsWHTT2bFJCR       )) selectionBits_jesDown |= kWHTT2bFJCR;
      if(passAllCuts( cut_jesDown, cutsWHTT1bFJCR       )) selectionBits_jesDown |= kWHTT1bFJCR;
      if(passAllCuts( cut_jesDown, cutsWHFJSR           )) selectionBits_jesDown |= kWHFJSR;
      // End WH Boosted Selection 
    }
    if(nMinusOneBits==0 && selectionBits==0 && selectionBits_jesUp==0 && selectionBits_jesDown==0) continue;

    // Weighting
    if(sample==kData) weight=1;
    else {
      bLoad(b["normalizedWeight"],ientry);
      bLoad(b["sf_npv"],ientry);
      if(sample==kWjets || sample==kZjets) {
        if(useHtBinnedVJetsKFactor) {
          bLoad(b["trueGenBosonPt"],ientry);
          bLoad(b["lhe_HT"],ientry);
          gt.sf_qcdV=qcdKFactor(sample, gt.lhe_HT);
          gt.sf_ewkV=theEWKCorrPars[0]+theEWKCorrPars[1]*(TMath::Power((gt.trueGenBosonPt+theEWKCorrPars[2]),theEWKCorrPars[3]));
        } else {
          bLoad(b["sf_qcdV"],ientry);
          bLoad(b["sf_ewkV"],ientry);
        }
      } else { 
        gt.sf_qcdV=1;
        gt.sf_ewkV=1;
      }
      
      // Compute Stitching weight for W+jets b quark and hadron enriched samples
      //if(sample==kTT) bLoad(b["sf_tt"],ientry); else gt.sf_tt=1;
      //weight = normalizedWeight * theLumi * gt.sf_npv * gt.sf_ewkV * gt.sf_qcdV * gt.sf_tt;
      weight = normalizedWeight * theLumi * gt.sf_npv * gt.sf_ewkV * gt.sf_qcdV;
      if(sample==kWjets) {
        float stitchWeight=1.;
        bLoad(b["trueGenBosonPt"],ientry);
        bLoad(b["nB"],ientry);
        bLoad(b["nBGenJets"],ientry);
        bool bQuarkEnriched  = gt.nB > 0;
        bool bHadronEnriched = gt.nBGenJets>0 && gt.nB==0;
        if((!bQuarkEnriched && !bHadronEnriched) || gt.trueGenBosonPt<100) {
          stitchWeight=1;
        } else if(bQuarkEnriched) {
          if(gt.trueGenBosonPt>=100 && gt.trueGenBosonPt<200) 
            stitchWeight = vhbbPlot::BenrichedVPT100;
          else if(gt.trueGenBosonPt>=200)
            stitchWeight = vhbbPlot::BenrichedVPT200;
        } else if(bHadronEnriched) {
          if(gt.trueGenBosonPt>=100 && gt.trueGenBosonPt<200) 
            stitchWeight = vhbbPlot::BfilterVPT100;
          else if(gt.trueGenBosonPt>=200)
            stitchWeight = vhbbPlot::BfilterVPT200;
        }
        weight *= stitchWeight;
      }
      if (typeLepSel==1) {
        bLoad(b["muonSfReco"],ientry);
        bLoad(b["muonSfTight"],ientry);
        bLoad(b["muonSfUnc"],ientry);
        bLoad(b["sf_muTrig"],ientry);
        weight *= gt.sf_muTrig * gt.muonSfReco[0] * gt.muonSfTight[0];
        //weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
        // need muon trigger efficiency
        //double theTrigEff = muTrigEff->GetBinContent( muTrigEff->FindBin( fabs(lepton1Eta), TMath::Max((float)26., TMath::Min((float)500., lepton1Pt)) ) );
        //weight *= theTrigEff;
      } else if(typeLepSel==2) {
        bLoad(b["electronSfReco"],ientry);
        bLoad(b["electronSfTight"],ientry);
        //bLoad(b["electronSfMvaWP80"],ientry);
        bLoad(b["electronSfUnc"],ientry);
        bLoad(b["sf_eleTrig"],ientry);
        weight *= gt.sf_eleTrig * gt.electronSfReco[0] * gt.electronSfTight[0];
        //weight *= gt.sf_eleTrig * gt.electronSfReco[0] * gt.electronSfMvaWP80[0];
      }
      float recorrect_vhEWK=1, recorrect_vhEWKUp=1, recorrect_vhEWKDown=1;
      if(sample==kVZ) {
        bLoad(b["sf_wz"],ientry);
        bLoad(b["sf_zz"],ientry);
        bLoad(b["sf_zzUnc"],ientry);
        weight *= gt.sf_wz * gt.sf_zz;
      } else if(sample==kVH) {
        bLoad(b["sf_vh"],ientry);
        bLoad(b["sf_vhUp"],ientry);
        bLoad(b["sf_vhDown"],ientry);
        TString vhChannel=lepton1Charge>0?"WplusH":"WminusH";
        recorrect_vhEWK     = vhEWKCorr(vhChannel,gt.sf_vh    );
        recorrect_vhEWKUp   = vhEWKCorr(vhChannel,gt.sf_vhUp  );
        recorrect_vhEWKDown = vhEWKCorr(vhChannel,gt.sf_vhDown);
        weight *= recorrect_vhEWK;
      }
      
      if(selection>=kWHLightFlavorCR && selection<=kWHPresel) {
        bLoad(b["sf_cmvaWeight_Cent"],ientry);
        weight *= gt.sf_csvWeights[GeneralTree::csvCent];
      }
      if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
        bLoad(b["sf_sjbtag0"],ientry);
        bLoad(b["sf_sjbtagGT0"],ientry);
        //bLoad(b["sf_sjbtag2"],ientry);
        if(selection==kWHLightFlavorFJCR) weight *= (*sf_sjbtag0);
        else                              weight *= (*sf_sjbtagGT0);
        for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
          double cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
          weight *= cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCent, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
        }
      }
      
      // #############################
      // # Variations of the Weights #
      // #############################
      
      // PDF and QCD scale
      bLoad(b["pdfUp"],ientry);
      bLoad(b["pdfDown"],ientry);
      bLoad(b["scale"],ientry);
      weight_pdfUp = weight * gt.pdfUp;
      weight_pdfDown = weight * gt.pdfDown;
      weight_QCDr1f2 = weight * gt.scale[0] * sfLHEweights_QCDr1f2;
      weight_QCDr1f5 = weight * gt.scale[1] * sfLHEweights_QCDr1f5;
      weight_QCDr2f1 = weight * gt.scale[2] * sfLHEweights_QCDr2f1;
      weight_QCDr2f2 = weight * gt.scale[3] * sfLHEweights_QCDr2f2;
      weight_QCDr5f1 = weight * gt.scale[4] * sfLHEweights_QCDr5f1;
      weight_QCDr5f5 = weight * gt.scale[5] * sfLHEweights_QCDr5f5;
      if(sample==kWjets || sample==kZjets) {
        bLoad(b["genBosonPt"],ientry);
        if(sample==kWjets) {
          weight_NLOQCDfUp   = weight * hWjets_relErr_QCDfScaleUp  ->GetBinContent(hWjets_relErr_QCDfScaleUp  ->FindBin(gt.genBosonPt));  
          weight_NLOQCDfDown = weight * hWjets_relErr_QCDfScaleDown->GetBinContent(hWjets_relErr_QCDfScaleDown->FindBin(gt.genBosonPt));  
          weight_NLOQCDrUp   = weight * hWjets_relErr_QCDrScaleUp  ->GetBinContent(hWjets_relErr_QCDrScaleUp  ->FindBin(gt.genBosonPt)); 
          weight_NLOQCDrDown = weight * hWjets_relErr_QCDrScaleDown->GetBinContent(hWjets_relErr_QCDrScaleDown->FindBin(gt.genBosonPt));  
        } else {
          weight_NLOQCDfUp   = weight * hZjets_relErr_QCDfScaleUp  ->GetBinContent(hZjets_relErr_QCDfScaleUp  ->FindBin(gt.genBosonPt));  
          weight_NLOQCDfDown = weight * hZjets_relErr_QCDfScaleDown->GetBinContent(hZjets_relErr_QCDfScaleDown->FindBin(gt.genBosonPt));  
          weight_NLOQCDrUp   = weight * hZjets_relErr_QCDrScaleUp  ->GetBinContent(hZjets_relErr_QCDrScaleUp  ->FindBin(gt.genBosonPt)); 
          weight_NLOQCDrDown = weight * hZjets_relErr_QCDrScaleDown->GetBinContent(hZjets_relErr_QCDrScaleDown->FindBin(gt.genBosonPt));  
        }
      }
      weight_VHCorrUp   = (sample==kVH)? weight * recorrect_vhEWKUp   / recorrect_vhEWK : weight;
      weight_VHCorrDown = (sample==kVH)? weight * recorrect_vhEWKDown / recorrect_vhEWK : weight;
      
      if(selection>=kWHLightFlavorCR && selection<=kWHFJSR) { // WH weights
        // Lepton ID & ISO SF uncertainty
        if(typeLepSel==1) weight_lepSFUp = weight * (1.+gt.muonSfUnc[0]);
        else weight_lepSFUp = weight * (1.+gt.electronSfUnc[0]);
      }
      if(selection>=kWHLightFlavorCR && selection<=kWHFJPresel) {
        // CMVA weight uncertainty
        // https://cmssdt.cern.ch/lxr/source/PhysicsTools/Heppy/python/physicsutils/BTagWeightCalculator.py?v=CMSSW_8_0_20
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors
        for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
          double cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
          double centralWeight = cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCent, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
          weight_cmvaJESUp       [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvJESup       , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaLFUp        [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvLFup        , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaHFUp        [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvHFup        , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaHFStats1Up  [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvHFStats1up  , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaHFStats2Up  [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvHFStats2up  , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaLFStats1Up  [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvLFStats1up  , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaLFStats2Up  [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvLFStats2up  , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaCErr1Up     [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCErr1up     , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaCErr2Up     [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCErr2up     , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaJESDown     [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvJESdown     , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaLFDown      [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvLFdown      , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaHFDown      [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvHFdown      , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaHFStats1Down[iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvHFStats1down, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaHFStats2Down[iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvHFStats2down, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaLFStats1Down[iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvLFStats1down, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaLFStats2Down[iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvLFStats2down, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaCErr1Down   [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCErr1down   , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
          weight_cmvaCErr2Down   [iPt][iEta] = weight*cmvaReweighter->getCSVWeight(jetPts[iPt][iEta], jetEtas[iPt][iEta], jetCMVAs[iPt][iEta], jetFlavors[iPt][iEta], GeneralTree::csvCErr2down   , cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/centralWeight; 
        }
      }
      if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
        // WH boosted weights
        if(selection==kWHLightFlavorFJCR) {
          bLoad(b["sf_sjbtag0"],ientry);
          bLoad(b["sf_sjbtag0BUp"],ientry);
          bLoad(b["sf_sjbtag0BDown"],ientry);
          bLoad(b["sf_sjbtag0MUp"],ientry);
          bLoad(b["sf_sjbtag0MDown"],ientry);
          weight_btagBUp   = weight * (*sf_sjbtag0BUp  ) / (*sf_sjbtag0);
          weight_btagBDown = weight * (*sf_sjbtag0BDown) / (*sf_sjbtag0);
          weight_btagMUp   = weight * (*sf_sjbtag0MUp  ) / (*sf_sjbtag0);
          weight_btagMDown = weight * (*sf_sjbtag0MDown) / (*sf_sjbtag0);
        } else {
          //bLoad(b["sf_sjbtag2"],ientry);
          //bLoad(b["sf_sjbtag2BUp"],ientry);
          //bLoad(b["sf_sjbtag2BDown"],ientry);
          //bLoad(b["sf_sjbtag2MUp"],ientry);
          //bLoad(b["sf_sjbtag2MDown"],ientry);
          //weight_btagBUp   = weight * (*sf_sjbtag2BUp  ) / (*sf_sjbtag2);
          //weight_btagBDown = weight * (*sf_sjbtag2BDown) / (*sf_sjbtag2);
          //weight_btagMUp   = weight * (*sf_sjbtag2MUp  ) / (*sf_sjbtag2);
          //weight_btagMDown = weight * (*sf_sjbtag2MDown) / (*sf_sjbtag2);
          bLoad(b["sf_sjbtagGT0"],ientry);
          bLoad(b["sf_sjbtagGT0BUp"],ientry);
          bLoad(b["sf_sjbtagGT0BDown"],ientry);
          bLoad(b["sf_sjbtagGT0MUp"],ientry);
          bLoad(b["sf_sjbtagGT0MDown"],ientry);
          weight_btagBUp   = weight * (*sf_sjbtagGT0BUp  ) / (*sf_sjbtagGT0);
          weight_btagBDown = weight * (*sf_sjbtagGT0BDown) / (*sf_sjbtagGT0);
          weight_btagMUp   = weight * (*sf_sjbtagGT0MUp  ) / (*sf_sjbtagGT0);
          weight_btagMDown = weight * (*sf_sjbtagGT0MDown) / (*sf_sjbtagGT0);
        }
      }

    }
    // Category Assignment
    switch(sample) {
      case kData:
        theCategory=kPlotData; break;
      case kQCD:
        theCategory=kPlotQCD;  break;
      case kVZ:
        //bLoad(b["nB"],ientry);
        //if(gt.nB>=2) theCategory=kPlotVZbb;
        //else      theCategory=kPlotVVLF;
        bLoad(b["nBGenJets"],ientry);
        if(gt.nBGenJets>=2) theCategory=kPlotVZbb;
        else theCategory=kPlotVVLF;
        break;
      case kWW:
        theCategory=kPlotVVLF; break;
      case kTT:
        theCategory=kPlotTT; break;
      case kTop:
        theCategory=kPlotTop; break;
      case kWjets:
        //bLoad(b["nB"],ientry);
        //if(gt.nB>=2) theCategory=kPlotWbb;
        //else if(gt.nB==1) theCategory=kPlotWb;
        //else theCategory=kPlotWLF; // light flavor
        //if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
        //} else {
          bLoad(b["nBGenJets"],ientry);
          if(gt.nBGenJets>=2) theCategory=kPlotWbb;
          else if(gt.nBGenJets==1) theCategory=kPlotWb;
          else theCategory=kPlotWLF; // light flavor
          break;
        //}
      case kZjets:
        //if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
        //} else {
          bLoad(b["nBGenJets"],ientry);
          if(gt.nBGenJets>=2) theCategory=kPlotZbb;
          else if(gt.nBGenJets==1) theCategory=kPlotZb;
          else theCategory=kPlotZLF; // light flavor
          break;
        //}
      case kVH:
        theCategory=kPlotVH; break;
      default:
        theCategory=nPlotCategories;
        break;
    }
    plotTree->Fill();
    //if(ientry%100000==0) printf("\tread %d bytes\n",nBytesRead);
  }
  plotTree->Write();
  outputFile->Close();
  inputFile->Close();
  printf("Done\n");
  return true;




}

