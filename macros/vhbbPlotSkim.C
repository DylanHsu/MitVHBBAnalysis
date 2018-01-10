#include "PandaAnalysis/Flat/interface/GeneralTree.h"
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

using namespace vhbbPlot;
bool vhbbPlotSkim(
  TString inputFileName="", 
  TString outputFileName="",
  vhbbPlot::sampleType sample=nSampleTypes,
  vhbbPlot::selectionType selection = kWHLightFlavorCR,
  bool debug=false,
  Long64_t maxEntries=0
) {
  const double theLumi=35900.;
  
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
   
  // temporary
  TFile *muTrigEffFile = TFile::Open("/home/dhsu/CMSSW_8_0_29/src/PandaAnalysis/data/trigger_eff/muon_trig_Run2016BtoF.root", "READ");
  if(!muTrigEffFile && muTrigEffFile->IsOpen()) { throw std::runtime_error(""); return false; }
  TH2D *muTrigEff = (TH2D*)muTrigEffFile->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");

  TFile *outputFile = TFile::Open(outputFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  TTree *plotTree = new TTree("plotTree","Tree for making plots");
  unsigned char typeLepSel, theCategory;
  unsigned selectionBits, selectionBits_jesUp, selectionBits_jesDown, nMinusOneBits;
  float lepton1Pt, lepton1Eta, lepton1Phi, lepton1RelIso;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1Eta, hbbJet1Phi;
  float hbbJet2Pt, hbbJet2Eta, hbbJet2Phi;
  float hbbJet1PtUp, hbbJet1PtDown, hbbJet2PtUp, hbbJet2PtDown;
  float hbbDijetPt, hbbDijetPtUp, hbbDijetPtDown;
  float hbbDijetMass, hbbDijetMassUp, hbbDijetMassDown;
  float bDiscrMin, bDiscrMax;
  float deltaPhiLep1Met, deltaPhiVH;
  float topWBosonPt, topWBosonPt_jesUp, topWBosonPt_jesDown, mT_jesUp, mT_jesDown; // need to be in PandaExpress ntuples
  
  // For boosted categories, the number of 30 GeV AK4 jets not in the fat jet
  int nIsojet=0, nIsojet_jesUp=0, nIsojet_jesDown=0; 
  vector<unsigned char> isojets, isojets_jesUp, isojets_jesDown;

  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_scaleUp, weight_scaleDown;
  float weight_lepSFUp;
  float weight_cmvaUp, weight_cmvaDown; 
  
  float weight_btag0BUp, weight_btag0BDown;
  float weight_btag0MUp, weight_btag0MDown;
  float *sf_btag0      = (float*)dummyTree->GetBranch("sf_btag0"     )->GetAddress();
  float *sf_btag0BUp   = (float*)dummyTree->GetBranch("sf_btag0BUp"  )->GetAddress();
  float *sf_btag0BDown = (float*)dummyTree->GetBranch("sf_btag0BDown")->GetAddress();
  float *sf_btag0MUp   = (float*)dummyTree->GetBranch("sf_btag0MUp"  )->GetAddress();
  float *sf_btag0MDown = (float*)dummyTree->GetBranch("sf_btag0MDown")->GetAddress();
  
  plotTree->Branch("selectionBits"   , &selectionBits   );
  plotTree->Branch("selectionBits_jesUp"   , &selectionBits_jesUp   );
  plotTree->Branch("selectionBits_jesDown" , &selectionBits_jesDown );
  plotTree->Branch("nMinusOneBits"   , &nMinusOneBits   );
  plotTree->Branch("theCategory"     , &theCategory     );
  plotTree->Branch("pfmet"           , &gt.pfmet        );
  plotTree->Branch("pfmetUp"         , &gt.pfmetUp      );
  plotTree->Branch("pfmetDown"       , &gt.pfmetDown    );
  plotTree->Branch("pfmetsig"        , &gt.pfmetsig     );
  plotTree->Branch("pfmetphi"        , &gt.pfmetphi     );
  plotTree->Branch("weight"          , &weight          );
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    plotTree->Branch("nJet"            , &gt.nJet         );
    plotTree->Branch("hbbDijetPt"      , &hbbDijetPt      );
    plotTree->Branch("hbbDijetPtUp"    , &hbbDijetPtUp    );
    plotTree->Branch("hbbDijetPtDown"  , &hbbDijetPtDown  );
    plotTree->Branch("hbbDijetPt"      , &hbbDijetPt      );
    plotTree->Branch("hbbDijetMass"    , &hbbDijetMass    );
    plotTree->Branch("hbbDijetMassUp"  , &hbbDijetMassUp  );
    plotTree->Branch("hbbDijetMassDown", &hbbDijetMassDown);
    plotTree->Branch("bDiscrMin"       , &bDiscrMin       );
    plotTree->Branch("bDiscrMax"       , &bDiscrMax       );
    plotTree->Branch("hbbJet1Pt"       , &hbbJet1Pt       );
    plotTree->Branch("hbbJet1Eta"      , &hbbJet1Eta      );
    plotTree->Branch("hbbJet1Phi"      , &hbbJet1Phi      );
    plotTree->Branch("hbbJet2Pt"       , &hbbJet2Pt       );
    plotTree->Branch("hbbJet2Eta"      , &hbbJet2Eta      );
    plotTree->Branch("hbbJet2Phi"      , &hbbJet2Phi      );
    plotTree->Branch("hbbJet1PtUp"     , &hbbJet1PtUp     );
    plotTree->Branch("hbbJet1PtDown"   , &hbbJet1PtDown   );
    plotTree->Branch("hbbJet2PtUp"     , &hbbJet2PtUp     );
    plotTree->Branch("hbbJet2PtDown"   , &hbbJet2PtDown   );
    plotTree->Branch("typeLepSel"            , &typeLepSel               );
    plotTree->Branch("lepton1Pt"             , &lepton1Pt                );
    plotTree->Branch("lepton1Eta"            , &lepton1Eta               );
    plotTree->Branch("lepton1Phi"            , &lepton1Phi               );
    plotTree->Branch("lepton1RelIso"         , &lepton1RelIso            );
    plotTree->Branch("topMassLep1Met"        , &gt.topMassLep1Met        );
    plotTree->Branch("topMassLep1Met_jesUp"  , &gt.topMassLep1Met_jesUp  );
    plotTree->Branch("topMassLep1Met_jesDown", &gt.topMassLep1Met_jesDown);
    plotTree->Branch("topWBosonPt"           , &gt.topWBosonPt           );
    plotTree->Branch("topWBosonPt_jesUp"     , &topWBosonPt_jesUp        );
    plotTree->Branch("topWBosonPt_jesDown"   , &topWBosonPt_jesDown      );
    plotTree->Branch("topWBosonEta"          , &gt.topWBosonEta          );
    plotTree->Branch("topWBosonPhi"          , &gt.topWBosonPhi          );
    plotTree->Branch("mT"                    , &gt.mT                    );
    plotTree->Branch("mT_jesUp"              , &mT_jesUp                 );
    plotTree->Branch("mT_jesDown"            , &mT_jesDown               );
    plotTree->Branch("sumEtSoft1"            , &gt.sumEtSoft1            );
    plotTree->Branch("nSoft2"                , &gt.nSoft2                );
    plotTree->Branch("nSoft5"                , &gt.nSoft5                );
    plotTree->Branch("nSoft10"               , &gt.nSoft10               );
    plotTree->Branch("topWBosonCosThetaCS"   , &gt.topWBosonCosThetaCS   );
    plotTree->Branch("hbbCosThetaJJ"         , &gt.hbbCosThetaJJ         );
    plotTree->Branch("hbbCosThetaCSJ1"       , &gt.hbbCosThetaCSJ1       );
    plotTree->Branch("deltaPhiLep1Met"       , &deltaPhiLep1Met          );
    plotTree->Branch("deltaPhiVH"            , &deltaPhiVH               );
    plotTree->Branch("weight_pdfUp"          , &weight_pdfUp             );
    plotTree->Branch("weight_pdfDown"        , &weight_pdfDown           );
    plotTree->Branch("weight_scaleUp"        , &weight_scaleUp           );
    plotTree->Branch("weight_scaleDown"      , &weight_scaleDown         );
    plotTree->Branch("weight_lepSFUp"        , &weight_lepSFUp           );
    plotTree->Branch("weight_cmvaUp"         , &weight_cmvaUp            );
    plotTree->Branch("weight_cmvaDown"       , &weight_cmvaDown          );
  } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
    plotTree->Branch("topWBosonPt"           , &gt.topWBosonPt           );
    plotTree->Branch("topWBosonPt_jesUp"     , &topWBosonPt_jesUp        );
    plotTree->Branch("topWBosonPt_jesDown"   , &topWBosonPt_jesDown      );
    plotTree->Branch("topWBosonPhi"          , &gt.topWBosonPhi          );
    plotTree->Branch("deltaPhiLep1Met"       , &deltaPhiLep1Met          );
    plotTree->Branch("deltaPhiVH"            , &deltaPhiVH               );
    plotTree->Branch("typeLepSel"            , &typeLepSel               );
    plotTree->Branch("lepton1Pt"             , &lepton1Pt                );
    plotTree->Branch("lepton1Eta"            , &lepton1Eta               );
    plotTree->Branch("lepton1Phi"            , &lepton1Phi               );
    plotTree->Branch("lepton1RelIso"         , &lepton1RelIso            );
    plotTree->Branch("nIsojet"               , &nIsojet                  );
    plotTree->Branch("mT"                    , &gt.mT                    );
    plotTree->Branch("mT_jesUp"              , &mT_jesUp                 );
    plotTree->Branch("mT_jesDown"            , &mT_jesDown               );
    plotTree->Branch("fj1Tau32"              , &gt.fj1Tau32              );
    plotTree->Branch("fj1Tau21"              , &gt.fj1Tau21              );
    plotTree->Branch("fj1Tau32SD"            , &gt.fj1Tau32SD            ); 
    plotTree->Branch("fj1Tau21SD"            , &gt.fj1Tau21SD            ); 
    plotTree->Branch("fj1MSD"                , &gt.fj1MSD                );
    plotTree->Branch("fj1MSDScaleUp"         , &gt.fj1MSDScaleUp         );    
    plotTree->Branch("fj1MSDScaleDown"       , &gt.fj1MSDScaleDown       );      
    plotTree->Branch("fj1MSDSmeared"         , &gt.fj1MSDSmeared         );    
    plotTree->Branch("fj1MSDSmearedUp"       , &gt.fj1MSDSmearedUp       );      
    plotTree->Branch("fj1MSDSmearedDown"     , &gt.fj1MSDSmearedDown     );        
    plotTree->Branch("fj1Pt"                 , &gt.fj1Pt                 );
    plotTree->Branch("fj1PtScaleUp"          , &gt.fj1PtScaleUp          );   
    plotTree->Branch("fj1PtScaleDown"        , &gt.fj1PtScaleDown        );     
    plotTree->Branch("fj1PtSmeared"          , &gt.fj1PtSmeared          );   
    plotTree->Branch("fj1PtSmearedUp"        , &gt.fj1PtSmearedUp        );     
    plotTree->Branch("fj1PtSmearedDown"      , &gt.fj1PtSmearedDown      );       
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
    plotTree->Branch("weight_pdfUp"          , &weight_pdfUp             );
    plotTree->Branch("weight_pdfDown"        , &weight_pdfDown           );
    plotTree->Branch("weight_scaleUp"        , &weight_scaleUp           );
    plotTree->Branch("weight_scaleDown"      , &weight_scaleDown         );
    plotTree->Branch("weight_lepSFUp"        , &weight_lepSFUp           );
    plotTree->Branch("weight_cmvaUp"         , &weight_cmvaUp            );
    plotTree->Branch("weight_cmvaDown"       , &weight_cmvaDown          );
  }
  

  delete dummyTree; // no longer needed
  Long64_t nentries = events->GetEntries();
  if(maxEntries!=0) nentries = TMath::Min((ULong64_t)maxEntries, (ULong64_t)nentries);
  for (Long64_t ientry=0; ientry<nentries; ientry++) {
    if(debug) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries);
    else if(ientry%100000==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries);

    int nBytesRead=0;
    theCategory=-1; // plot category 
    typeLepSel=99; // 0: mixed e-mu, 1: all mu, 2: all e, 99: undefined
    
    // Vectors to keep around in memory
    TVector2 vectorBosonV2;
    TLorentzVector hbbJet1P4, hbbJet2P4, hbbDijetP4;
    
    // Analysis Preselection
    if(selection>=kWHLightFlavorCR && selection<=kWHSR) { // WH Resolved Category
      {
        
        // Lepton multiplicity
        nBytesRead+=bLoad(b["nLooseLep"],ientry);
        if(gt.nLooseLep<1) continue; //N_al = 0
        if(debug) printf("Passed lepton multiplicity\n");
        
        // Trigger
        if(sample==kData) {
          nBytesRead+=bLoad(b["trigger"],ientry);
          if( (gt.trigger & (1<<3 | 1<<1))==0 ) continue;
        }

        nBytesRead+=bLoad(b["metFilter"],ientry);
        if(gt.metFilter!=1) continue;
        if(debug) printf("Passed MET filters\n");
     
        // Jet multiplicity
        nBytesRead+=bLoad(b["nJet"],ientry);
        nBytesRead+=bLoad(b["nFatjet"],ientry);
        if     (gt.nJet<2) continue;
        if     (gt.nFatjet!=0) continue;
        // Jet kinematics
        nBytesRead+=bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
        nBytesRead+=bLoad(b["jetPt"],ientry);
        if(debug) printf("hbb jet1 pt %.2f, jet2 pt %.2f\n", gt.jetPt[gt.hbbjtidx[0]], gt.jetPt[gt.hbbjtidx[1]]);
        if(gt.jetPt[gt.hbbjtidx[0]]<25 || gt.jetPt[gt.hbbjtidx[1]]<25) continue;
        
        nBytesRead+=bLoad(b["hbbpt_reg"],ientry);
        nBytesRead+=bLoad(b["hbbeta"],ientry);
        nBytesRead+=bLoad(b["hbbphi"],ientry);
        nBytesRead+=bLoad(b["hbbm_reg"],ientry);
        nBytesRead+=bLoad(b["jetEta"],ientry);
        nBytesRead+=bLoad(b["jetPhi"],ientry);
        hbbDijetPt = gt.hbbpt_reg;
        hbbDijetMass = gt.hbbm_reg;
        hbbJet1Pt=gt.jetPt[gt.hbbjtidx[0]]; hbbJet1Eta=gt.jetEta[gt.hbbjtidx[0]]; hbbJet1Phi=gt.jetPhi[gt.hbbjtidx[0]];
        hbbJet2Pt=gt.jetPt[gt.hbbjtidx[1]]; hbbJet2Eta=gt.jetEta[gt.hbbjtidx[1]]; hbbJet2Phi=gt.jetPhi[gt.hbbjtidx[1]];
        if(debug) printf("hbbDijetPt %.2f, hbbDijetMass %.2f\n", hbbDijetPt, hbbDijetMass); 
        if(hbbDijetPt<50.) continue;
        if(hbbDijetMass<0 || hbbDijetMass>250.) continue;
        if(debug) printf("passed jet kinematics\n");

        // Lepton ID and isolation
        nBytesRead+=bLoad(b["nTightMuon"],ientry);
        if(gt.nTightMuon==1 && (gt.trigger & 1<<3)!=0) typeLepSel=1; else {
          nBytesRead+=bLoad(b["nTightElectron"],ientry);
          if(gt.nTightElectron==1 && (gt.trigger & 1<<1)!=0) typeLepSel=2;
        } if(typeLepSel!=1 && typeLepSel!=2) continue;
        if(debug) printf("Passed lepton ID/iso multiplicity\n");

        // Lepton kinematics
        if     (typeLepSel==1) {
          nBytesRead+=bLoad(b["muonPt"],ientry);
          nBytesRead+=bLoad(b["muonEta"],ientry);
          nBytesRead+=bLoad(b["muonPhi"],ientry);
          nBytesRead+=bLoad(b["muonCombIso"],ientry);
          lepton1Pt = gt.muonPt[0]; lepton1Eta = gt.muonEta[0]; lepton1Phi = gt.muonPhi[0]; lepton1RelIso = gt.muonCombIso[0]/gt.muonPt[0];
        } else if(typeLepSel==2) {
          nBytesRead+=bLoad(b["electronPt"],ientry);
          nBytesRead+=bLoad(b["electronEta"],ientry);
          nBytesRead+=bLoad(b["electronPhi"],ientry);
          lepton1Pt = gt.electronPt[0]; lepton1Eta = gt.electronEta[0]; lepton1Phi = gt.electronPhi[0]; lepton1RelIso = gt.electronCombIso[0]/gt.electronPt[0];
        } else continue;
        if(debug) printf("Passed lepton kinematics\n");
        
        // W reconstruction
        nBytesRead+=bLoad(b["topWBosonPt"],ientry);
        if(gt.topWBosonPt<50.) continue;
        if(debug) printf("passed W reco\n");
      } // end preselection
      
      // Met
      nBytesRead+=bLoad(b["pfmet"],ientry);
      nBytesRead+=bLoad(b["pfmetUp"],ientry);
      nBytesRead+=bLoad(b["pfmetDown"],ientry);
      nBytesRead+=bLoad(b["pfmetphi"],ientry);
      nBytesRead+=bLoad(b["pfmetsig"],ientry);
      
      // Jet B-tagging
      nBytesRead+=bLoad(b["jetCMVA"],ientry);
      bDiscrMax=gt.jetCMVA[gt.hbbjtidx[0]];
      bDiscrMin=gt.jetCMVA[gt.hbbjtidx[1]];
      
      // Top reconstruction
      nBytesRead+=bLoad(b["topMassLep1Met"],ientry);
      nBytesRead+=bLoad(b["topMassLep1Met_jesUp"],ientry);
      nBytesRead+=bLoad(b["topMassLep1Met_jesDown"],ientry);
      nBytesRead+=bLoad(b["topWBosonPt"],ientry);
      nBytesRead+=bLoad(b["topWBosonEta"],ientry);
      nBytesRead+=bLoad(b["topWBosonPhi"],ientry);
      nBytesRead+=bLoad(b["mT"],ientry);
      
      // mT, W pT up down -- temporary
      { 
        TVector2 metV2Up, metV2Down, lepton1V2;
        metV2Up.SetMagPhi(gt.pfmetUp, gt.pfmetphi);
        metV2Down.SetMagPhi(gt.pfmetDown, gt.pfmetphi);
        lepton1V2.SetMagPhi(lepton1Pt, lepton1Phi);
        mT_jesUp   = TMath::Sqrt(2.*gt.pfmetUp*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        mT_jesDown = TMath::Sqrt(2.*gt.pfmetDown*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        topWBosonPt_jesUp   = (metV2Up+lepton1V2).Mod();
        topWBosonPt_jesDown = (metV2Down+lepton1V2).Mod();
        //if(debug) printf("topWBosonPt %.1f, topWBosonPt_jesUp %.1f, topWBosonPt_jesDown %.1f\n", topWBosonPt, topWBosonPt_jesUp, topWBosonPt_jesDown);
      }
      // Other stuff
      nBytesRead+=bLoad(b["hbbphi"],ientry);
      deltaPhiVH = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.hbbphi));
      deltaPhiLep1Met = fabs(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi));
      
      nBytesRead+=bLoad(b["hbbpt_reg_jesUp"],ientry);
      nBytesRead+=bLoad(b["hbbpt_reg_jesDown"],ientry);
      nBytesRead+=bLoad(b["hbbm_reg_jesUp"],ientry);
      nBytesRead+=bLoad(b["hbbm_reg_jesDown"],ientry);
      hbbDijetPtUp = gt.hbbpt_reg_jesUp;
      hbbDijetPtDown = gt.hbbpt_reg_jesDown;
      hbbDijetMassUp = gt.hbbm_reg_jesUp;
      hbbDijetMassDown = gt.hbbm_reg_jesDown;
      nBytesRead+=bLoad(b["jetPtUp"],ientry);
      nBytesRead+=bLoad(b["jetPtDown"],ientry);
      hbbJet1PtUp = gt.jetPtUp[gt.hbbjtidx[0]];
      hbbJet1PtDown = gt.jetPtDown[gt.hbbjtidx[0]];
      hbbJet2PtUp = gt.jetPtUp[gt.hbbjtidx[1]];
      hbbJet2PtDown = gt.jetPtDown[gt.hbbjtidx[1]];
      nBytesRead+=bLoad(b["topWBosonCosThetaCS"],ientry);
      nBytesRead+=bLoad(b["hbbCosThetaJJ"],ientry);
      nBytesRead+=bLoad(b["hbbCosThetaCSJ1"],ientry);

      // Jet multiplicity for JES
      nBytesRead+=bLoad(b["nJet_jesUp"],ientry);
      nBytesRead+=bLoad(b["nJet_jesDown"],ientry);
      nBytesRead+=bLoad(b["sumEtSoft1"],ientry);
      nBytesRead+=bLoad(b["nSoft2"],ientry);
      nBytesRead+=bLoad(b["nSoft5"],ientry);
      nBytesRead+=bLoad(b["nSoft10"],ientry);
      
      // End WH Resolved Category 
    } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) { 
      // WH Boosted Category
      {
        // Lepton multiplicity
        nBytesRead+=bLoad(b["nLooseLep"],ientry);
        if(gt.nLooseLep<1) continue; //N_al = 0
        if(debug) printf("Passed lepton multiplicity\n");
        
        // Trigger
        if(sample==kData) {
          nBytesRead+=bLoad(b["trigger"],ientry);
          if( (gt.trigger & (1<<3 | 1<<1))==0 ) continue;
        }

        nBytesRead+=bLoad(b["metFilter"],ientry);
        if(gt.metFilter!=1) continue;
        if(debug) printf("Passed MET filters\n");
     
        // Jet multiplicity
        nBytesRead+=bLoad(b["nJet"],ientry);
        nBytesRead+=bLoad(b["nFatjet"],ientry);
        if     (gt.nFatjet==0) continue;
        // Jet kinematics
        nBytesRead+=bLoad(b["fj1Pt"],ientry);
        if(fabs(gt.fj1Eta)>2.4) continue;
        if(debug) printf("passed jet kinematics\n");

        // Lepton ID and isolation
        nBytesRead+=bLoad(b["nTightMuon"],ientry);
        if(gt.nTightMuon==1 && (gt.trigger & 1<<3)!=0) typeLepSel=1; else {
          nBytesRead+=bLoad(b["nTightElectron"],ientry);
          if(gt.nTightElectron==1 && (gt.trigger & 1<<1)!=0) typeLepSel=2;
        } if(typeLepSel!=1 && typeLepSel!=2) continue;
        if(debug) printf("Passed lepton ID/iso multiplicity\n");

        // Lepton kinematics
        if     (typeLepSel==1) {
          nBytesRead+=bLoad(b["muonPt"],ientry);
          nBytesRead+=bLoad(b["muonEta"],ientry);
          nBytesRead+=bLoad(b["muonPhi"],ientry);
          nBytesRead+=bLoad(b["muonCombIso"],ientry);
          lepton1Pt = gt.muonPt[0]; lepton1Eta = gt.muonEta[0]; lepton1Phi = gt.muonPhi[0]; lepton1RelIso = gt.muonCombIso[0]/gt.muonPt[0];
        } else if(typeLepSel==2) {
          nBytesRead+=bLoad(b["electronPt"],ientry);
          nBytesRead+=bLoad(b["electronEta"],ientry);
          nBytesRead+=bLoad(b["electronPhi"],ientry);
          lepton1Pt = gt.electronPt[0]; lepton1Eta = gt.electronEta[0]; lepton1Phi = gt.electronPhi[0]; lepton1RelIso = gt.electronCombIso[0]/gt.electronPt[0];
        } else continue;
        if(debug) printf("Passed lepton kinematics\n");
        
        // W reconstruction
        nBytesRead+=bLoad(b["topWBosonPt"],ientry);
        if(gt.topWBosonPt<50.) continue;
        if(debug) printf("passed W reco\n");
      } // end preselection
      
      // Met
      nBytesRead+=bLoad(b["pfmet"],ientry);
      nBytesRead+=bLoad(b["pfmetUp"],ientry);
      nBytesRead+=bLoad(b["pfmetDown"],ientry);
      nBytesRead+=bLoad(b["pfmetphi"],ientry);
      nBytesRead+=bLoad(b["pfmetsig"],ientry);
      
      // Fatjet Properties
      nBytesRead+=bLoad(b["fj1Tau32"],ientry);
      nBytesRead+=bLoad(b["fj1Tau21"],ientry);
      nBytesRead+=bLoad(b["fj1Tau32SD"],ientry); 
      nBytesRead+=bLoad(b["fj1Tau21SD"],ientry); 
      nBytesRead+=bLoad(b["fj1MSDScaleUp"],ientry);    
      nBytesRead+=bLoad(b["fj1MSDScaleDown"],ientry);      
      nBytesRead+=bLoad(b["fj1MSDSmeared"],ientry);    
      nBytesRead+=bLoad(b["fj1MSDSmearedUp"],ientry);      
      nBytesRead+=bLoad(b["fj1MSDSmearedDown"],ientry);        
      nBytesRead+=bLoad(b["fj1PtScaleUp"],ientry);   
      nBytesRead+=bLoad(b["fj1PtScaleDown"],ientry);     
      nBytesRead+=bLoad(b["fj1PtSmeared"],ientry);   
      nBytesRead+=bLoad(b["fj1PtSmearedUp"],ientry);     
      nBytesRead+=bLoad(b["fj1PtSmearedDown"],ientry);       
      nBytesRead+=bLoad(b["fj1Phi"],ientry);
      nBytesRead+=bLoad(b["fj1MSD"],ientry);
      nBytesRead+=bLoad(b["fj1Eta"],ientry);
      nBytesRead+=bLoad(b["fj1M"],ientry);
      nBytesRead+=bLoad(b["fj1MaxCSV"],ientry);
      nBytesRead+=bLoad(b["fj1MinCSV"],ientry);
      nBytesRead+=bLoad(b["fj1DoubleCSV"],ientry);   
      nBytesRead+=bLoad(b["fj1HTTMass"],ientry); 
      nBytesRead+=bLoad(b["fj1HTTFRec"],ientry); 
      nBytesRead+=bLoad(b["fj1SDEFrac100"],ientry);    
       
      // Ak4 jet properties
      nBytesRead+=bLoad(b["jetPt"],ientry);
      nBytesRead+=bLoad(b["jetEta"],ientry);
      nBytesRead+=bLoad(b["jetPhi"],ientry);
      nBytesRead+=bLoad(b["jetPtUp"],ientry);
      nBytesRead+=bLoad(b["jetPtDown"],ientry);
      nBytesRead+=bLoad(b["jetCMVA"],ientry);
      //nBytesRead+=bLoad(b["jetIso"],ientry); // broken
      //nBytesRead+=bLoad(b["isojetNBtags"],ientry); //broken

      isojets.clear();
      isojets_jesUp.clear();
      isojets_jesDown.clear();
      unsigned char isojetNBtags=0;
      for(unsigned char iJ=0; iJ<gt.nJet; iJ++) {
        // Cannot use jetIso for counting central jets right now
        //if(!gt.jetIso[iJ]) continue; 
        if(fabs(gt.jetEta[iJ])>2.4) continue;
        if(get.jetPt[iJ]>30 && gt.jetCMVA[iJ]>bDiscrLoose) isojetNBtags++;
        float dRJetFatjet=sqrt(pow(gt.jetEta[iJ]-gt.fj1Eta,2)+pow(TVector2::Phi_mpi_pi(gt.jetPhi[iJ]-gt.fj1Phi),2));
        if(dRJetFatjet<2.25) continue;
        if(gt.jetPt[iJ]>30) isojets.push_back(iJ);
        if(gt.jetPtUp[iJ]>30) isojets_jesUp.push_back(iJ);
        if(gt.jetPtDown[iJ]>30) isojets_jesDown.push_back(iJ);
      }
      nIsojet=isojets.size();
      nIsojet_jesUp=isojets_jesUp.size();
      nIsojet_jesDown=isojets_jesDown.size();
 
      // Energy Correlation Functions
      nBytesRead+=bLoad(b["fj1ECFN_1_1_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_1_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_1_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_2_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_2_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_2_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_3_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_3_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_3_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_4_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_4_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_4_05"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_1_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_1_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_1_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_2_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_2_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_2_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_3_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_3_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_3_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_4_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_4_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_4_10"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_1_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_1_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_1_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_2_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_2_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_2_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_3_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_3_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_3_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_4_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_4_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_4_20"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_1_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_1_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_1_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_2_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_2_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_2_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_3_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_3_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_3_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_1_4_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_2_4_40"],ientry);
      nBytesRead+=bLoad(b["fj1ECFN_3_4_40"],ientry);

      // W reconstruction
      nBytesRead+=bLoad(b["topWBosonPt"],ientry);
      nBytesRead+=bLoad(b["topWBosonPhi"],ientry);
      nBytesRead+=bLoad(b["mT"],ientry);
      
      // mT, W pT up down
      { 
        TVector2 metV2Up, metV2Down, lepton1V2;
        metV2Up.SetMagPhi(gt.pfmetUp, gt.pfmetphi);
        metV2Down.SetMagPhi(gt.pfmetDown, gt.pfmetphi);
        lepton1V2.SetMagPhi(lepton1Pt, lepton1Phi);
        mT_jesUp   = TMath::Sqrt(2.*gt.pfmetUp*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        mT_jesDown = TMath::Sqrt(2.*gt.pfmetDown*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi))));
        topWBosonPt_jesUp   = (metV2Up+lepton1V2).Mod();
        topWBosonPt_jesDown = (metV2Down+lepton1V2).Mod();
        //if(debug) printf("topWBosonPt %.1f, topWBosonPt_jesUp %.1f, topWBosonPt_jesDown %.1f\n", topWBosonPt, topWBosonPt_jesUp, topWBosonPt_jesDown);
      }
      
      // Other stuff
      deltaPhiVH = fabs(TVector2::Phi_mpi_pi(gt.topWBosonPhi-gt.fj1Phi));
      deltaPhiLep1Met = fabs(TVector2::Phi_mpi_pi(lepton1Phi-gt.pfmetphi));
      // End WH Boosted Category
    }
 
    // Set Selection Bits
    selectionBits=0; nMinusOneBits=0; selectionBits_jesUp=0; selectionBits_jesDown=0;
    std::map<TString, bool> cut, cut_jesUp, cut_jesDown;
    
    if(selection>=kWHLightFlavorCR && selection<=kWHSR) { // Begin WH Resolved Selection

      cut["2ndLepVeto" ] = gt.nLooseLep==1;
      cut["ultraLepIso"] = lepton1RelIso<0.06;
      cut["WpT"        ] = gt.topWBosonPt>100;
      cut["pTjj"       ] = hbbDijetPt>100; 
      cut["lepton1Pt"  ] = ((typeLepSel==1 && lepton1Pt>25) || (typeLepSel==2 && lepton1Pt>30));
      cut["dPhiVH"     ] = deltaPhiVH > 2.5;
      cut["dPhiLep1Met"] = deltaPhiLep1Met < 2;
      cut["tightBTag"  ] = (bDiscrMax >= bDiscrTight);
      cut["tightBVeto" ] = (bDiscrMax < bDiscrTight);
      cut["looseBTag"  ] = (bDiscrMax >= bDiscrLoose);
      cut["looseBTag2" ] = (bDiscrMin >= bDiscrLoose);
      cut["metSig"     ] = gt.pfmetsig>2;
      cut["nJet_WHTT"  ] = gt.nJet>=4;
      cut["nJet_WHHF"  ] = gt.nJet==2;
      cut["nJet_WHSR"  ] = gt.nJet<4;
      cut["mH_WHSR"    ] = ((hbbDijetMass >= 90) && (hbbDijetMass <  150));
      cut["mH_Flip"    ] = ((hbbDijetMass < 90) || (hbbDijetMass >= 150));

      vector<TString> cutsWHLightFlavorCR, cutsWHHeavyFlavorCR, cutsWH2TopCR, cutsWHSR;
      cutsWHLightFlavorCR ={"ultraLepIso", "2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiLep1Met","looseBTag","tightBVeto","metSig"};
      cutsWHHeavyFlavorCR ={"ultraLepIso", "2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiLep1Met","nJet_WHHF","tightBTag","mH_Flip","metSig"};
      cutsWH2TopCR        ={"ultraLepIso", "2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiLep1Met","nJet_WHTT","tightBTag"};
      cutsWHSR            ={"ultraLepIso", "2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiLep1Met","nJet_WHSR","tightBTag","looseBTag2","mH_WHSR"};

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
      cut_jesUp["nJet_WHTT"  ] = gt.nJet_jesUp>=4;
      cut_jesUp["nJet_WHHF"  ] = gt.nJet_jesUp==2;
      cut_jesUp["nJet_WHSR"  ] = gt.nJet_jesUp<4;
      cut_jesUp["mH_WHSR"    ] = ((hbbDijetMassUp >= 90) && (hbbDijetMassUp <  150));
      cut_jesUp["mH_Flip"    ] = ((hbbDijetMassUp < 90) || (hbbDijetMassUp >= 150 && hbbDijetMassUp < 250));
      if(passAllCuts( cut_jesUp, cutsWHLightFlavorCR)) selectionBits_jesUp |= kWHLightFlavorCR;
      if(passAllCuts( cut_jesUp, cutsWHHeavyFlavorCR)) selectionBits_jesUp |= kWHHeavyFlavorCR;
      if(passAllCuts( cut_jesUp, cutsWH2TopCR       )) selectionBits_jesUp |= kWH2TopCR;
      if(passAllCuts( cut_jesUp, cutsWHSR           )) selectionBits_jesUp |= kWHSR;
      
      cut_jesDown=cut;
      cut_jesDown["WpT"        ] = topWBosonPt_jesDown>100;
      cut_jesDown["pTjj"       ] = hbbDijetPtDown>100; 
      cut_jesDown["nJet_WHTT"  ] = gt.nJet_jesDown>=4;
      cut_jesDown["nJet_WHHF"  ] = gt.nJet_jesDown==2;
      cut_jesDown["nJet_WHSR"  ] = gt.nJet_jesDown<4;
      cut_jesDown["mH_WHSR"    ] = ((hbbDijetMassDown >= 90) && (hbbDijetMassDown <  150));
      cut_jesDown["mH_Flip"    ] = ((hbbDijetMassDown < 90) || (hbbDijetMassDown >= 150 && hbbDijetMassDown < 250));
      if(passAllCuts( cut_jesDown, cutsWHLightFlavorCR)) selectionBits_jesDown |= kWHLightFlavorCR;
      if(passAllCuts( cut_jesDown, cutsWHHeavyFlavorCR)) selectionBits_jesDown |= kWHHeavyFlavorCR;
      if(passAllCuts( cut_jesDown, cutsWH2TopCR       )) selectionBits_jesDown |= kWH2TopCR;
      if(passAllCuts( cut_jesDown, cutsWHSR           )) selectionBits_jesDown |= kWHSR;
      // End WH Resolved Selection

    } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
      // Begin WH Boosted Selection
      float MSDmin=90, MSDmax=150;
      cut["2ndLepVeto" ] = gt.nLooseLep==1;
      cut["ultraLepIso"] = lepton1RelIso<0.06;
      cut["WpT"        ] = gt.topWBosonPt>100;
      cut["pTfj"       ] = gt.fj1Pt>220; 
      cut["lepton1Pt"  ] = ((typeLepSel==1 && lepton1Pt>25) || (typeLepSel==2 && lepton1Pt>30));
      cut["dPhiVH"     ] = deltaPhiVH > 2.5;
      cut["dPhiLep1Met"] = deltaPhiLep1Met < 2;
      cut["nJet_WHFJTT"] = nIsojet>=2;
      cut["nJet_WHFJHF"] = nIsojet==0;
      cut["nJet_WHFJSR"] = nIsojet<2;
      cut["mH_WHFJSR"  ] = ((gt.fj1MSD>=MSDmin && gt.fj1MSD<MSDmax));
      cut["mH_FJFlip"  ] = ((gt.fj1MSD<MSDmin || gt.fj1MSD>=MSDmax));
      cut["DoubleB"    ] = gt.fj1DoubleCSV >= 0.8;
      cut["DoubleBVeto"] = gt.fj1DoubleCSV < 0.8;
      cut["isojet0Btag"] = isojetNBtags==0;
      cut["metSig"     ] = gt.pfmetsig>2;
      
      vector<TString> cutsWHLightFlavorFJCR, cutsWHHeavyFlavorFJCR, cutsWH2TopFJCR, cutsWHFJSR;
      cutsWHLightFlavorFJCR ={"ultraLepIso", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","dPhiLep1Met",              "isojet0Btag","DoubleBVeto",               "metSig"};
      cutsWHHeavyFlavorFJCR ={"ultraLepIso", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","dPhiLep1Met","nJet_WHFJHF","isojet0Btag","DoubleB"    ,"mH_FJFlip"   ,"metSig"};
      cutsWH2TopFJCR        ={"ultraLepIso", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","dPhiLep1Met","nJet_WHFJTT","isojet0Btag","DoubleB"                            };
      cutsWHFJSR            ={"ultraLepIso", "2ndLepVeto","WpT","pTfj","lepton1Pt","dPhiVH","dPhiLep1Met","nJet_WHFJSR","isojet0Btag","DoubleB"    ,"mH_WHFJSR"            };

      if(passAllCuts( cut, cutsWHLightFlavorFJCR))   selectionBits |= kWHLightFlavorFJCR;
      if(passAllCuts( cut, cutsWHHeavyFlavorFJCR))   selectionBits |= kWHHeavyFlavorFJCR;
      if(passAllCuts( cut, cutsWH2TopFJCR       ))   selectionBits |= kWH2TopFJCR;
      if(passAllCuts( cut, cutsWHFJSR           ))   selectionBits |= kWHFJSR;
      if(passNMinusOne( cut, cutsWHLightFlavorFJCR)) nMinusOneBits |= kWHLightFlavorFJCR;
      if(passNMinusOne( cut, cutsWHHeavyFlavorFJCR)) nMinusOneBits |= kWHHeavyFlavorFJCR;
      if(passNMinusOne( cut, cutsWH2TopFJCR       )) nMinusOneBits |= kWH2TopFJCR;
      if(passNMinusOne( cut, cutsWHFJSR           )) nMinusOneBits |= kWHFJSR;

      cut_jesUp=cut;
      cut_jesUp["WpT"        ] = topWBosonPt_jesUp>100;
      cut_jesUp["pTfj"       ] = gt.fj1PtScaleUp>220; 
      cut_jesUp["nJet_WHFJTT"] = nIsojet_jesUp>=2;
      cut_jesUp["nJet_WHFJHF"] = nIsojet_jesUp==0;
      cut_jesUp["nJet_WHFJSR"] = nIsojet_jesUp<2;
      cut_jesUp["mH_WHFJSR"  ] = ((gt.fj1MSDScaleUp>=MSDmin && gt.fj1MSDScaleUp<MSDmax));
      cut_jesUp["mH_FJFlip"  ] = ((gt.fj1MSDScaleUp<MSDmin || gt.fj1MSDScaleUp>=MSDmax));
      if(passAllCuts( cut_jesUp, cutsWHLightFlavorFJCR)) selectionBits_jesUp |= kWHLightFlavorFJCR;
      if(passAllCuts( cut_jesUp, cutsWHHeavyFlavorFJCR)) selectionBits_jesUp |= kWHHeavyFlavorFJCR;
      if(passAllCuts( cut_jesUp, cutsWH2TopFJCR       )) selectionBits_jesUp |= kWH2TopFJCR;
      if(passAllCuts( cut_jesUp, cutsWHFJSR           )) selectionBits_jesUp |= kWHFJSR;
      
      cut_jesDown=cut;
      cut_jesDown["WpT"        ] = topWBosonPt_jesDown>100;
      cut_jesDown["pTfj"       ] = gt.fj1PtScaleDown>220; 
      cut_jesDown["nJet_WHFJTT"] = nIsojet_jesDown>=2;
      cut_jesDown["nJet_WHFJHF"] = nIsojet_jesDown==0;
      cut_jesDown["nJet_WHFJSR"] = nIsojet_jesDown<2;
      cut_jesDown["mH_WHFJSR"  ] = ((gt.fj1MSDScaleDown>=MSDmin && gt.fj1MSDScaleDown<MSDmax));
      cut_jesDown["mH_FJFlip"  ] = ((gt.fj1MSDScaleDown<MSDmin || gt.fj1MSDScaleDown>=MSDmax));
      if(passAllCuts( cut_jesDown, cutsWHLightFlavorFJCR)) selectionBits_jesDown |= kWHLightFlavorFJCR;
      if(passAllCuts( cut_jesDown, cutsWHHeavyFlavorFJCR)) selectionBits_jesDown |= kWHHeavyFlavorFJCR;
      if(passAllCuts( cut_jesDown, cutsWH2TopFJCR       )) selectionBits_jesDown |= kWH2TopFJCR;
      if(passAllCuts( cut_jesDown, cutsWHFJSR           )) selectionBits_jesDown |= kWHFJSR;
      // End WH Boosted Selection 
    }
    if(nMinusOneBits==0 && selectionBits==0 && selectionBits_jesUp==0 && selectionBits_jesDown==0) continue;

    // Weighting
    if(sample==kData) weight=1;
    else {
      nBytesRead+=bLoad(b["normalizedWeight"],ientry);
      nBytesRead+=bLoad(b["sf_npv"],ientry);
      if(sample==kWjets || sample==kZjets) { nBytesRead+=bLoad(b["sf_qcdV"],ientry); nBytesRead+=bLoad(b["sf_ewkV"],ientry); } 
      else { gt.sf_qcdV=1; gt.sf_ewkV=1; }
      //nBytesRead+=bLoad(sf_btag1,ientry);
      if(sample==kTT) nBytesRead+=bLoad(b["sf_tt"],ientry); else gt.sf_tt=1;
      weight = normalizedWeight * theLumi * gt.sf_npv * gt.sf_ewkV * gt.sf_qcdV * gt.sf_tt;
      if (typeLepSel==1) {
        nBytesRead+=bLoad(b["muonSfReco"],ientry);
        nBytesRead+=bLoad(b["muonSfTight"],ientry);
        nBytesRead+=bLoad(b["muonSfUnc"],ientry);
        nBytesRead+=bLoad(b["sf_muTrig"],ientry);
        weight *= gt.sf_muTrig * gt.muonSfReco[0] * gt.muonSfTight[0];
        //weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
        // need muon trigger efficiency
        //double theTrigEff = muTrigEff->GetBinContent( muTrigEff->FindBin( fabs(lepton1Eta), TMath::Max((float)26., TMath::Min((float)500., lepton1Pt)) ) );
        //weight *= theTrigEff;
      } else if(typeLepSel==2) {
        nBytesRead+=bLoad(b["electronSfReco"],ientry);
        nBytesRead+=bLoad(b["electronSfTight"],ientry);
        nBytesRead+=bLoad(b["electronSfUnc"],ientry);
        nBytesRead+=bLoad(b["sf_eleTrig"],ientry);
        weight *= gt.sf_eleTrig * gt.electronSfReco[0] * gt.electronSfTight[0];
      }
      if(sample==kVZ) {
        nBytesRead+=bLoad(b["sf_wz"],ientry);
        nBytesRead+=bLoad(b["sf_zz"],ientry);
        nBytesRead+=bLoad(b["sf_zzUnc"],ientry);
        weight *= gt.sf_wz * gt.sf_zz;
      } else if(sample==kVH) {
        nBytesRead+=bLoad(b["sf_zh"],ientry);
        nBytesRead+=bLoad(b["sf_zhUp"],ientry);
        nBytesRead+=bLoad(b["sf_zhDown"],ientry);
        weight *= gt.sf_zh;
      }
      // #############################
      // # Variations of the Weights #
      // #############################
      // PDF and QCD scale
      nBytesRead+=bLoad(b["pdfUp"],ientry);
      nBytesRead+=bLoad(b["pdfDown"],ientry);
      //nBytesRead+=bLoad(b["scaleUp"],ientry); // this is broken right now, thanks sid
      nBytesRead+=bLoad(b["scale"],ientry);
      nBytesRead+=bLoad(b["scaleDown"],ientry);
      weight_pdfUp = weight * gt.pdfUp;
      weight_pdfDown = weight * gt.pdfDown;
      weight_scaleUp = weight * (1.+std::max({gt.scale[0], gt.scale[1], gt.scale[2], gt.scale[3], gt.scale[4], gt.scale[5]}));
      weight_scaleDown = weight * (1.+gt.scaleDown);
      
      if(selection>=kWHLightFlavorCR && selection<=kWHFJSR) { // WH weights
        // Lepton ID & ISO SF uncertainty
        if(typeLepSel==1) weight_lepSFUp = weight * (1.+gt.muonSfUnc[0]);
        else weight_lepSFUp = weight * (1.+gt.electronSfUnc[0]);
      //}
      //if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
        // WH resolved weights

        // CMVA weight
        nBytesRead+=bLoad(b["sf_cmvaWeight_Cent"],ientry);
        weight *= gt.sf_csvWeights[GeneralTree::csvCent];
      
        // CMVA weight uncertainty
        // https://cmssdt.cern.ch/lxr/source/PhysicsTools/Heppy/python/physicsutils/BTagWeightCalculator.py?v=CMSSW_8_0_20
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors
        nBytesRead+=bLoad(b["sf_cmvaWeight_LFup"        ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_LFdown"      ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_HFup"        ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_HFdown"      ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_HFStats1up"  ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_HFStats1down"],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_HFStats2up"  ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_HFStats2down"],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_LFStats1up"  ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_LFStats1down"],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_LFStats2up"  ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_LFStats2down"],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_CErr1up"     ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_CErr1down"   ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_CErr2up"     ],ientry);
        nBytesRead+=bLoad(b["sf_cmvaWeight_CErr2down"   ],ientry);
        weight_cmvaUp = weight*(1.+sqrt(
          pow(gt.sf_csvWeights[GeneralTree::csvLFup       ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvHFup       ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvHFStats1up ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvHFStats2up ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvLFStats1up ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvLFStats2up ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvCErr1up    ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvCErr2up    ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2)
        ));
        weight_cmvaDown = weight*(1.-sqrt(
          pow(gt.sf_csvWeights[GeneralTree::csvLFdown       ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvHFdown       ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvHFStats1down ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvHFStats2down ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvLFStats1down ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvLFStats2down ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvCErr1down    ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2) +
          pow(gt.sf_csvWeights[GeneralTree::csvCErr2down    ] / gt.sf_csvWeights[GeneralTree::csvCent] - 1., 2)
        ));
      }
      /*if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
        // WH boosted weights
        nBytesRead+=bLoad(b["sf_btag0"],ientry);
        nBytesRead+=bLoad(b["sf_btag0BUp"],ientry);
        nBytesRead+=bLoad(b["sf_btag0BDown"],ientry);
        nBytesRead+=bLoad(b["sf_btag0MUp"],ientry);
        nBytesRead+=bLoad(b["sf_btag0MDown"],ientry);
        weight_btag0BUp = weight * (*sf_btag0BUp);
        weight *= (*sf_btag0);
      }*/

    }
    // Category Assignment
    switch(sample) {
      case kData:
        theCategory=kPlotData; break;
      case kQCD:
        theCategory=kPlotQCD;  break;
      case kVZ:
        nBytesRead+=bLoad(b["nB"],ientry);
        if(gt.nB>=2) theCategory=kPlotVZbb;
        else      theCategory=kPlotVVLF;
        break;
      case kWW:
        theCategory=kPlotVVLF; break;
      case kTT:
        theCategory=kPlotTT; break;
      case kTop:
        theCategory=kPlotTop; break;
      case kWjets:
        nBytesRead+=bLoad(b["nB"],ientry);
        if(gt.nB>=2) theCategory=kPlotWbb;
        else if(gt.nB==1) theCategory=kPlotWb;
        else theCategory=kPlotWLF; // light flavor
        break;
      case kZjets:
        nBytesRead+=bLoad(b["nB"],ientry);
        if(gt.nB>=2) theCategory=kPlotZbb;
        else if(gt.nB==1) theCategory=kPlotZb;
        else theCategory=kPlotZLF; // light flavor
        break;
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

