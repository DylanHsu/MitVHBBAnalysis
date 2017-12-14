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
  delete dummyTree;
   
  // temporary
  TFile *muTrigEffFile = TFile::Open("/home/dhsu/CMSSW_8_0_29/src/PandaAnalysis/data/trigger_eff/muon_trig_Run2016BtoF.root", "READ");
  if(!muTrigEffFile && muTrigEffFile->IsOpen()) { throw std::runtime_error(""); return false; }
  TH2D *muTrigEff = (TH2D*)muTrigEffFile->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");

  TFile *outputFile = TFile::Open(outputFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  TTree *plotTree = new TTree("plotTree","Tree for making plots");
  unsigned char typeLepSel, theCategory, nJet, nJet_jesUp, nJet_jesDown;
  unsigned selectionBits, nMinusOneBits, selectionBits_jesUp, selectionBits_jesDown;
  float pfmet, pfmetphi, pfmetsig, pfmetUp, pfmetDown;
  float lepton1Pt, lepton1Eta, lepton1Phi;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1PtUp, hbbJet1PtDown, hbbJet1Eta, hbbJet1Phi;
  float hbbJet2Pt, hbbJet2PtUp, hbbJet2PtDown, hbbJet2Eta, hbbJet2Phi;
  float vectorBosonPt, hbbDijetPt, hbbDijetPtUp, hbbDijetPtDown, hbbDijetMass, hbbDijetMassUp, hbbDijetMassDown;
  float bDiscrMin, bDiscrMax;
  float deltaPhiMetLep1, deltaPhiVH, deltaPhiVH_jesUp, deltaPhiVH_jesDown;
  float topMass, topMass_jesUp, topMass_jesDown;
  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_scaleUp, weight_scaleDown;
  float weight_lepSFUp;
  float weight_cmvaUp, weight_cmvaDown; 
  
  plotTree->Branch("selectionBits"         , &selectionBits         );
  plotTree->Branch("selectionBits_jesUp"   , &selectionBits_jesUp   );
  plotTree->Branch("selectionBits_jesDown" , &selectionBits_jesDown );
  plotTree->Branch("nMinusOneBits"         , &nMinusOneBits         );
  plotTree->Branch("theCategory"           , &theCategory           );
  plotTree->Branch("nJet"                  , &nJet                  );
  plotTree->Branch("nJet_jesUp"            , &nJet_jesUp            );
  plotTree->Branch("nJet_jesDown"          , &nJet_jesDown          );
  plotTree->Branch("pfmet"                 , &pfmet                 );
  plotTree->Branch("pfmetUp"               , &pfmetUp               );
  plotTree->Branch("pfmetDown"             , &pfmetDown             );
  plotTree->Branch("pfmetsig"              , &pfmetsig              );
  plotTree->Branch("pfmetphi"              , &pfmetphi              );
  plotTree->Branch("bDiscrMin"             , &bDiscrMin             );
  plotTree->Branch("bDiscrMax"             , &bDiscrMax             );
  plotTree->Branch("hbbDijetPt"            , &hbbDijetPt            );
  plotTree->Branch("hbbDijetPtUp"          , &hbbDijetPtUp          );
  plotTree->Branch("hbbDijetPtDown"        , &hbbDijetPtDown        );
  plotTree->Branch("hbbDijetMass"          , &hbbDijetMass          );
  plotTree->Branch("hbbDijetMassUp"        , &hbbDijetMassUp        );
  plotTree->Branch("hbbDijetMassDown"      , &hbbDijetMassDown      );
  plotTree->Branch("hbbJet1Pt"             , &hbbJet1Pt             );
  plotTree->Branch("hbbJet1PtUp"           , &hbbJet1PtUp           );
  plotTree->Branch("hbbJet1PtDown"         , &hbbJet1PtDown         );
  plotTree->Branch("hbbJet1Eta"            , &hbbJet1Eta            );
  plotTree->Branch("hbbJet1Phi"            , &hbbJet1Phi            );
  plotTree->Branch("hbbJet2Pt"             , &hbbJet2Pt             );
  plotTree->Branch("hbbJet2PtUp"           , &hbbJet2PtUp           );
  plotTree->Branch("hbbJet2PtDown"         , &hbbJet2PtDown         );
  plotTree->Branch("hbbJet2Eta"            , &hbbJet2Eta            );
  plotTree->Branch("hbbJet2Phi"            , &hbbJet2Phi            );
  plotTree->Branch("weight"                , &weight                );
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    plotTree->Branch("typeLepSel"              , &typeLepSel               );
    plotTree->Branch("lepton1Pt"               , &lepton1Pt                );
    plotTree->Branch("lepton1Eta"              , &lepton1Eta               );
    plotTree->Branch("lepton1Phi"              , &lepton1Phi               );
    plotTree->Branch("vectorBosonPt"           , &vectorBosonPt            );
    plotTree->Branch("deltaPhiMetLep1"         , &deltaPhiMetLep1          );
    plotTree->Branch("deltaPhiVH"              , &deltaPhiVH               );
    plotTree->Branch("deltaPhiVH_jesUp"        , &deltaPhiVH_jesUp         );
    plotTree->Branch("deltaPhiVH_jesDown"      , &deltaPhiVH_jesDown       );
    plotTree->Branch("topMass"                 , &topMass                  );
    plotTree->Branch("topMass_jesUp"           , &topMass_jesUp            );
    plotTree->Branch("topMass_jesDown"         , &topMass_jesDown          );
    plotTree->Branch("weight_pdfUp"            , &weight_pdfUp             );
    plotTree->Branch("weight_pdfDown"          , &weight_pdfDown           );
    plotTree->Branch("weight_scaleUp"          , &weight_scaleUp           );
    plotTree->Branch("weight_scaleDown"        , &weight_scaleDown         );
    plotTree->Branch("weight_lepSFUp"          , &weight_lepSFUp           );
    plotTree->Branch("weight_cmvaUp"           , &weight_cmvaUp            );
    plotTree->Branch("weight_cmvaDown"         , &weight_cmvaDown          );
  }
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
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
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
        if     (gt.nJet<2) continue;
        nJet=gt.nJet; nJet_jesUp=gt.nJet; nJet_jesDown=gt.nJet;
        // Jet kinematics
        nBytesRead+=bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
        nBytesRead+=bLoad(b["jetPt"],ientry);
        if(debug) printf("hbb jet1 pt %.2f, jet2 pt %.2f\n", gt.jetPt[gt.hbbjtidx[0]], gt.jetPt[gt.hbbjtidx[1]]);
        if(gt.jetPt[gt.hbbjtidx[0]]<25 || gt.jetPt[gt.hbbjtidx[1]]<25) continue;
        
        nBytesRead+=bLoad(b["hbbpt_reg"],ientry);
        nBytesRead+=bLoad(b["hbbpt_jesUp"],ientry);
        nBytesRead+=bLoad(b["hbbpt_jesDown"],ientry);
        nBytesRead+=bLoad(b["hbbeta"],ientry);
        nBytesRead+=bLoad(b["hbbphi"],ientry);
        nBytesRead+=bLoad(b["hbbm_reg"],ientry);
        nBytesRead+=bLoad(b["hbbm_jesUp"],ientry);
        nBytesRead+=bLoad(b["hbbm_jesDown"],ientry);
        nBytesRead+=bLoad(b["jetEta"],ientry);
        nBytesRead+=bLoad(b["jetPhi"],ientry);
        nBytesRead+=bLoad(b["jetPtUp"],ientry);
        nBytesRead+=bLoad(b["jetPtDown"],ientry);
        nBytesRead+=bLoad(b["jetRegFac"],ientry);
        hbbDijetPt = gt.hbbpt_reg;
        hbbDijetPtUp = gt.hbbpt_jesUp;
        hbbDijetPtDown = gt.hbbpt_jesDown;
        hbbDijetMass = gt.hbbm_reg;
        hbbDijetMassUp = gt.hbbm_jesUp;
        hbbDijetMassDown = gt.hbbm_jesDown;
        hbbJet1Pt=gt.jetRegFac[gt.hbbjtidx[0]] * gt.jetPt[gt.hbbjtidx[0]]; 
        hbbJet1PtUp=gt.jetRegFac[gt.hbbjtidx[0]] * gt.jetPtUp[gt.hbbjtidx[0]]; 
        hbbJet1PtDown=gt.jetRegFac[gt.hbbjtidx[0]] * gt.jetPtDown[gt.hbbjtidx[0]]; 
        hbbJet1Eta=gt.jetEta[gt.hbbjtidx[0]];
        hbbJet1Phi=gt.jetPhi[gt.hbbjtidx[0]];
        hbbJet2Pt=gt.jetRegFac[gt.hbbjtidx[1]] * gt.jetPt[gt.hbbjtidx[1]];
        hbbJet2PtUp=gt.jetRegFac[gt.hbbjtidx[1]] * gt.jetPtUp[gt.hbbjtidx[1]]; 
        hbbJet2PtDown=gt.jetRegFac[gt.hbbjtidx[1]] * gt.jetPtDown[gt.hbbjtidx[1]]; 
        hbbJet1Eta=gt.jetEta[gt.hbbjtidx[1]];
        hbbJet1Phi=gt.jetPhi[gt.hbbjtidx[1]];
        if(debug) printf("hbbDijetPt %.2f, hbbDijetMass %.2f\n", hbbDijetPt, hbbDijetMass); 
        if(hbbDijetPt<50. && hbbDijetPtUp<50.) continue;
        if((hbbDijetMass<0 && hbbDijetMassDown<0)|| (hbbDijetMass>250. && hbbDijetMassUp>250.)) continue;
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
          lepton1Pt = gt.muonPt[0]; lepton1Eta = gt.muonEta[0]; lepton1Phi = gt.muonPhi[0];
        } else if(typeLepSel==2) {
          nBytesRead+=bLoad(b["electronPt"],ientry);
          nBytesRead+=bLoad(b["electronEta"],ientry);
          nBytesRead+=bLoad(b["electronPhi"],ientry);
          lepton1Pt = gt.electronPt[0]; lepton1Eta = gt.electronEta[0]; lepton1Phi = gt.electronPhi[0];
        } else continue;
        if(debug) printf("Passed lepton kinematics\n");
        
        // Met
        nBytesRead+=bLoad(b["pfmet"],ientry);
        nBytesRead+=bLoad(b["pfmetUp"],ientry);
        nBytesRead+=bLoad(b["pfmetDown"],ientry);
        pfmet = gt.pfmet; pfmetUp = gt.pfmetUp; pfmetDown = gt.pfmetDown;
        
        //if(gt.pfmet<30) continue;
        //if(debug) printf("passed MET\n");
        
        // W reconstruction
        nBytesRead+=bLoad(b["pfmetphi"],ientry);
        pfmetphi = gt.pfmetphi;
        { 
          TVector2 lepton1V2, pfmetV2;
          lepton1V2.SetMagPhi(lepton1Pt,lepton1Phi);
          pfmetV2.SetMagPhi(pfmet,pfmetphi);
          vectorBosonV2 = lepton1V2+pfmetV2; 
          deltaPhiMetLep1 = pfmetV2.DeltaPhi(lepton1V2);
        } vectorBosonPt=vectorBosonV2.Mod();
        if(vectorBosonPt<50.) continue;
        if(debug) printf("passed W reco\n");
      } // end preselection
      
      // Jet B-tagging
      nBytesRead+=bLoad(b["jetCMVA"],ientry);
      bDiscrMax=gt.jetCMVA[gt.hbbjtidx[0]];
      bDiscrMin=gt.jetCMVA[gt.hbbjtidx[1]];

      // Top reconstruction
      {
        nBytesRead+=bLoad(b["pfmet"],ientry);
        nBytesRead+=bLoad(b["pfmetUp"],ientry);
        nBytesRead+=bLoad(b["pfmetDown"],ientry);
        nBytesRead+=bLoad(b["jetE"],ientry);
        nBytesRead+=bLoad(b["jetPt"],ientry);
        nBytesRead+=bLoad(b["jetEta"],ientry);
        nBytesRead+=bLoad(b["jetPhi"],ientry);
        nBytesRead+=bLoad(b["jetPtUp"],ientry);
        nBytesRead+=bLoad(b["jetPtDown"],ientry);
        nBytesRead+=bLoad(b["jetRegFac"],ientry);
        nBytesRead+=bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
        TLorentzVector leptonP4, metP4, nuP4, jet1P4, jet2P4, WP4, topP4;
        bool jet1IsCloser; float jet1Mass, jet2Mass, dRJet1W, dRJet2W;
        leptonP4.SetPtEtaPhiM( lepton1Pt, lepton1Eta, lepton1Phi, (typeLepSel==1? 0.106:511e-6));
        // Nominal jet energies
        metP4.SetPtEtaPhiM( gt.pfmet, 0, gt.pfmetphi, 0 );
        nuP4 = getNu4Momentum( leptonP4, metP4 );
        WP4 = leptonP4 + nuP4;
        jet1P4.SetPtEtaPhiE(gt.jetPt[gt.hbbjtidx[0]],gt.jetEta[gt.hbbjtidx[0]],gt.jetPhi[gt.hbbjtidx[0]],gt.jetE[gt.hbbjtidx[0]]); 
        jet2P4.SetPtEtaPhiE(gt.jetPt[gt.hbbjtidx[1]],gt.jetEta[gt.hbbjtidx[1]],gt.jetPhi[gt.hbbjtidx[1]],gt.jetE[gt.hbbjtidx[1]]);
        jet1Mass=jet1P4.M(); jet2Mass=jet2P4.M();
        jet1P4.SetPtEtaPhiM(gt.jetPt[gt.hbbjtidx[0]]*gt.jetRegFac[gt.hbbjtidx[0]],gt.jetEta[gt.hbbjtidx[0]],gt.jetPhi[gt.hbbjtidx[0]], jet1Mass); 
        jet2P4.SetPtEtaPhiM(gt.jetPt[gt.hbbjtidx[1]]*gt.jetRegFac[gt.hbbjtidx[1]],gt.jetEta[gt.hbbjtidx[1]],gt.jetPhi[gt.hbbjtidx[1]], jet2Mass);
        dRJet1W=jet1P4.DeltaR(leptonP4); dRJet2W=jet2P4.DeltaR(leptonP4);
        jet1IsCloser = (dRJet1W < dRJet2W);
        topP4 = jet1IsCloser? jet1P4+WP4 : jet2P4+WP4;
        topMass = topP4.M();
        if(debug) printf("reconstructed a top mass of %.1f GeV from lepton (pT, eta, phi)=(%.1f,%.2f,%.2f), MET (pT,phi)=(%.1f,%.2f), jet (pT,eta,phi,M)=(%.1f,%.2f,%.2f,%.1f), nu (pT,eta,phi)=(%.1f, %.2f, %.2f), W (pT,eta,phi,M)=(%.1f,%.2f,%.2f,%.1f)\n", topMass, lepton1Pt, lepton1Eta, lepton1Phi, gt.pfmet, gt.pfmetphi, jet1IsCloser? jet1P4.Pt():jet2P4.Pt(), jet1IsCloser? jet1P4.Eta():jet2P4.Eta(), jet1IsCloser? jet1P4.Phi():jet2P4.Phi(), jet1IsCloser? jet1P4.M():jet2P4.M(), nuP4.Pt(),nuP4.Eta(),nuP4.Phi(), WP4.Pt(), WP4.Eta(), WP4.Phi(), WP4.M());
        // Top mass with jet energies varied up
        metP4.SetPtEtaPhiM( gt.pfmetUp, 0, gt.pfmetphi, 0 );
        nuP4 = getNu4Momentum( leptonP4, metP4 );
        WP4 = leptonP4 + nuP4;
        jet1P4.SetPtEtaPhiM(
          gt.jetPtUp[gt.hbbjtidx[0]]*gt.jetRegFac[gt.hbbjtidx[0]],
          gt.jetEta[gt.hbbjtidx[0]],
          gt.jetPhi[gt.hbbjtidx[0]],
          jet1Mass
        ); 
        jet2P4.SetPtEtaPhiM(
          gt.jetPtUp[gt.hbbjtidx[1]]*gt.jetRegFac[gt.hbbjtidx[1]],
          gt.jetEta[gt.hbbjtidx[1]],
          gt.jetPhi[gt.hbbjtidx[1]],
          jet2Mass
        );
        dRJet1W=jet1P4.DeltaR(leptonP4); dRJet2W=jet2P4.DeltaR(leptonP4);
        jet1IsCloser = (dRJet1W < dRJet2W);
        topP4 = jet1IsCloser? jet1P4+WP4 : jet2P4+WP4;
        topMass_jesUp = topP4.M();
        // Top mass with jet energies varied down
        metP4.SetPtEtaPhiM( gt.pfmetDown, 0, gt.pfmetphi, 0 );
        nuP4 = getNu4Momentum( leptonP4, metP4 );
        WP4 = leptonP4 + nuP4;
        jet1P4.SetPtEtaPhiM(
          gt.jetPtDown[gt.hbbjtidx[0]]*gt.jetRegFac[gt.hbbjtidx[0]],
          gt.jetEta[gt.hbbjtidx[0]],
          gt.jetPhi[gt.hbbjtidx[0]],
          jet1Mass
        ); 
        jet2P4.SetPtEtaPhiM(
          gt.jetPtDown[gt.hbbjtidx[1]]*gt.jetRegFac[gt.hbbjtidx[1]],
          gt.jetEta[gt.hbbjtidx[1]],
          gt.jetPhi[gt.hbbjtidx[1]],
          jet2Mass
        );
        dRJet1W=jet1P4.DeltaR(leptonP4); dRJet2W=jet2P4.DeltaR(leptonP4);
        jet1IsCloser = (dRJet1W < dRJet2W);
        topP4 = jet1IsCloser? jet1P4+WP4 : jet2P4+WP4;
        topMass_jesDown = topP4.M();
      }

      // VH reconstruction
      { 
        nBytesRead+=bLoad(b["hbbpt_reg"    ],ientry);
        nBytesRead+=bLoad(b["hbbpt_jesUp"  ],ientry);
        nBytesRead+=bLoad(b["hbbpt_jesDown"],ientry);
        nBytesRead+=bLoad(b["hbbphi"],ientry);
        TVector2 hbbDijetV2;
        hbbDijetV2.SetMagPhi(gt.hbbpt_reg, gt.hbbphi);
        deltaPhiVH = vectorBosonV2.DeltaPhi(hbbDijetV2);
        hbbDijetV2.SetMagPhi(gt.hbbpt_jesUp, gt.hbbphi);
        deltaPhiVH_jesUp = vectorBosonV2.DeltaPhi(hbbDijetV2);
        hbbDijetV2.SetMagPhi(gt.hbbpt_jesDown, gt.hbbphi);
        deltaPhiVH_jesDown = vectorBosonV2.DeltaPhi(hbbDijetV2);
      }
      nBytesRead+=bLoad(b["pfmetsig"],ientry);
      pfmetsig = gt.pfmetsig;

    }
    
    // Set Selection Bits
    selectionBits=0; nMinusOneBits=0;
    std::map<TString, bool> cut, cut_jesUp, cut_jesDown;;
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
      vector<TString> cutsWHLightFlavorCR, cutsWHHeavyFlavorCR, cutsWH2TopCR, cutsWHSR;
      cutsWHLightFlavorCR ={"2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiMetLep1","looseBTag","tightBVeto","metSig"};
      cutsWHHeavyFlavorCR ={"2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiMetLep1","nJet_WHHF","tightBTag","mH_Flip","metSig"};
      cutsWH2TopCR        ={"2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiMetLep1","nJet_WHTT","tightBTag"};
      cutsWHSR            ={"2ndLepVeto","WpT","pTjj","lepton1Pt","dPhiVH","dPhiMetLep1","nJet_WHSR","tightBTag","looseBTag2","mH_WHSR"};

      cut["2ndLepVeto" ] = gt.nLooseLep==1;
      cut["WpT"        ] = vectorBosonPt>100;
      cut["pTjj"       ] = hbbDijetPt>100; 
      cut["lepton1Pt"  ] = ((typeLepSel==1 && lepton1Pt>25) || (typeLepSel==2 && lepton1Pt>30));
      cut["dPhiVH"     ] = deltaPhiVH > 2.5;
      cut["dPhiMetLep1"] = deltaPhiMetLep1 < 2;
      cut["tightBTag"  ] = (bDiscrMax >= bDiscrTight);
      cut["tightBVeto" ] = (bDiscrMax < bDiscrTight);
      cut["looseBTag"  ] = (bDiscrMax >= bDiscrLoose);
      cut["looseBTag2" ] = (bDiscrMin >= bDiscrLoose);
      cut["metSig"     ] = pfmetsig>2;
      cut["nJet_WHTT"  ] = nJet>=4;
      cut["nJet_WHHF"  ] = nJet==2;
      cut["nJet_WHSR"  ] = nJet<4;
      cut["mH_WHSR"    ] = ((hbbDijetMass >= 90) && (hbbDijetMass <  150));
      cut["mH_Flip"    ] = ((hbbDijetMass < 90) || (hbbDijetMass >= 150));

      cut_jesUp=cut;
      cut_jesUp["pTjj"       ] = hbbDijetPtUp>100; 
      cut_jesUp["dPhiVH"     ] = deltaPhiVH_jesUp > 2.5;
      cut_jesUp["nJet_WHTT"  ] = nJet_jesUp>=4;
      cut_jesUp["nJet_WHHF"  ] = nJet_jesUp==2;
      cut_jesUp["nJet_WHSR"  ] = nJet_jesUp<4;
      cut_jesUp["mH_WHSR"    ] = ((hbbDijetMassUp >= 90) && (hbbDijetMassUp <  150));
      cut_jesUp["mH_Flip"    ] = ((hbbDijetMassUp < 90) || (hbbDijetMassUp >= 150));
      
      cut_jesDown=cut; 
      cut_jesDown["pTjj"       ] = hbbDijetPtDown>100; 
      cut_jesDown["dPhiVH"     ] = deltaPhiVH_jesDown > 2.5;
      cut_jesDown["nJet_WHTT"  ] = nJet_jesDown>=4;
      cut_jesDown["nJet_WHHF"  ] = nJet_jesDown==2;
      cut_jesDown["nJet_WHSR"  ] = nJet_jesDown<4;
      cut_jesDown["mH_WHSR"    ] = ((hbbDijetMassDown >= 90) && (hbbDijetMassDown <  150));
      cut_jesDown["mH_Flip"    ] = ((hbbDijetMassDown < 90) || (hbbDijetMassDown >= 150));
      
      if(passAllCuts( cut, cutsWHLightFlavorCR))   selectionBits |= kWHLightFlavorCR;
      if(passAllCuts( cut, cutsWHHeavyFlavorCR))   selectionBits |= kWHHeavyFlavorCR;
      if(passAllCuts( cut, cutsWH2TopCR       ))   selectionBits |= kWH2TopCR;
      if(passAllCuts( cut, cutsWHSR           ))   selectionBits |= kWHSR;
        
      if(passAllCuts( cut_jesUp, cutsWHLightFlavorCR))   selectionBits_jesUp |= kWHLightFlavorCR;
      if(passAllCuts( cut_jesUp, cutsWHHeavyFlavorCR))   selectionBits_jesUp |= kWHHeavyFlavorCR;
      if(passAllCuts( cut_jesUp, cutsWH2TopCR       ))   selectionBits_jesUp |= kWH2TopCR;
      if(passAllCuts( cut_jesUp, cutsWHSR           ))   selectionBits_jesUp |= kWHSR;
        
      if(passAllCuts( cut_jesDown, cutsWHLightFlavorCR))   selectionBits_jesDown |= kWHLightFlavorCR;
      if(passAllCuts( cut_jesDown, cutsWHHeavyFlavorCR))   selectionBits_jesDown |= kWHHeavyFlavorCR;
      if(passAllCuts( cut_jesDown, cutsWH2TopCR       ))   selectionBits_jesDown |= kWH2TopCR;
      if(passAllCuts( cut_jesDown, cutsWHSR           ))   selectionBits_jesDown |= kWHSR;
        
      if(passNMinusOne( cut, cutsWHLightFlavorCR)) nMinusOneBits |= kWHLightFlavorCR;
      if(passNMinusOne( cut, cutsWHHeavyFlavorCR)) nMinusOneBits |= kWHHeavyFlavorCR;
      if(passNMinusOne( cut, cutsWH2TopCR       )) nMinusOneBits |= kWH2TopCR;
      if(passNMinusOne( cut, cutsWHSR           )) nMinusOneBits |= kWHSR;
    }
    if(nMinusOneBits==0) continue;

    // Weighting
    if(sample==kData) weight=1;
    else if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
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
        //weight *= gt.sf_muTrig * gt.muonSfReco[0] * gt.muonSfTight[0];
        weight *= gt.muonSfReco[0] * gt.muonSfTight[0];
        // need muon trigger efficiency
        double theTrigEff = muTrigEff->GetBinContent( muTrigEff->FindBin( fabs(lepton1Eta), TMath::Max((float)26., TMath::Min((float)500., lepton1Pt)) ) );
        weight *= theTrigEff;
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
      
      // CMVA weight
      nBytesRead+=bLoad(b["sf_cmvaWeight_Cent"],ientry);
      weight *= gt.sf_csvWeights[GeneralTree::csvCent];
      
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
      
      // Lepton ID & ISO SF uncertainty
      if(typeLepSel==1) weight_lepSFUp = weight * (1.+gt.muonSfUnc[0]);
      else weight_lepSFUp = weight * (1.+gt.electronSfUnc[0]);
      
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
