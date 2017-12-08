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
#include "PandaAnalysis/Flat/interface/GeneralTree.h"
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
  //const float bDiscrLoose = 0.5426, bDiscrMedium = 0.8484, bDiscrTight  = 0.9535; //csv
  const float bDiscrLoose = -0.5884, bDiscrMedium = 0.4432, bDiscrTight  = 0.9432; //cmva
  
  // Load Shared Objects for ACLIC
  bool loadPandaAnalysisFlat=(0==gSystem->Load("libPandaAnalysisFlat.so"));
  if(!loadPandaAnalysisFlat) { throw std::runtime_error("Error loading shared object libPandaAnalysisFlat.so"); return false; }

  // Handle File IO
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
  //TFile *muTrigEffFile = TFile::Open("/home/dhsu/CMSSW_8_0_29/src/PandaAnalysis/data/trigger_eff/muon_trig_Run2016BtoF.root", "READ");
  //if(!muTrigEffFile && muTrigEffFile->IsOpen()) { throw std::runtime_error(""); return false; }
  //TH2D *muTrigEff = (TH2D*)muTrigEffFile->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");

  TFile *outputFile = TFile::Open(outputFileName,"RECREATE");
  TTree *plotTree = new TTree("plotTree","Tree for making plots");
  unsigned char typeLepSel, theCategory, nJet;
  unsigned selectionBits;
  float pfmet, pfmetphi, pfmetsig;
  float lepton1Pt, lepton1Eta, lepton1Phi;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1Eta, hbbJet1Phi;
  float hbbJet2Pt, hbbJet2Eta, hbbJet2Phi;
  float vectorBosonPt, hbbDijetPt, hbbDijetMass;
  float bDiscrMin, bDiscrMax;
  float deltaPhiMetLep1, deltaPhiVH;
  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_r1f2, weight_r1f5, weight_r2f1, weight_r2f2, weight_r5f1, weight_r5f5;
  
  plotTree->Branch("selectionBits"   , &selectionBits   );
  plotTree->Branch("theCategory"     , &theCategory     );
  plotTree->Branch("nJet"            , &nJet            );
  plotTree->Branch("weight"          , &weight          );
  plotTree->Branch("pfmet"           , &pfmet           );
  plotTree->Branch("pfmetsig"        , &pfmetsig        );
  plotTree->Branch("pfmetphi"        , &pfmetphi        );
  plotTree->Branch("hbbDijetPt"      , &hbbDijetPt      );
  plotTree->Branch("hbbDijetMass"    , &hbbDijetMass    );
  plotTree->Branch("bDiscrMin"       , &bDiscrMin       );
  plotTree->Branch("bDiscrMax"       , &bDiscrMax       );
  plotTree->Branch("hbbJet1Pt"       , &hbbJet1Pt       );
  plotTree->Branch("hbbJet1Eta"      , &hbbJet1Eta      );
  plotTree->Branch("hbbJet1Phi"      , &hbbJet1Phi      );
  plotTree->Branch("hbbJet2Pt"       , &hbbJet2Pt       );
  plotTree->Branch("hbbJet2Eta"      , &hbbJet2Eta      );
  plotTree->Branch("hbbJet2Phi"      , &hbbJet2Phi      );
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    plotTree->Branch("typeLepSel"      , &typeLepSel      );
    plotTree->Branch("lepton1Pt"       , &lepton1Pt       );
    plotTree->Branch("lepton1Eta"      , &lepton1Eta      );
    plotTree->Branch("lepton1Phi"      , &lepton1Phi      );
    plotTree->Branch("lepton2Pt"       , &lepton2Pt       );
    plotTree->Branch("lepton2Eta"      , &lepton2Eta      );
    plotTree->Branch("lepton2Phi"      , &lepton2Phi      );
    plotTree->Branch("vectorBosonPt"   , &vectorBosonPt   );
    plotTree->Branch("deltaPhiMetLep1" , &deltaPhiMetLep1 );
    plotTree->Branch("deltaPhiVH"      , &deltaPhiVH      );
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
    
    // Selection
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
      // Preselection
      {
        // Lepton multiplicity
        nBytesRead+=bLoad(b["nLooseLep"],ientry);
        if(gt.nLooseLep<1) continue; //N_al = 0
        if(debug) printf("Passed lepton multiplicity\n");

        nBytesRead+=bLoad(b["metFilter"],ientry);
        if(gt.metFilter!=1) continue;
        if(debug) printf("Passed MET filters\n");
     
        // Jet multiplicity
        nBytesRead+=bLoad(b["nJet"],ientry);
        if     (gt.nJet<2) continue;
        // Jet kinematics
        nBytesRead+=bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
        nBytesRead+=bLoad(b["jetPt"],ientry);
        if(debug) printf("hbb jet1 pt %.2f, jet2 pt %.2f\n", gt.jetPt[gt.hbbjtidx[0]], gt.jetPt[gt.hbbjtidx[1]]);
        if(gt.jetPt[gt.hbbjtidx[0]]<25 || gt.jetPt[gt.hbbjtidx[1]]<25) continue;
        
        nBytesRead+=bLoad(b["hbbpt"],ientry);
        nBytesRead+=bLoad(b["hbbeta"],ientry);
        nBytesRead+=bLoad(b["hbbphi"],ientry);
        nBytesRead+=bLoad(b["hbbm"],ientry);
        nBytesRead+=bLoad(b["jetEta"],ientry);
        nBytesRead+=bLoad(b["jetPhi"],ientry);
        hbbDijetPt = gt.hbbpt;
        hbbDijetMass = gt.hbbm;
        hbbJet1Pt=gt.jetPt[gt.hbbjtidx[0]]; hbbJet1Eta=gt.jetEta[gt.hbbjtidx[0]]; hbbJet1Phi=gt.jetPhi[gt.hbbjtidx[0]];
        hbbJet2Pt=gt.jetPt[gt.hbbjtidx[1]]; hbbJet1Eta=gt.jetEta[gt.hbbjtidx[1]]; hbbJet1Phi=gt.jetPhi[gt.hbbjtidx[1]];
        if(debug) printf("hbbDijetPt %.2f, hbbDijetMass %.2f\n", hbbDijetPt, hbbDijetMass); 
        if(hbbDijetPt<50.) continue;
        //bool passDijetMass=false;
        if(hbbDijetMass<0 || hbbDijetMass>250.) continue;
        //if(!passDijetMass) continue;
        if(debug) printf("passed jet kinematics\n");

        // Lepton ID and isolation
        nBytesRead+=bLoad(b["nTightMuon"],ientry);
        if(gt.nTightMuon==1) typeLepSel=1; else {
          nBytesRead+=bLoad(b["nTightElectron"],ientry);
          if(gt.nTightElectron==1) typeLepSel=2;
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
        pfmet = gt.pfmet;
        
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
      
      // Other stuff
      nBytesRead+=bLoad(b["hbbpt"],ientry);
      nBytesRead+=bLoad(b["hbbphi"],ientry);
      { 
        TVector2 hbbDijetV2; hbbDijetV2.SetMagPhi(gt.hbbpt, gt.hbbphi);
        deltaPhiVH = vectorBosonV2.DeltaPhi(hbbDijetV2);
      }
      nBytesRead+=bLoad(b["pfmetsig"],ientry);
      pfmetsig = gt.pfmetsig;
      


    }
    
    // Set Selection Bits
    selectionBits=0;
    bool whCommonCuts =
      /* 2nd lepton veto     */ gt.nLooseLep==1 &&
      /* W pT > 100          */ vectorBosonPt>100 &&
      /* pTjj > 100          */ hbbDijetPt>100 && 
      /* pTe(m) > 30(25)     */ ((typeLepSel==1 && lepton1Pt>25) || (typeLepSel==2 && lepton1Pt>30)) &&
      /* deltaPhiVH > 2.5    */ deltaPhiVH > 2.5 &&
      /* deltaPhiMetLep1 < 2 */ deltaPhiMetLep1 < 2
    ;

    if(
      pfmetsig > 2 &&
      (bDiscrMax >= bDiscrLoose) && (bDiscrMax < bDiscrMedium)
    ) selectionBits |= kWHLightFlavorCR;
    if(
      nJet==2 &&
      (bDiscrMax >= bDiscrTight) && //(bDiscrMin >= bDiscrLoose ) &&
      ((hbbDijetMass < 90) || (hbbDijetMass >= 150))
    ) selectionBits |= kWHHeavyFlavorCR;
    if(
      nJet>=4 &&
      pfmetsig > 2 &&
      (bDiscrMax >= bDiscrTight) //&& (bDiscrMin >= bDiscrLoose ) 
    ) selectionBits |= kWH2TopCR;
    if(
      nJet<4 &&
      (bDiscrMax >= bDiscrTight) && (bDiscrMin >= bDiscrLoose ) &&
      ((hbbDijetMass >= 90) && (hbbDijetMass <  150))
    ) selectionBits |= kWHSR;

    // Weighting
    if(sample==kData)
      weight=1;
    else if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
      nBytesRead+=bLoad(b["normalizedWeight"],ientry);
	  nBytesRead+=bLoad(b["sf_pu"],ientry);
	  if(sample==kWjets || sample==kZjets) { nBytesRead+=bLoad(b["sf_qcdV"],ientry); nBytesRead+=bLoad(b["sf_ewkV"],ientry); } 
      else { gt.sf_qcdV=1; gt.sf_ewkV=1; }
	  //nBytesRead+=bLoad(sf_btag1,ientry);
 	  if(sample==kTT) nBytesRead+=bLoad(b["sf_tt"],ientry); else gt.sf_tt=1;
      weight = normalizedWeight * gt.sf_pu * gt.sf_ewkV * gt.sf_qcdV * gt.sf_tt;
      if (typeLepSel==1) {
        nBytesRead+=bLoad(b["muonSfReco"],ientry);
        nBytesRead+=bLoad(b["muonSfTight"],ientry);
        nBytesRead+=bLoad(b["sf_muTrig"],ientry);
        weight *= gt.sf_muTrig * gt.muonSfReco[0] * gt.muonSfTight[0];
        // need muon trigger efficiency
        //double theTrigEff = muTrigEff->GetBinContent( muTrigEff->FindBin( fabs(lepton1Eta), TMath::Max((float)26., TMath::Min((float)500., lepton1Pt)) ) );
        //weight *= theTrigEff;
      } else if(typeLepSel==2) {
        nBytesRead+=bLoad(b["electronSfReco"],ientry);
        nBytesRead+=bLoad(b["electronSfTight"],ientry);
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
        theCategory=-1;
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
