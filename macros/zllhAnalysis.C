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
// useBoostedCategory - This controls whether to cut on the number of fatjets in the resolved regions
// False means you will use all possible events to do the resolved VH
// True means the events with a fatjet are reserved for boosted VH

const bool useHtBinnedVJetsKFactor=true;
const int NJES = (int)shiftjes::N;

using namespace vhbbPlot;
void vhbbPlotSkim(
  TString dataCardDir,
  bool useBoostedCategory=false,
  int MVAVarType=3,
  unsigned year=2016,
  bool debug=false
) {
  assert(dataCardDir!="");
  if(dataCardDir!="") system(Form("mkdir -p MitVHBBAnalysis/datacards/%s",dataCardDir.Data()));
  const double theLumi=35900.;
  
  // Choice of the MVA variable type, binning, and the name
  // This can be different for each control region
  std::map<selectionType, vector<float>> MVAbins;
  std::map<selectionType, TString> MVAVarName, shapeType; // title of histos, and simple name
  
  vector<selectionType> selections = { // List of enums for the selections we care about, defined in vhbbPlot.h
    kZllHLightFlavorCR   ,
    kZllHHeavyFlavorCR   ,
    kZllH2TopCR          ,
    kZllHSR              ,
    kZllHPresel          ,
    kZllHLightFlavorFJCR ,
    kZllHHeavyFlavorFJCR ,
    kZllHTT2bFJCR        ,
    kZllHTT1bFJCR        ,
    kZllHFJSR            ,
    kZllHFJPresel        
  };

  for(unsigned iSel=0; iSel<selections.size(); iSel++) { // Define the shape variable
    selectionType sel = selections[iSel];
    if(MVAVarType==1) {
      // 1 - simple pT variable
      MVAVarName[sel]="Higgs p_{T} classifier [GeV]";
      if(sel>=kZllHLightFlavorCR && sel<=kZllHPresel) {
        MVAbins[sel]={100,120,140,160,180,200,250,300,350};
        MVAVarName[sel]="H(bb) pT";
        shapeType[sel]="ptShape";
      } else if(sel>=kZllHLightFlavorFJCR && sel<=kZllHFJPresel) {
        MVAbins[sel]={250,300,350,400,450,500,550,600};
        MVAVarName[sel]="H(bb) pT";
        shapeType[sel]="ptShape";
      } 
    //} else if(MVAVarType==2) {
      // 2 - multiclass BDT in SR, subleading CMVA in CR (not implemented)
    } else if(MVAVarType==3) {
      // 3 - normal BDT in SR, subleading CMVA in CR
      if(sel==kZllHLightFlavorCR) {
        MVAbins[sel]={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000};
        MVAVarName[sel]="Subleading H(bb) CMVA";
        shapeType[sel]="lesserCMVAShape";
      } else if(sel==kZllHHeavyFlavorCR || sel==kZllH2TopCR) {
        MVAbins[sel]={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000, 0.7333, 0.8667, 1.0000};
        MVAVarName[sel]="Subleading H(bb) CMVA";
        shapeType[sel]="lesserCMVAShape";
      } else if(sel==kZllHSR) {
        MVAbins[sel]={-1,-0.5,0, 0.20,0.40,0.60,0.70,0.80,0.9,1.00};
        MVAVarName[sel]="BDT Output";
        shapeType[sel]="singleClassBDTShape"; 
      } else if(sel>=kZllHLightFlavorFJCR && sel<kZllHFJSR) {
        if(sel==kZllHHeavyFlavorFJCR)
          MVAbins[sel]={40,45,50,55,60,65,70,75,80};
        else
          MVAbins[sel]={40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,160,180,200};
        MVAVarName[sel]="Fatjet soft drop mass [GeV]";
        shapeType[sel]="softDropMassShape";
      } else if(sel==kZllHFJSR) {
        MVAbins[sel]={-0.6,-0.2,0,0.2,0.3,0.6};
        MVAVarName[sel]="BDT Output";
        shapeType[sel]="singleClassBDTShape"; 
      }
    } else throw std::runtime_error("bad MVAVarType");
  }
  
  ////////////////////////////////////////////////////////////////////////
  // Instantiate GeneralTree object for reading the ntuples
  GeneralTree gt;
  gt.is_hbb         = true;
  gt.is_fatjet      = true;
  gt.is_leptonic    = true;
  gt.btagWeights    = true;
  gt.useCMVA        = true;
  // Branches not in GeneralTree;
  std::map<TString, void*> extraAddresses;
  float normalizedWeight; extraAddresses["normalizedWeight"] = &normalizedWeight;
  // Map of the branches
  std::map<TString, TBranch*> b;
  // Dummy tree for a list of branch names
  TTree *dummyTree = new TTree("dummyTree","dummyTree");
  gt.WriteTree(dummyTree);
  // Done making GeneralTree object
  ////////////////////////////////////////////////////////////////////////
  // Load corrections to apply offline
  // CMVA reweighting for 2016
  CSVHelper *cmvaReweighter = new CSVHelper("PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root"   , "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root"   , 5);

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
  // Done loading offline corrections
  ////////////////////////////////////////////////////////////////////////
  // Declare local analysis variables
  // Isojets: For boosted categories, the number of 30 GeV AK4 jets not in the fat jet
  vector<int> nIsojet(NJES);
  vector<vector<unsigned char>> isojets(NJES);
  vector<unsigned char> isojetNBtags(NJES);
  // Selection bits
  unsigned selectionBits[NJES], nMinusOneBits;
  unsigned char typeLepSel, theCategory;
  // Declare variables for the variations of the weights
  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_QCDr1f2, weight_QCDr1f5, weight_QCDr2f1, weight_QCDr2f2, weight_QCDr5f1, weight_QCDr5f5;
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
  // Done declaring weight variables
  ////////////////////////////////////////////////////////////////////////

  float isojetBtagCut = (year==2017)? deepcsvLoose : cmvaLoose;
  std::map<selectionType,vector<TString>> cuts;
  cuts[kZllHLightFlavorCR  ] = {"ZpT","pTjj","bveto","Zmass"                       };
  cuts[kZllHHeavyFlavorCR  ] = {"ZpT","pTjj","btag" ,"ZmassTight","lowMET","dPhiZH"};
  cuts[kZllH2TopCR         ] = {"ZpT","pTjj","btag" ,"ZmassSB"                     };
  cuts[kZllHSR             ] = {"ZpT","pTjj","btag" ,"Zmass"                       };
  cuts[kZllHPresel         ] = {"ZpT","pTjj"                                       };
  cuts[kZllHLightFlavorFJCR] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "0ijb", "Zmass"               };
  cuts[kZllHHeavyFlavorFJCR] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD_SB", "0ijb", "ZmassTight", "lowMET"};
  cuts[kZllHTT1bFJCR       ] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "1ijb", "Zmass"               };
  cuts[kZllHTT2bFJCR       ] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD"   , "1ijb", "Zmass"               };
  cuts[kZllHFJSR           ] = {"ZpTFJ","pTFJ","dPhiZHFJ","mSD_SR", "0ijb", "Zmass"               };
  cuts[kZllHFJPresel       ] = {"ZpTFJ","pTFJ","dPhiZHFJ"                                         };

  // Begin Chain Loop
  //for(placeholder) { 
    // Nasty hack code to get the tree branches and their addresses, have to do it for each file
    sampleType sample = kData;
    TFile *inputFile=0;

    TTree *events = (TTree*)inputFile->Get("events");
    if(!events) { throw std::runtime_error("Problem loading tree"); }
    TObjArray *listOfBranches = events->GetListOfBranches();
    for(unsigned iB=0; iB<(unsigned)listOfBranches->GetEntries(); iB++) {
      TBranch *branch = (TBranch*)listOfBranches->At(iB);
      TString branchName = branch->GetName();
      bool isExtraBranch=false;
      auto x = extraAddresses.find(branchName);
      isExtraBranch= (x!=extraAddresses.end());
      if(isExtraBranch) {
        b[x->first] = events->GetBranch(x->first);
        if(!b[x->first]) { throw std::runtime_error(Form("Extra branch %s could not be found in events tree\n", x->first.Data())); }
        b[x->first]->SetAddress(x->second);
        if(debug) printf("Booking extra branch \"%s\"\n", x->first.Data());
        continue;
      }
      TBranch *dummyBranch = dummyTree->GetBranch(branchName);
      if(!dummyBranch) { throw std::runtime_error(Form("WARNING: Couldn't find branch \"%s\" in the dummyTree, skipping", branchName.Data())); continue; }
      void* address=dummyBranch->GetAddress();
      b[branchName] = events->GetBranch(branchName);
      if(!b[branchName]) { throw std::runtime_error(Form("Branch \"%s\" could not be found in events tree", x->first.Data())); }
      b[branchName]->SetAddress(address);
      if(debug) printf("Booking GeneralTree branch \"%s\"\n", branchName.Data());
    }
 
    Long64_t nentries = events->GetEntries();
    // Begin Event Loop
    for (Long64_t ientry=0; ientry<nentries; ientry++) {
      if(debug || ientry%100000==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries);
      //////////////////////
      // Clear variables
      theCategory=-1; // plot category 
      typeLepSel=99; // 0: mixed e-mu, 1: all mu, 2: all e, 99: undefined

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
      if(gt.nLooseLep<2) continue; 
      if(debug) printf("Passed lepton multiplicity\n");

      // Trigger
      if(sample==kData) {
        bLoad(b["trigger"],ientry);
        if( (gt.trigger & (1<<4 | 1<<5))==0 ) continue;
      }

      // Lepton ID and isolation
      bLoad(b["nLooseElectron"],ientry);
      bLoad(b["nTightElectron"],ientry);
      bLoad(b["nLooseMuon"],ientry);
      bLoad(b["nTightMuon"],ientry);
      bLoad(b["muonSelBit"],ientry);
      bLoad(b["muonPt"],ientry);
      bLoad(b["electronSelBit"],ientry);
      bLoad(b["electronPt"],ientry);

      if(gt.muonPt[0]>20 && gt.muonPt[1]>10 && 
        gt.muonPdgId[0]+gt.muonPdgId[1]==0 &&
        (gt.muonSelBit[0] & 1<<3)!=0 &&
        (gt.muonSelBit[1] & 1<<3)!=0 &&
        (gt.trigger & 1<<4)!=0
      ) typeLepSel=1;
      else if(gt.electronPt[0]>25 && gt.electronPt[1]>15 && 
        gt.electronPdgId[0]+gt.electronPdgId[1]==0 &&
        (gt.electronSelBit[0] & 1<<5)!=0 &&
        (gt.electronSelBit[1] & 1<<5)!=0 &&
        (gt.trigger & 1<<5)!=0
      ) typeLepSel=2;
      // E-Mu Selection Not implemented yet!!!
      if(typeLepSel!=1 && typeLepSel!=2) continue;
      if(debug) printf("Passed lepton ID/iso multiplicity\n");

      // Lepton kinematics
      float lepton1Pt,lepton1Eta,lepton1Phi,lepton1RelIso,lepton1D0,lepton1DZ,
            lepton2Pt,lepton2Eta,lepton2Phi,lepton2RelIso,lepton2D0,lepton2DZ;
      if (typeLepSel==1) {
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
        lepton2Pt     = gt.electronPt[1]; 
        lepton2Eta    = gt.electronEta[1];
        lepton2Phi    = gt.electronPhi[1];
        lepton2RelIso = gt.electronCombIso[1]/gt.electronPt[1];
        lepton2D0     = gt.electronD0[1];
        lepton2DZ     = gt.electronDZ[1];
      }
      if(debug) printf("Passed lepton kinematics\n");

      bLoad(b["ZBosonPt"],ientry);
      bLoad(b["ZBosonM"],ientry);
      if(gt.ZBosonPt<30) continue;
      if(gt.ZBosonM<10) continue;
      if(debug) printf("Passed Z boson reconstruction\n");
      
      // Jet multiplicity
      bLoad(b["nJet"],ientry);
      if     (gt.nJet[0]<2) continue;
      bool isBoostedCategory=false;
      if(useBoostedCategory) { 
        bLoad(b["nFatjet"],ientry);
        bLoad(b["fjMSD_corr"],ientry);
        bLoad(b["fjPt"],ientry);
        bLoad(b["fjEta"],ientry);
        bLoad(b["fjPhi"],ientry);
        float temp_deltaPhiVH = fabs(TVector2::Phi_mpi_pi(gt.ZBosonPhi-gt.fjPhi));
        if(
          gt.nFatjet != 0 && 
          gt.fjPt[0] >= 250 && 
          gt.fjMSD_corr[0] >= 40 &&
          fabs(gt.fjEta) < 2.4 &&
          gt.ZBosonPt >= 250 &&
          temp_deltaPhiVH >= 2.5
        ) isBoostedCategory=true;
      }

      // Jet kinematics
      if(isBoostedCategory) {
        // No checks here? 
      } else { 
        bLoad(b["nJot"],ientry);
        bLoad(b["hbbjtidx"],ientry); // indices of Higgs daughter jets
        bLoad(b["hbbpt"],ientry);
        bLoad(b["hbbm_reg"],ientry);
        if(
          gt.hbbpt[0]<50 || 
          gt.hbbm_reg[0]<0 ||
          gt.hbbm_reg[0]>250) 
          continue;
      }
      if(debug) printf("passed jet kinematics\n");

      // Met
      bLoad(b["pfmet"],ientry);
      //bLoad(b["pfmetphi"],ientry);
      
      // Jet B-tagging
      bool bjet1IsTight, bjet2IsLoose;
      if(year==2016) {
        bLoad(b["jotCMVA"],ientry);
        bjet1IsTight = gt.jotCMVA[gt.hbbjtidx[0][0]] > cmvaTight;
        bjet2IsLoose = gt.jotCMVA[gt.hbbjtidx[0][1]] > cmvaLoose;
      } else if(year==2017) {
        bLoad(b["jotCSV"],ientry);
        bjet1IsTight = gt.jotCSV[gt.hbbjtidx[0][0]] > deepcsvTight;
        bjet2IsLoose = gt.jotCSV[gt.hbbjtidx[0][1]] > deepcsvLoose;
      }
      
      // Isojets for boosted category
      if(isBoostedCategory) {
        bLoad(b["jotPt"],ientry);
        bLoad(b["jotEta"],ientry);
        bLoad(b["jotPhi"],ientry);
        bLoad(b["fjEta"],ientry);
        bLoad(b["fjPhi"],ientry);
        for(unsigned iJES=0; iJES<NJES; iJES++) 
        for(unsigned char iJ=0; iJ<gt.nJot[iJES]; iJ++) {
          if(fabs(gt.jotEta[iJ])>2.4) continue;
          float dR2JetFatjet=pow(gt.jotEta[iJ]-gt.fjEta,2)+pow(TVector2::Phi_mpi_pi(gt.jotPhi[iJ]-gt.fjPhi),2);
          if(dR2JetFatjet<0.64) continue;
          
          float isojetBtag = (year==2017)? gt.jotCMVA[iJ] : gt.jotCSV[iJ];
          if(gt.jotPt[iJES][iJ]>30) {
            isojets[iJES].push_back(iJ);
            if(isojetBtag>isojetBtagCut) isojetNBtags[iJES]++;
          }
          
          // CMVA jet kinematic decorrelated weight nuisances for the isojets
          // Nominal JES scenario only
          if(iJES!=0) continue;
          int iPt=-1, iEta=-1;
          double jetAbsEta=fabs(gt.jotEta[iJ]);
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
            jetPts    [iPt][iEta].push_back(gt.jotPt[iJES][iJ]);
            jetEtas   [iPt][iEta].push_back(gt.jotEta[iJ]);
            jetFlavors[iPt][iEta].push_back(gt.jotFlav[iJ]);
            jetBtags  [iPt][iEta].push_back(isojetBtag);
          }
        }
      }
      for(unsigned iJES=0; iJES<NJES; iJES++)
        nIsojet[iJES]=isojets[iJES].size();
      
      // Load branches for the cuts
      bLoad(b["ZBosonPt"],ientry);
      bLoad(b["ZBosonPhi"],ientry);
      bLoad(b["hbbphi"],ientry);
      bLoad(b["fjPhi"],ientry);
      bLoad(b["hbbpt_reg"],ientry);
      bLoad(b["hbbm_reg"],ientry);
      bLoad(b["fjPt"],ientry);
      bLoad(b["fjMSD_corr"],ientry);
      
      // Set Selection Bits
      std::map<TString, bool> cut;
      for(unsigned iJES=0; iJES<NJES; iJES++) {
        if(iJES==0) {
          cut["ZpT"       ] = gt.ZBosonPt > 50;
          cut["ZpTFJ"     ] = gt.ZBosonPt > 250;
          cut["dPhiZHFJ"  ] = fabs(gt.fjPhi - gt.ZBosonPhi) > 2.5;
          cut["btag"      ] = bjet1IsTight && bjet2IsLoose;
          cut["bveto"     ] = !bjet1IsTight && !bjet2IsLoose;
          cut["btagFJ"    ] = gt.fjDoubleCSV > doubleBCut;
          cut["bvetoFJ"   ] = gt.fjDoubleCSV < doubleBCut;
          cut["Zmass"     ] = gt.ZBosonM >= 75 && gt.ZBosonM < 105;
          cut["ZmassTight"] = gt.ZBosonM >= 85 && gt.ZBosonM < 97;
          cut["ZmassSB"   ] = (gt.ZBosonM >= 10 && gt.ZBosonM < 75) || (gt.ZBosonM >= 105 && gt.ZBosonM < 120);
          cut["lowMET"    ] = gt.pfmet[0] < 60;
        } 
        if(iJES==(int)shiftjes::kJESTotalUp  ) cut["lowMET"] = gt.pfmet[1] < 60;
        if(iJES==(int)shiftjes::kJESTotalDown) cut["lowMET"] = gt.pfmet[2] < 60;
        cut["dPhiZH"  ] = fabs(gt.hbbphi[iJES] - gt.ZBosonPhi) > 2.5;
        cut["pTjj"    ] = gt.hbbpt_reg[iJES] > 100;
        cut["mjj"     ] = gt.hbbm_reg[iJES] >= 90 && gt.hbbm_reg[iJES] < 150;
        cut["mjjSB"   ] = !cut["mjjSB"] && gt.hbbm_reg[iJES]<250;
        cut["mSD"     ] = gt.fjMSD_corr[iJES] >= 40;
        cut["mSD_SR"  ] = gt.fjMSD_corr[iJES] >= 80 && gt.fjMSD_corr[iJES]<150;
        cut["mSD_SB"  ] = cut["mSD"] && gt.fjMSD_corr[iJES]<80;
        cut["pTFJ"    ] = gt.fjPt[iJES] > 250;
        cut["0ijb"    ] = isojetNBtags[iJES]==0;
        cut["1ijb"    ] = isojetNBtags[iJES]>0;
        
        selectionBits[iJES]=0; nMinusOneBits=0;
        for(unsigned iSel=0; iSel<selections.size(); iSel++) { 
          selectionType sel = selections[iSel];
          if(passAllCuts(cut, cuts[sel])) selectionBits[iJES] |= sel;
          if(iJES==0 && passNMinusOne(cut, cuts[sel]))
            nMinusOneBits |= sel;
        }
      
      }
 
    } // End Event Loop
  //} // End Chain Loop
  delete dummyTree; // no longer needed
}
