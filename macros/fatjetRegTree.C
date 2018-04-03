#include <map>
#include "PandaTree/Objects/interface/Event.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "vhbbPlot.h"
#include <sstream>
#include "PandaCore/Tools/interface/JERReader.h"
#include "TMVA/Reader.h"
#include "TMath.h"

using namespace panda;

const bool onlyNonPromptLeptons=true;

void fatjetRegTree(
  TString inputFileName,
  TString outputFileName,
  Long64_t maxEntries=-1,
  int debug=0
) {
  //vhbbPlot::sampleType sample = vhbbPlot::kVH,
  gSystem->Load("libPandaTreeObjects.so");
  gSystem->Load("libPandaCoreTools.so");
  
  // Apply MSD correction after the smearing
  TFile *MSDcorr = TFile::Open("PandaAnalysis/data/puppiCorr.root"); assert(MSDcorr);
  TF1 *puppisd_corrGEN = (TF1*)MSDcorr->Get("puppiJECcorr_gen");
  TF1 *puppisd_corrRECO_cen = (TF1*)MSDcorr->Get("puppiJECcorr_reco_0eta1v3");
  TF1 *puppisd_corrRECO_for = (TF1*)MSDcorr->Get("puppiJECcorr_reco_1v3eta2v5");
  JERReader *ak8JERReader = new JERReader("PandaAnalysis/data/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt","PandaAnalysis/data/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");
  
  TFile *inputFile=0;
  int retries=0;
  while(true) {
    inputFile = TFile::Open(inputFileName,"read");
    if(inputFile && inputFile->IsOpen()) break;
    retries++;
    if(retries>100) { throw std::runtime_error("Error opening input file"); return; }
  }
  TTree* tree = (TTree*)inputFile->Get("events"); // get the tree object from the file
  Event event; // create an Event object
  event.setStatus(*tree, {"!*"});
  event.setAddress(*tree, {
    "runNumber", "lumiNumber", "eventNumber", "rho", "weight",
    "muons", "electrons", "taus",
    "genParticles",
    "puppiAK8Jets", "chsAK4Jets", "ak8GenJets", 
    "pfCandidates","vertices","tracks","secondaryVertices"
  }); 

  Long64_t sum_mc_weights=0;  
  TH1D* all_tree=(TH1D*)inputFile->FindObjectAny("hSumW");
  sum_mc_weights = all_tree->GetBinContent(1);
  
 
  TFile *outputFile = new TFile(outputFileName, "recreate");
  TTree *regTree = new TTree("regTree","regTree");
  // Notes for myself:
  // fj1_emFraction, 
  
  float weight                  ; regTree->Branch("weight"                , &weight                  , "weight/F"                  );
  float fj1_pt                  ; regTree->Branch("fj1_pt"                , &fj1_pt                  , "fj1_pt/F"                  );
  float fj1_eta                 ; regTree->Branch("fj1_eta"               , &fj1_eta                 , "fj1_eta/F"                 );
  float fj1_phi                 ; regTree->Branch("fj1_phi"               , &fj1_phi                 , "fj1_phi/F"                 );
  float fj1_doubleB             ; regTree->Branch("fj1_doubleB"           , &fj1_doubleB             , "fj1_doubleB/F"             );
  float fj1_MSD                 ; regTree->Branch("fj1_MSD"               , &fj1_MSD                 , "fj1_MSD/F"                 );
  float fj1_Mnu                 ; regTree->Branch("fj1_Mnu"               , &fj1_Mnu                 , "fj1_Mnu/F"                 );
  float fj1_genM                ; regTree->Branch("fj1_genM"              , &fj1_genM                , "fj1_genM/F"                );
  float fj1_genMnu              ; regTree->Branch("fj1_genMnu"            , &fj1_genMnu              , "fj1_genMnu/F"              );
  float fj1_nhf                 ; regTree->Branch("fj1_nhf"               , &fj1_nhf                 , "fj1_nhf/F"                 );
  float fj1_chf                 ; regTree->Branch("fj1_chf"               , &fj1_chf                 , "fj1_chf/F"                 );
  float fj1_nef                 ; regTree->Branch("fj1_nef"               , &fj1_nef                 , "fj1_nef/F"                 );
  float fj1_cef                 ; regTree->Branch("fj1_cef"               , &fj1_cef                 , "fj1_cef/F"                 );
  float fj1_emFraction          ; regTree->Branch("fj1_emFraction"        , &fj1_emFraction          , "fj1_emFraction/F"          );
  int   fj1_nLep                ; regTree->Branch("fj1_nLep"              , &fj1_nLep                , "fj1_nLep/I"                );
  float sumNu_pt                ; regTree->Branch("sumNu_pt"              , &sumNu_pt                , "sumNu_pt/F"                );
  float sumNu_dEtaFj1           ; regTree->Branch("sumNu_dEtaFj1"         , &sumNu_dEtaFj1           , "sumNu_dEtaFj1/F"           );
  float sumNu_dPhiFj1           ; regTree->Branch("sumNu_dPhiFj1"         , &sumNu_dPhiFj1           , "sumNu_dPhiFj1/F"           );
  float sumNu_E                 ; regTree->Branch("sumNu_E"               , &sumNu_E                 , "sumNu_E/F"                 );
  float ele1_pt                 ; regTree->Branch("ele1_pt"               , &ele1_pt                 , "ele1_pt/F"                 );
  float ele1_ptRelFj1           ; regTree->Branch("ele1_ptRelFj1"         , &ele1_ptRelFj1           , "ele1_ptRelFj1/F"           );
  float ele1_dRFj1              ; regTree->Branch("ele1_dRFj1"            , &ele1_dRFj1              , "ele1_dRFj1/F"              );
  float ele1_eta                ; regTree->Branch("ele1_eta"              , &ele1_eta                , "ele1_eta/F"                );
  float ele1_phi                ; regTree->Branch("ele1_phi"              , &ele1_phi                , "ele1_phi/F"                );
  float ele1_dxy                ; regTree->Branch("ele1_dxy"              , &ele1_dxy                , "ele1_dxy/F"                );
  float ele1_dz                 ; regTree->Branch("ele1_dz"               , &ele1_dz                 , "ele1_dz/F"                 );
  float ele1_hfMva              ; regTree->Branch("ele1_hfMva"            , &ele1_hfMva              , "ele1_hfMva/F"              );
  float ele1_sv1_pt             ; regTree->Branch("ele1_sv1_pt"           , &ele1_sv1_pt             , "ele1_sv1_pt/F"             ); 
  float ele1_sv1_ptRelFj1       ; regTree->Branch("ele1_sv1_ptRelFj1"     , &ele1_sv1_ptRelFj1       , "ele1_sv1_ptRelFj1/F"       ); 
  float ele1_sv1_eta            ; regTree->Branch("ele1_sv1_eta"          , &ele1_sv1_eta            , "ele1_sv1_eta/F"            );  
  float ele1_sv1_phi            ; regTree->Branch("ele1_sv1_phi"          , &ele1_sv1_phi            , "ele1_sv1_phi/F"            );  
  float ele1_sv1_dRFj1          ; regTree->Branch("ele1_sv1_dRFj1"        , &ele1_sv1_dRFj1          , "ele1_sv1_dRFj1/F"          );  
  float ele1_sv1_m              ; regTree->Branch("ele1_sv1_m"            , &ele1_sv1_m              , "ele1_sv1_m/F"              ); 
  float ele1_sv1_logMRel        ; regTree->Branch("ele1_sv1_logMRel"      , &ele1_sv1_logMRel        , "ele1_sv1_logMRel/F"        ); 
  float ele1_sv1_vtx3DVal       ; regTree->Branch("ele1_sv1_vtx3DVal"     , &ele1_sv1_vtx3DVal       , "ele1_sv1_vtx3DVal/F"       ); 
  float ele1_sv1_vtx3DeVal      ; regTree->Branch("ele1_sv1_vtx3DeVal"    , &ele1_sv1_vtx3DeVal      , "ele1_sv1_vtx3DeVal/F"      ); 
  int   ele1_sv1_ntrk           ; regTree->Branch("ele1_sv1_ntrk"         , &ele1_sv1_ntrk           , "ele1_sv1_ntrk/I"           ); 
  float ele1_sv1_chi2           ; regTree->Branch("ele1_sv1_chi2"         , &ele1_sv1_chi2           , "ele1_sv1_chi2/F"           );  
  float ele1_sv1_significance   ; regTree->Branch("ele1_sv1_significance" , &ele1_sv1_significance   , "ele1_sv1_significance/F"   );  
  float mu1_pt                  ; regTree->Branch("mu1_pt"                , &mu1_pt                  , "mu1_pt/F"                  );
  float mu1_ptRelFj1            ; regTree->Branch("mu1_ptRelFj1"          , &mu1_ptRelFj1            , "mu1_ptRelFj1/F"            );
  float mu1_dRFj1               ; regTree->Branch("mu1_dRFj1"             , &mu1_dRFj1               , "mu1_dRFj1/F"               );
  float mu1_eta                 ; regTree->Branch("mu1_eta"               , &mu1_eta                 , "mu1_eta/F"                 );
  float mu1_phi                 ; regTree->Branch("mu1_phi"               , &mu1_phi                 , "mu1_phi/F"                 );
  float mu1_dxy                 ; regTree->Branch("mu1_dxy"               , &mu1_dxy                 , "mu1_dxy/F"                 );
  float mu1_dz                  ; regTree->Branch("mu1_dz"                , &mu1_dz                  , "mu1_dz/F"                  );
  float mu1_sv1_pt              ; regTree->Branch("mu1_sv1_pt"            , &mu1_sv1_pt              , "mu1_sv1_pt/F"              ); 
  float mu1_sv1_ptRelFj1        ; regTree->Branch("mu1_sv1_ptRelFj1"      , &mu1_sv1_ptRelFj1        , "mu1_sv1_ptRelFj1/F"        ); 
  float mu1_sv1_dRFj1           ; regTree->Branch("mu1_sv1_dRFj1"         , &mu1_sv1_dRFj1           , "mu1_sv1_dRFj1/F"           ); 
  float mu1_sv1_eta             ; regTree->Branch("mu1_sv1_eta"           , &mu1_sv1_eta             , "mu1_sv1_eta/F"             );  
  float mu1_sv1_phi             ; regTree->Branch("mu1_sv1_phi"           , &mu1_sv1_phi             , "mu1_sv1_phi/F"             );  
  float mu1_sv1_m               ; regTree->Branch("mu1_sv1_m"             , &mu1_sv1_m               , "mu1_sv1_m/F"               ); 
  float mu1_sv1_logMRel         ; regTree->Branch("mu1_sv1_logMRel"       , &mu1_sv1_logMRel         , "mu1_sv1_logMRel/F"         ); 
  float mu1_sv1_vtx3DVal        ; regTree->Branch("mu1_sv1_vtx3DVal"      , &mu1_sv1_vtx3DVal        , "mu1_sv1_vtx3DVal/F"        ); 
  float mu1_sv1_vtx3DeVal       ; regTree->Branch("mu1_sv1_vtx3DeVal"     , &mu1_sv1_vtx3DeVal       , "mu1_sv1_vtx3DeVal/F"       ); 
  int   mu1_sv1_ntrk            ; regTree->Branch("mu1_sv1_ntrk"          , &mu1_sv1_ntrk            , "mu1_sv1_ntrk/I"            ); 
  float mu1_sv1_chi2            ; regTree->Branch("mu1_sv1_chi2"          , &mu1_sv1_chi2            , "mu1_sv1_chi2/F"            );  
  float mu1_sv1_significance    ; regTree->Branch("mu1_sv1_significance"  , &mu1_sv1_significance    , "mu1_sv1_significance/F"    );  
  float fj1_sv1_pt              ; regTree->Branch("fj1_sv1_pt"            , &fj1_sv1_pt              , "fj1_sv1_pt/F"              ); 
  float fj1_sv1_ptRelFj1        ; regTree->Branch("fj1_sv1_ptRelFj1"      , &fj1_sv1_ptRelFj1        , "fj1_sv1_ptRelFj1/F"        ); 
  float fj1_sv1_dRFj1           ; regTree->Branch("fj1_sv1_dRFj1"         , &fj1_sv1_dRFj1           , "fj1_sv1_dRFj1/F"           ); 
  float fj1_sv1_eta             ; regTree->Branch("fj1_sv1_eta"           , &fj1_sv1_eta             , "fj1_sv1_eta/F"            );  
  float fj1_sv1_phi             ; regTree->Branch("fj1_sv1_phi"           , &fj1_sv1_phi             , "fj1_sv1_phi/F"            );  
  float fj1_sv1_m               ; regTree->Branch("fj1_sv1_m"             , &fj1_sv1_m               , "fj1_sv1_m/F"               ); 
  float fj1_sv1_logMRel         ; regTree->Branch("fj1_sv1_logMRel"       , &fj1_sv1_logMRel         , "fj1_sv1_logMRel/F"         ); 
  float fj1_sv1_vtx3DVal        ; regTree->Branch("fj1_sv1_vtx3DVal"      , &fj1_sv1_vtx3DVal        , "fj1_sv1_vtx3DVal/F"        ); 
  float fj1_sv1_vtx3DeVal       ; regTree->Branch("fj1_sv1_vtx3DeVal"     , &fj1_sv1_vtx3DeVal       , "fj1_sv1_vtx3DeVal/F"       ); 
  int   fj1_sv1_ntrk            ; regTree->Branch("fj1_sv1_ntrk"          , &fj1_sv1_ntrk            , "fj1_sv1_ntrk/I"            ); 
  float fj1_sv1_chi2            ; regTree->Branch("fj1_sv1_chi2"          , &fj1_sv1_chi2            , "fj1_sv1_chi2/F"            );  
  float fj1_sv1_significance    ; regTree->Branch("fj1_sv1_significance"  , &fj1_sv1_significance    , "fj1_sv1_significance/F"    );  
  float fj1_sv2_pt              ; regTree->Branch("fj1_sv2_pt"            , &fj1_sv2_pt              , "fj1_sv2_pt/F"              ); 
  float fj1_sv2_ptRelFj1        ; regTree->Branch("fj1_sv2_ptRelFj1"      , &fj1_sv2_ptRelFj1        , "fj1_sv2_ptRelFj1/F"        ); 
  float fj1_sv2_dRFj1           ; regTree->Branch("fj1_sv2_dRFj1"         , &fj1_sv2_dRFj1           , "fj1_sv2_dRFj1/F"           ); 
  float fj1_sv2_eta             ; regTree->Branch("fj1_sv2_eta"           , &fj1_sv2_eta             , "fj1_sv2_eta/F"            );  
  float fj1_sv2_phi             ; regTree->Branch("fj1_sv2_phi"           , &fj1_sv2_phi             , "fj1_sv2_phi/F"            );  
  float fj1_sv2_m               ; regTree->Branch("fj1_sv2_m"             , &fj1_sv2_m               , "fj1_sv2_m/F"               ); 
  float fj1_sv2_logMRel         ; regTree->Branch("fj1_sv2_logMRel"       , &fj1_sv2_logMRel         , "fj1_sv2_logMRel/F"         ); 
  float fj1_sv2_vtx3DVal        ; regTree->Branch("fj1_sv2_vtx3DVal"      , &fj1_sv2_vtx3DVal        , "fj1_sv2_vtx3DVal/F"        ); 
  float fj1_sv2_vtx3DeVal       ; regTree->Branch("fj1_sv2_vtx3DeVal"     , &fj1_sv2_vtx3DeVal       , "fj1_sv2_vtx3DeVal/F"       ); 
  int   fj1_sv2_ntrk            ; regTree->Branch("fj1_sv2_ntrk"          , &fj1_sv2_ntrk            , "fj1_sv2_ntrk/I"            ); 
  float fj1_sv2_chi2            ; regTree->Branch("fj1_sv2_chi2"          , &fj1_sv2_chi2            , "fj1_sv2_chi2/F"            );  
  float fj1_sv2_significance    ; regTree->Branch("fj1_sv2_significance"  , &fj1_sv2_significance    , "fj1_sv2_significance/F"    );  
 
  // Objects we don't want to reallocate over and over
  TLorentzVector sumNuV4, fatjetV4, nuV4, fatjetNuV4, genFatjetV4, genFatjetNuV4;
  FatJet *fj1=0;
  // Instantiate the HF Electron MVA Reader
  TMVA::Reader *hfEleMvaReader=new TMVA::Reader();
  float hfEleMvaInputs[18];
  hfEleMvaReader->AddVariable("ele_eta"            , &hfEleMvaInputs[ 0]);
  hfEleMvaReader->AddVariable("ele_ecalIso/ele_pt" , &hfEleMvaInputs[ 1]);
  hfEleMvaReader->AddVariable("ele_hcalIso/ele_pt" , &hfEleMvaInputs[ 2]);
  hfEleMvaReader->AddVariable("ele_trackIso/ele_pt", &hfEleMvaInputs[ 3]);
  hfEleMvaReader->AddVariable("ele_chIso/ele_pt"   , &hfEleMvaInputs[ 4]);
  hfEleMvaReader->AddVariable("ele_nhIso/ele_pt"   , &hfEleMvaInputs[ 5]);
  hfEleMvaReader->AddVariable("ele_phIso/ele_pt"   , &hfEleMvaInputs[ 6]);
  hfEleMvaReader->AddVariable("ele_dxy"            , &hfEleMvaInputs[ 7]);
  hfEleMvaReader->AddVariable("ele_dz"             , &hfEleMvaInputs[ 8]);
  hfEleMvaReader->AddVariable("ele_sieie"          , &hfEleMvaInputs[ 9]);
  hfEleMvaReader->AddVariable("ele_sipip"          , &hfEleMvaInputs[10]);
  hfEleMvaReader->AddVariable("ele_r9"             , &hfEleMvaInputs[11]);
  hfEleMvaReader->AddVariable("ele_dEtaInSeed"     , &hfEleMvaInputs[12]);
  hfEleMvaReader->AddVariable("ele_dPhiIn"         , &hfEleMvaInputs[13]);
  hfEleMvaReader->AddVariable("ele_hOverE"         , &hfEleMvaInputs[14]);
  hfEleMvaReader->AddVariable("ele_ooEmooP"        , &hfEleMvaInputs[15]);
  hfEleMvaReader->AddVariable("ele_nMissingHits"   , &hfEleMvaInputs[16]);
  hfEleMvaReader->AddVariable("ele_tripleCharge"   , &hfEleMvaInputs[17]);
  hfEleMvaReader->BookMVA("BDT","weights/hfEleMva_hfEleMva_hp0_test4.weights.xml");

  Long64_t nEntries = tree->GetEntries();
  if(maxEntries>0) nEntries=TMath::Min(maxEntries,nEntries);
  Long64_t oneTenth = nEntries/10;
  for (Long64_t iEntry = 0; iEntry<nEntries; ++iEntry) {
    if(debug || iEntry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",iEntry,nEntries); 
    //printf("hello %d\n",__LINE__);
    event.getEntry(*tree, iEntry);
    // Reset output tree branches
    sumNuV4.SetPtEtaPhiM(0,0,0,0);
    fj1=0;
    weight                =    0;
    fj1_pt                =   -1;
    fj1_eta               =   -9;
    fj1_phi               =   -9;
    fj1_doubleB           =   -9;
    fj1_MSD               =   -1;
    fj1_Mnu               =   -1;
    fj1_genM              =   -1;
    fj1_genMnu            =   -1;
    fj1_nhf               =   -1;
    fj1_chf               =   -1;
    fj1_nef               =   -1;
    fj1_cef               =   -1;
    sumNu_pt              =   -1;
    sumNu_dEtaFj1         =   -9;
    sumNu_dPhiFj1         =   -9;
    sumNu_E               =   -1;
    ele1_pt               =   -1;
    ele1_ptRelFj1         =   -1;
    ele1_dRFj1            =   -1;
    ele1_eta              =   -9;
    ele1_phi              =   -9;
    ele1_dxy              =  -99;
    ele1_dz               =  -99;
    ele1_hfMva            =   -1;
    ele1_sv1_pt           =   -1;
    ele1_sv1_ptRelFj1     =   -1;
    ele1_sv1_eta          =   -9;
    ele1_sv1_phi          =   -9;
    ele1_sv1_dRFj1        =   -1;
    ele1_sv1_m            =   -1;
    ele1_sv1_logMRel      =   -9;
    ele1_sv1_vtx3DVal     =   -1;
    ele1_sv1_vtx3DeVal    =   -1;
    ele1_sv1_ntrk         =   -1;
    ele1_sv1_chi2         =   -1;
    ele1_sv1_significance =   -1;
    mu1_pt                =   -1;
    mu1_ptRelFj1          =   -1;
    mu1_dRFj1             =   -1;
    mu1_eta               =   -9;
    mu1_phi               =   -9;
    mu1_dxy               =  -99;
    mu1_dz                =  -99;
    mu1_sv1_pt            =   -1;
    mu1_sv1_ptRelFj1      =   -1;
    mu1_sv1_dRFj1         =   -1;
    mu1_sv1_eta           =   -9;
    mu1_sv1_phi           =   -9;
    mu1_sv1_m             =   -1;
    mu1_sv1_logMRel       =  -99;
    mu1_sv1_vtx3DVal      =   -1;
    mu1_sv1_vtx3DeVal     =   -1;
    mu1_sv1_ntrk          =   -1;
    mu1_sv1_chi2          =   -1;
    mu1_sv1_significance  =   -1;
    fj1_sv1_pt            =   -1;
    fj1_sv1_ptRelFj1      =   -1;
    fj1_sv1_dRFj1         =   -1;
    fj1_sv1_eta           =   -9;
    fj1_sv1_phi           =   -9;
    fj1_sv1_m             =   -1;
    fj1_sv1_logMRel       =  -99;
    fj1_sv1_vtx3DVal      =   -1;
    fj1_sv1_vtx3DeVal     =   -1;
    fj1_sv1_ntrk          =   -1;
    fj1_sv1_chi2          =   -1;
    fj1_sv1_significance  =   -1;
    fj1_sv2_pt            =   -1;
    fj1_sv2_ptRelFj1      =   -1;
    fj1_sv2_dRFj1         =   -1;
    fj1_sv2_eta           =   -9;
    fj1_sv2_phi           =   -9;
    fj1_sv2_m             =   -1;
    fj1_sv2_logMRel       =  -99;
    fj1_sv2_vtx3DVal      =   -1;
    fj1_sv2_vtx3DeVal     =   -1;
    fj1_sv2_ntrk          =   -1;
    fj1_sv2_chi2          =   -1;
    fj1_sv2_significance  =   -1;
    
    for (unsigned nJ = 0; nJ<event.puppiAK8Jets.size(); nJ++){
      panda::FatJet* fatjet = &event.puppiAK8Jets[nJ];
      // Fatjet preselection
      if(fatjet->pt()<200) continue;
      if(fabs(fatjet->eta())>2.4) continue;
      if(!fatjet->loose) continue;

      // MSD correction 
      float genCorr  = 1.; 
      float recoCorr = 1.; 
      float totalCorr = 1.; 
      genCorr = puppisd_corrGEN->Eval( fatjet->pt() );
      if ( fabs(fatjet->eta()) <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( fatjet->pt() );
      else recoCorr = puppisd_corrRECO_for->Eval( fatjet->pt() );
      totalCorr = genCorr * recoCorr;
      float fj1_MSD_unsmeared = fatjet->mSD * totalCorr;
      // Smearing
      double smear=1, smearUp=1, smearDown=1;
      ak8JERReader->getStochasticSmear(fatjet->pt(),fatjet->eta(),event.rho,smear,smearUp,smearDown);

      // Populate output tree info
      if(!fj1) {
        fj1_pt      = fatjet->pt() * smear; 
        fj1_eta     = fatjet->eta();
        fj1_phi     = fatjet->phi();
        fj1_doubleB = fatjet->double_sub;
        fj1_MSD     = fj1_MSD_unsmeared * smear;
        fj1_nhf     = fatjet->nhf;
        fj1_chf     = fatjet->chf; 
        fj1_nef     = fatjet->nef; 
        fj1_cef     = fatjet->cef; 
        fj1_emFraction = fatjet->nef + fatjet->cef;
        fj1=fatjet; // keep the pointer "fj1" in scope
        break;
      }
    }
    if(!fj1) continue;

    // Final state neutrinos near the fatjet
    // No parentage calculation for now
    unsigned char nB=0,nC=0;
    std::vector<const GenParticle*> validGenP;
    for (unsigned iG = 0; iG<event.genParticles.size(); iG++) {
      auto& genParticle = event.genParticles[iG];
      if (genParticle.finalState != 1) continue;
      unsigned absid = abs(genParticle.pdgid);
      if (absid!=12 && absid!=14 && absid!=16 && absid!=5) continue;
      if (fj1->dR2((GenParticle)genParticle)>0.64) continue;
      
      bool isDuplicate=false;
      for (auto* vgp : validGenP)
        if(genParticle.pdgid==vgp->pdgid &&
          fabs(genParticle.pt()/vgp->pt()-1)<0.10 &&
          sqrt(pow(genParticle.eta()-vgp->eta(),2)+pow(TVector2::Phi_mpi_pi(genParticle.phi()-vgp->phi()),2)) < 0.00001
        ) { isDuplicate=true; break; }
      if(isDuplicate) continue;
      validGenP.push_back(&genParticle);
      
      short parentId=0; if(genParticle.parent.isValid()) parentId = int(genParticle.parent.get()->pdgid);
      if(absid==12 || absid==14 || absid==16) {
        // Make sure neutrino is from a HF hadron
        // Allowed codes: 400-499 (charmed meson), 500-599 (B meson),
        // above 4000-5999 (Lambda_C^+ or something), and some weird D/B mesons over 10000
        unsigned short absParentIdMod10k = unsigned(abs(parentId))%10000;
        bool nuFromHFHadron = (absParentIdMod10k>=400 && absParentIdMod10k<600) || (absParentIdMod10k>=4000 && absParentIdMod10k<6000);
        if(absParentIdMod10k==15) { // go up the chain for taus
          bool tauFromHFHadron=false;
          GenParticle *familyMember = (GenParticle*)genParticle.parent.get();
          GenParticle *ancestor=0;
          while(!tauFromHFHadron) {
            if(!familyMember->parent.isValid()) break;
            ancestor = (GenParticle*) familyMember->parent.get();
            unsigned short absAncestorIdMod10k = abs(ancestor->pdgid)%10000;
            if(
              (absAncestorIdMod10k>=400 && absAncestorIdMod10k<600) || 
              (absAncestorIdMod10k>=4000 && absAncestorIdMod10k<6000)
            ) {  tauFromHFHadron=true; break;  }
            else familyMember=ancestor;
          }
          if(tauFromHFHadron) nuFromHFHadron = true; 
        }
        if(!nuFromHFHadron) continue;
        nuV4.SetPtEtaPhiM(genParticle.pt(), genParticle.eta(), genParticle.phi(), 0);
        sumNuV4+=nuV4;
      } else if(absid==5) { // B quark counting, this is useless right now
        if(genParticle.parent.isValid() && genParticle.parent->pdgid==genParticle.pdgid) continue;
        nB++; 
      }
    }
    fatjetV4.SetPtEtaPhiM(fj1_pt, fj1_eta, fj1_phi, fj1_MSD);
    fatjetNuV4 = sumNuV4 + fatjetV4;
    fj1_Mnu = fatjetNuV4.M();
    if(fj1->matchedGenJet.isValid()) {
      GenJet *theGenJet = (GenJet*)fj1->matchedGenJet.get();
      genFatjetV4.SetPtEtaPhiM(theGenJet->pt(), theGenJet->eta(), theGenJet->phi(), theGenJet->m());
      genFatjetNuV4 = genFatjetV4 + sumNuV4;
      fj1_genM = theGenJet->m();
      fj1_genMnu = genFatjetNuV4.M();
    } 
    sumNu_pt                = sumNuV4.Pt() ;
    sumNu_dEtaFj1           = fabs(fj1->eta()-sumNuV4.Eta());
    sumNu_dPhiFj1           = fabs(TVector2::Phi_mpi_pi(fj1->phi()-sumNuV4.Phi()));
    sumNu_E                 = sumNuV4.E()  ;

    // Properties of the leading electrons in the fatjet
    std::vector<Electron*> hfElectrons;
    std::vector<float> hfEleMvaScores;
    for(unsigned iE = 0; iE<event.electrons.size(); iE++) {
      auto& electron  = event.electrons[iE];
      if(electron.pt()<5) break;
      if(fabs(electron.eta()) > 2.5) continue;
      if(!electron.conversionVeto) continue;
      if(fj1->dR2(electron)>0.64) continue;

      // Check promptness
      bool isPrompt=true;
      if(onlyNonPromptLeptons && electron.matchedPF.isValid()) for(unsigned iSV=0; iSV<event.secondaryVertices.size(); iSV++) {
        auto& sv = event.secondaryVertices[iSV];
        if (sv.chi2 <= 0 || sv.ndof<=0 || sv.ndof!=sv.ndof) continue;
        for (UShort_t iD=0; iD<sv.daughters.size(); iD++) {
          if(!sv.daughters.at(iD).isValid()) continue;
          if(sv.daughters.at(iD).get()==electron.matchedPF.get()) {
            isPrompt=false;
            break;
          }
        }
      }
      if(onlyNonPromptLeptons && isPrompt) continue;
      
      // Check MVA score
      float trackPt = electron.trackP/TMath::CosH(electron.eta());
      float ooEmooP=fabs(1.-electron.eseed/max(.0001f,electron.trackP))/electron.ecalE;
      hfEleMvaInputs[ 0] =electron.eta()           ;
      hfEleMvaInputs[ 1] =electron.ecalIso/trackPt ;
      hfEleMvaInputs[ 2] =electron.hcalIso/trackPt ;
      hfEleMvaInputs[ 3] =electron.trackIso/trackPt;
      hfEleMvaInputs[ 4] =electron.chIso/trackPt   ;
      hfEleMvaInputs[ 5] =electron.nhIso/trackPt   ;
      hfEleMvaInputs[ 6] =electron.phIso/trackPt   ;
      hfEleMvaInputs[ 7] =electron.dxy             ;
      hfEleMvaInputs[ 8] =electron.dz              ;
      hfEleMvaInputs[ 9] =electron.sieie           ;
      hfEleMvaInputs[10] =electron.sipip           ;
      hfEleMvaInputs[11] =electron.r9              ;
      hfEleMvaInputs[12] =electron.dEtaInSeed      ;
      hfEleMvaInputs[13] =electron.dPhiIn          ;
      hfEleMvaInputs[14] =electron.hOverE          ;
      hfEleMvaInputs[15] =ooEmooP                  ;
      hfEleMvaInputs[16] =(float)electron.nMissingHits;
      hfEleMvaInputs[17] =(float)electron.tripleCharge;
      float hfEleMvaScore = hfEleMvaReader->EvaluateMVA("BDT");
      if(debug>=3) { printf("  Electron %d, HF MVA score = %.3f\n",iE,hfEleMvaScore);
        printf("  hfEleMvaInputs[ 0] = electron.eta()            = %.2e\n",electron.eta()           );
        printf("  hfEleMvaInputs[ 1] = electron.ecalIso/trackPt  = %.2e\n",electron.ecalIso/trackPt );
        printf("  hfEleMvaInputs[ 2] = electron.hcalIso/trackPt  = %.2e\n",electron.hcalIso/trackPt );
        printf("  hfEleMvaInputs[ 3] = electron.trackIso/trackPt = %.2e\n",electron.trackIso/trackPt);
        printf("  hfEleMvaInputs[ 4] = electron.chIso/trackPt    = %.2e\n",electron.chIso/trackPt   );
        printf("  hfEleMvaInputs[ 5] = electron.nhIso/trackPt    = %.2e\n",electron.nhIso/trackPt   );
        printf("  hfEleMvaInputs[ 6] = electron.phIso/trackPt    = %.2e\n",electron.phIso/trackPt   );
        printf("  hfEleMvaInputs[ 7] = electron.dxy              = %.2e\n",electron.dxy             );
        printf("  hfEleMvaInputs[ 8] = electron.dz               = %.2e\n",electron.dz              );
        printf("  hfEleMvaInputs[ 9] = electron.sieie            = %.2e\n",electron.sieie           );
        printf("  hfEleMvaInputs[10] = electron.sipip            = %.2e\n",electron.sipip           );
        printf("  hfEleMvaInputs[11] = electron.r9               = %.2e\n",electron.r9              );
        printf("  hfEleMvaInputs[12] = electron.dEtaInSeed       = %.2e\n",electron.dEtaInSeed      );
        printf("  hfEleMvaInputs[13] = electron.dPhiIn           = %.2e\n",electron.dPhiIn          );
        printf("  hfEleMvaInputs[14] = electron.hOverE           = %.2e\n",electron.hOverE          );
        printf("  hfEleMvaInputs[15] = ooEmooP                   = %.2e\n",ooEmooP                  );
        printf("  hfEleMvaInputs[16] = electron.nMissingHits     = %.2e\n",(float)electron.nMissingHits    );
        printf("  hfEleMvaInputs[17] = electron.tripleCharge     = %.2e\n",(float)electron.tripleCharge    );
      }

      if(hfEleMvaScore<0) continue;
      hfElectrons.push_back(&electron);
      hfEleMvaScores.push_back(hfEleMvaScore);

    }

    // Properties of the leading muons in the fatjet
    std::vector<Muon*> hfMuons;
    for(unsigned iM = 0; iM<event.muons.size(); iM++){
      auto& muon  = event.muons[iM];
      if(muon.pt()<5) break;
      if(!muon.soft) continue;
      if(fabs(muon.eta())>2.4) continue;
      if(fj1->dR2(muon)>0.64) continue;
      // Check promptness
      bool isPrompt=true;
      if(onlyNonPromptLeptons && muon.matchedPF.isValid()) for(unsigned iSV=0; iSV<event.secondaryVertices.size(); iSV++) {
        auto& sv = event.secondaryVertices[iSV];
        if (sv.chi2 <= 0 || sv.ndof<=0 || sv.ndof!=sv.ndof) continue;
        for (UShort_t iD=0; iD<sv.daughters.size(); iD++) {
          if(!sv.daughters.at(iD).isValid()) continue;
          if(sv.daughters.at(iD).get()==muon.matchedPF.get()) {
            isPrompt=false;
            break;
          }
        }
      }
      if(onlyNonPromptLeptons && isPrompt) continue;
      hfMuons.push_back(&muon);
    }
    // See if there are secondary vertices who are the matriarch of these leptons, and take the hardest matriarch of each one
    // Initialize vectors of pointers to SVs with null pointers
    std::vector<SecondaryVertex*> hfMuonSVs(hfMuons.size());
    std::vector<SecondaryVertex*> hfElectronSVs(hfElectrons.size());
    // Looping over all secondary vertices
    std::vector<SecondaryVertex *> fj1_SVs;
    for(unsigned iSV=0; iSV<event.secondaryVertices.size(); iSV++) {
      auto& sv = event.secondaryVertices[iSV];
      if (sv.chi2 <= 0 || sv.ndof<=0 || sv.ndof!=sv.ndof) continue;
      // Loop over this SV's daughters
      for (UShort_t iD=0; iD<sv.daughters.size(); iD++) {
        if(!sv.daughters.at(iD).isValid()) continue;
        // See if this SV is the matriarch of any of our HF muons
        for(unsigned iHFM=0; iHFM<hfMuons.size(); iHFM++) {
          if(!hfMuons[iHFM]->matchedPF.isValid()) continue;
          if(sv.daughters.at(iD).get()!=hfMuons[iHFM]->matchedPF.get()) continue;
          // Save this SV for this muon if we don't already have it, or replace the existing one
          // if it's harder.
          float bestPt = hfMuonSVs[iHFM]? hfMuonSVs[iHFM]->pt():0;
          if(sv.pt() > bestPt) hfMuonSVs[iHFM] = &sv;
        }
        // Same calculation for the electrons
        for(unsigned iHFE=0; iHFE<hfElectrons.size(); iHFE++) {
          if(!hfElectrons[iHFE]->matchedPF.isValid()) continue;
          if(sv.daughters.at(iD).get()!=hfElectrons[iHFE]->matchedPF.get()) continue;
          float bestPt = hfElectronSVs[iHFE]? hfElectronSVs[iHFE]->pt():0;
          if(sv.pt() > bestPt) hfElectronSVs[iHFE] = &sv;
        }
        // Save all the secondary vertices for the fatjet and sort by significance later
        for(UShort_t iFJC=0; iFJC<fj1->constituents.size(); iFJC++) {
          if(!fj1->constituents.at(iFJC).isValid()) continue;
          if(sv.daughters.at(iD).get()!=fj1->constituents.at(iFJC).get()) continue;
          if(std::find(fj1_SVs.begin(), fj1_SVs.end(), &sv) == fj1_SVs.end()) fj1_SVs.push_back(&sv);
        }
      }
    }
    std::sort(fj1_SVs.begin(), fj1_SVs.end(), [](const SecondaryVertex *sv1, const SecondaryVertex *sv2) -> bool { 
      //return sv1->significance > sv2->significance; 
      return sv1->pt() > sv2->pt(); 
    });

    fj1_nLep = hfElectrons.size() + hfMuons.size();

    if(hfElectrons.size()>0) {
      ele1_pt        =hfElectrons[0]->pt() ;
      ele1_ptRelFj1  =hfElectrons[0]->pt()/fj1->pt();
      ele1_dRFj1     =hfElectrons[0]->dR(*fj1);
      ele1_eta       =hfElectrons[0]->eta();
      ele1_phi       =hfElectrons[0]->phi();
      ele1_dxy       =hfElectrons[0]->dxy  ;
      ele1_dz        =hfElectrons[0]->dz   ;
      ele1_hfMva     =hfEleMvaScores[0]    ;
      if(hfElectronSVs[0]) {
        ele1_sv1_pt           = hfElectronSVs[0]->pt()         ;
        ele1_sv1_ptRelFj1     = hfElectronSVs[0]->pt()/fj1->pt();
        ele1_sv1_dRFj1        = hfElectronSVs[0]->dR(*fj1)     ;
        ele1_sv1_eta          = hfElectronSVs[0]->eta()          ;
        ele1_sv1_phi          = hfElectronSVs[0]->phi()          ;
        ele1_sv1_m            = hfElectronSVs[0]->m()          ;
        ele1_sv1_logMRel      = (hfElectronSVs[0]->m()>0 && fj1_MSD>0) ? TMath::Log10(hfElectronSVs[0]->m()/fj1_MSD) : -99;
        ele1_sv1_vtx3DVal     = hfElectronSVs[0]->vtx3DVal     ;
        ele1_sv1_vtx3DeVal    = hfElectronSVs[0]->vtx3DeVal    ;
        ele1_sv1_ntrk         = hfElectronSVs[0]->ntrk         ;
        ele1_sv1_chi2         = hfElectronSVs[0]->chi2         ;
        ele1_sv1_significance = hfElectronSVs[0]->significance ;
      }
    } if(hfMuons.size()>0) {
      mu1_pt        =hfMuons[0]->pt() ;
      mu1_ptRelFj1  =hfMuons[0]->pt()/fj1->pt();
      mu1_dRFj1     =hfMuons[0]->dR(*fj1);
      mu1_eta       =hfMuons[0]->eta();
      mu1_phi       =hfMuons[0]->phi();
      mu1_dxy       =hfMuons[0]->dxy  ;
      mu1_dz        =hfMuons[0]->dz   ;
      if(hfMuonSVs[0]) {
        mu1_sv1_pt           = hfMuonSVs[0]->pt()         ;
        mu1_sv1_ptRelFj1     = hfMuonSVs[0]->pt()/fj1->pt();
        mu1_sv1_dRFj1        = hfMuonSVs[0]->dR(*fj1)     ;
        mu1_sv1_eta          = hfMuonSVs[0]->eta()          ;
        mu1_sv1_phi          = hfMuonSVs[0]->phi()          ;
        mu1_sv1_m            = hfMuonSVs[0]->m()          ;
        mu1_sv1_logMRel      = (hfMuonSVs[0]->m()>0 && fj1_MSD>0) ? TMath::Log10(hfMuonSVs[0]->m()/fj1_MSD) : -99;
        mu1_sv1_vtx3DVal     = hfMuonSVs[0]->vtx3DVal     ;
        mu1_sv1_vtx3DeVal    = hfMuonSVs[0]->vtx3DeVal    ;
        mu1_sv1_ntrk         = hfMuonSVs[0]->ntrk         ;
        mu1_sv1_chi2         = hfMuonSVs[0]->chi2         ;
        mu1_sv1_significance = hfMuonSVs[0]->significance ;
      }
    }
    if(fj1_SVs.size()>0) {
      fj1_sv1_pt           = fj1_SVs[0]->pt()         ;
      fj1_sv1_ptRelFj1     = fj1_SVs[0]->pt()/fj1->pt();
      fj1_sv1_dRFj1        = fj1_SVs[0]->dR(*fj1)     ;
      fj1_sv1_eta          = fj1_SVs[0]->eta()          ;
      fj1_sv1_phi          = fj1_SVs[0]->phi()          ;
      fj1_sv1_m            = fj1_SVs[0]->m()          ;
      fj1_sv1_logMRel      = (fj1_SVs[0]->m()>0 && fj1_MSD>0) ? TMath::Log10(fj1_SVs[0]->m()/fj1_MSD) : -99;
      fj1_sv1_vtx3DVal     = fj1_SVs[0]->vtx3DVal     ;
      fj1_sv1_vtx3DeVal    = fj1_SVs[0]->vtx3DeVal    ;
      fj1_sv1_ntrk         = fj1_SVs[0]->ntrk         ;
      fj1_sv1_chi2         = fj1_SVs[0]->chi2         ;
      fj1_sv1_significance = fj1_SVs[0]->significance ;
    } if(fj1_SVs.size()>1) {
      fj1_sv2_pt           = fj1_SVs[1]->pt()         ;
      fj1_sv2_ptRelFj1     = fj1_SVs[1]->pt()/fj1->pt();
      fj1_sv2_dRFj1        = fj1_SVs[1]->dR(*fj1)     ;
      fj1_sv2_eta          = fj1_SVs[1]->eta()          ;
      fj1_sv2_phi          = fj1_SVs[1]->phi()          ;
      fj1_sv2_m            = fj1_SVs[1]->m()          ;
      fj1_sv2_logMRel      = (fj1_SVs[1]->m()>0 && fj1_MSD>0) ? TMath::Log10(fj1_SVs[1]->m()/fj1_MSD) : -99;
      fj1_sv2_vtx3DVal     = fj1_SVs[1]->vtx3DVal     ;
      fj1_sv2_vtx3DeVal    = fj1_SVs[1]->vtx3DeVal    ;
      fj1_sv2_ntrk         = fj1_SVs[1]->ntrk         ;
      fj1_sv2_chi2         = fj1_SVs[1]->chi2         ;
      fj1_sv2_significance = fj1_SVs[1]->significance ;
    }
    weight = event.weight;
    regTree->Fill();
  }
  regTree->Write("regTree",TObject::kOverwrite);
  TH1D *sum_MCW=new TH1D("sum_weights","sum_weights", 1, 0, 2);
  sum_MCW->SetBinContent(1,sum_mc_weights);
  sum_MCW->Write("sum_weights");
  outputFile->Close();
  inputFile->Close();
 
}

