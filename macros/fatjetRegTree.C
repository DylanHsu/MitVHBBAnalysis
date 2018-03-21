#include <map>
#include "PandaTree/Objects/interface/Event.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "vhbbPlot.h"
#include <sstream>
//#include "RoccoR.cc"
//#include <TRandom3.h>

void fatjetRegTree(
  TString inputFileName,
  TString outputFileName,
  vhbbPlot::sampleType = vhbbPlot::kVH,
  Long64_t maxEntries=-1
) {
  gSystem->Load("libPandaTreeObjects.so");
  TFile *inputFile = TFile::Open(inputFileName,"read");
  TTree* tree = (TTree*)inputFile->Get("events"); // get the tree object from the file
  panda::Event event; // create an Event object
  event.setStatus(*tree, {"!*"});
  event.setAddress(*tree, {
    "runNumber", "lumiNumber", "eventNumber",
    "muons", "electrons", "taus",
    "genParticles",
    "puppiAK8Jets", "chsAK4Jets",
    "pfCandidates","vertices","weight","tracks","secondaryVertices"
  }); 

  Long64_t sum_mc_weights=0;  
  TH1D* all_tree=(TH1D*)inputFile->FindObjectAny("hSumW");
  sum_mc_weights = all_tree->GetBinContent(1);
  

  TFile *outputFile = new TFile(outputFileName, "recreate");
  TTree *regTree = new TTree("regTree","regTree");
  float fatjet_pt      ; regTree->Branch("fatjet_pt"      , &fatjet_pt      , "fatjet_pt/F"     );
  float fatjet_eta     ; regTree->Branch("fatjet_eta"     , &fatjet_eta     , "fatjet_eta/F"    );
  float fatjet_phi     ; regTree->Branch("fatjet_phi"     , &fatjet_phi     , "fatjet_phi/F"    );
  float fatjet_MSD     ; regTree->Branch("fatjet_MSD"     , &fatjet_MSD     , "fatjet_MSD/F"    );
  float fatjet_Mnu     ; regTree->Branch("fatjet_Mnu"     , &fatjet_Mnu     , "fatjet_Mnu/F"    );
  float fatjet_nhf     ; regTree->Branch("fatjet_nhf"     , &fatjet_nhf     , "fatjet_nhf/F"    );
  float fatjet_chf     ; regTree->Branch("fatjet_chf"     , &fatjet_chf     , "fatjet_chf/F"    );
  float fatjet_nef     ; regTree->Branch("fatjet_nef"     , &fatjet_nef     , "fatjet_nef/F"    );
  float fatjet_cef     ; regTree->Branch("fatjet_cef"     , &fatjet_cef     , "fatjet_cef/F"    );
  
  float sumNu_pt       ; regTree->Branch("sumNu_pt"       , &sumNu_pt       , "sumNu_pt/F"     );
  float sumNu_eta      ; regTree->Branch("sumNu_eta"      , &sumNu_eta      , "sumNu_eta/F"    );
  float sumNu_phi      ; regTree->Branch("sumNu_phi"      , &sumNu_phi      , "sumNu_phi/F"    );
  float sumNu_E        ; regTree->Branch("sumNu_E"        , &sumNu_E        , "sumNu_E/F"      );
  
  float          ele1_pt                 ; regTree->Branch("ele1_pt"        , &ele1_pt        , "ele1_pt/F"      );
  float          ele1_eta                ; regTree->Branch("ele1_eta"       , &ele1_eta       , "ele1_eta/F"     );
  float          ele1_phi                ; regTree->Branch("ele1_phi"       , &ele1_phi       , "ele1_phi/F"     );
  float          ele1_chIso              ; regTree->Branch("ele1_chIso"     , &ele1_chIso     , "ele1_chIso/F"   ); 
  float          ele1_nhIso              ; regTree->Branch("ele1_nhIso"     , &ele1_nhIso     , "ele1_nhIso/F"   );
  float          ele1_phIso              ; regTree->Branch("ele1_phIso"     , &ele1_phIso     , "ele1_phIso/F"   );
  float          ele1_dxy                ; regTree->Branch("ele1_dxy"       , &ele1_dxy       , "ele1_dxy/F"     );
  float          ele1_dz                 ; regTree->Branch("ele1_dz"        , &ele1_dz        , "ele1_dz/F"      );
  float          ele1_sieie              ; regTree->Branch("ele1_sieie"     , &ele1_sieie     , "ele1_sieie"     ); 
  float          ele1_sipip              ; regTree->Branch("ele1_sipip"     , &ele1_sipip     , "ele1_sipip"     ); 
  float          ele1_r9                 ; regTree->Branch("ele1_r9"        , &ele1_r9        , "ele1_r9"        );
  float          ele1_dEtaInSeed         ; regTree->Branch("ele1_dEtaInSeed", &ele1_dEtaInSeed, "ele1_dEtaInSeed");
  float          ele1_dPhiIn             ; regTree->Branch("ele1_dPhiIn"    , &ele1_dPhiIn    , "ele1_dPhiIn"    );
  float          ele1_eseed              ; regTree->Branch("ele1_eseed"     , &ele1_eseed     , "ele1_eseed"     );
  float          ele1_hOverE             ; regTree->Branch("ele1_hOverE"    , &ele1_hOverE    , "ele1_hOverE"    );
  float          ele1_ecalE              ; regTree->Branch("ele1_ecalE"     , &ele1_ecalE     , "ele1_ecalE"     );
  float          ele1_trackP             ; regTree->Branch("ele1_trackP"    , &ele1_trackP    , "ele1_trackP"    );
  float          ele2_pt                 ; regTree->Branch("ele2_pt"        , &ele2_pt        , "ele2_pt/F"      );
  float          ele2_eta                ; regTree->Branch("ele2_eta"       , &ele2_eta       , "ele2_eta/F"     );
  float          ele2_phi                ; regTree->Branch("ele2_phi"       , &ele2_phi       , "ele2_phi/F"     );
  float          ele2_chIso              ; regTree->Branch("ele2_chIso"     , &ele2_chIso     , "ele2_chIso/F"   ); 
  float          ele2_nhIso              ; regTree->Branch("ele2_nhIso"     , &ele2_nhIso     , "ele2_nhIso/F"   );
  float          ele2_phIso              ; regTree->Branch("ele2_phIso"     , &ele2_phIso     , "ele2_phIso/F"   );
  float          ele2_dxy                ; regTree->Branch("ele2_dxy"       , &ele2_dxy       , "ele2_dxy/F"     );
  float          ele2_dz                 ; regTree->Branch("ele2_dz"        , &ele2_dz        , "ele2_dz/F"      );
  float          ele2_sieie              ; regTree->Branch("ele2_sieie"     , &ele2_sieie     , "ele2_sieie"     ); 
  float          ele2_sipip              ; regTree->Branch("ele2_sipip"     , &ele2_sipip     , "ele2_sipip"     ); 
  float          ele2_r9                 ; regTree->Branch("ele2_r9"        , &ele2_r9        , "ele2_r9"        );
  float          ele2_dEtaInSeed         ; regTree->Branch("ele2_dEtaInSeed", &ele2_dEtaInSeed, "ele2_dEtaInSeed");
  float          ele2_dPhiIn             ; regTree->Branch("ele2_dPhiIn"    , &ele2_dPhiIn    , "ele2_dPhiIn"    );
  float          ele2_eseed              ; regTree->Branch("ele2_eseed"     , &ele2_eseed     , "ele2_eseed"     );
  float          ele2_hOverE             ; regTree->Branch("ele2_hOverE"    , &ele2_hOverE    , "ele2_hOverE"    );
  float          ele2_ecalE              ; regTree->Branch("ele2_ecalE"     , &ele2_ecalE     , "ele2_ecalE"     );
  float          ele2_trackP             ; regTree->Branch("ele2_trackP"    , &ele2_trackP    , "ele2_trackP"    );
  
  float          mu1_pt                  ; regTree->Branch("mu1_pt"                    , &mu1_pt                  , "mu1_pt/F"                  );
  float          mu1_eta                 ; regTree->Branch("mu1_eta"                   , &mu1_eta                 , "mu1_eta/F"                 );
  float          mu1_phi                 ; regTree->Branch("mu1_phi"                   , &mu1_phi                 , "mu1_phi/F"                 );
  float          mu1_chIso               ; regTree->Branch("mu1_chIso"                 , &mu1_chIso               , "mu1_chIso/F"               ); 
  float          mu1_nhIso               ; regTree->Branch("mu1_nhIso"                 , &mu1_nhIso               , "mu1_nhIso/F"               );
  float          mu1_phIso               ; regTree->Branch("mu1_phIso"                 , &mu1_phIso               , "mu1_phIso/F"               );
  float          mu1_dxy                 ; regTree->Branch("mu1_dxy"                   , &mu1_dxy                 , "mu1_dxy/F"                 );
  float          mu1_dz                  ; regTree->Branch("mu1_dz"                    , &mu1_dz                  , "mu1_dz/F"                  );
  float          mu1_validFraction       ; regTree->Branch("mu1_validFraction"         , &mu1_validFraction       , "mu1_validFraction/F"       ); 
  unsigned short mu1_nValidMuon          ; regTree->Branch("mu1_nValidMuon"            , &mu1_nValidMuon          , "mu1_nValidMuon/s"          ); 
  unsigned short mu1_nValidPixel         ; regTree->Branch("mu1_nValidPixel"           , &mu1_nValidPixel         , "mu1_nValidPixel/s"         ); 
  unsigned short mu1_trkLayersWithMmt    ; regTree->Branch("mu1_trkLayersWithMmt"      , &mu1_trkLayersWithMmt    , "mu1_trkLayersWithMmt/s"    ); 
  unsigned short mu1_pixLayersWithMmt    ; regTree->Branch("mu1_pixLayersWithMmt"      , &mu1_pixLayersWithMmt    , "mu1_pixLayersWithMmt/s"    ); 
  unsigned short mu1_nMatched            ; regTree->Branch("mu1_nMatched"              , &mu1_nMatched            , "mu1_nMatched/s"            ); 
  float          mu1_normChi2            ; regTree->Branch("mu1_normChi2"              , &mu1_normChi2            , "mu1_normChi2/F"            ); 
  unsigned short mu1_chi2LocalPosition   ; regTree->Branch("mu1_chi2LocalPosition"     , &mu1_chi2LocalPosition   , "mu1_chi2LocalPosition/s"   ); 
  unsigned short mu1_trkKink             ; regTree->Branch("mu1_trkKink"               , &mu1_trkKink             , "mu1_trkKink/s"             ); 
  float          mu1_segmentCompatibility; regTree->Branch("mu1_segmentCompatibility"  , &mu1_segmentCompatibility, "mu1_segmentCompatibility/F"); 
  float          mu2_pt                  ; regTree->Branch("mu2_pt"                    , &mu2_pt                  , "mu2_pt/F"                  );
  float          mu2_eta                 ; regTree->Branch("mu2_eta"                   , &mu2_eta                 , "mu2_eta/F"                 );
  float          mu2_phi                 ; regTree->Branch("mu2_phi"                   , &mu2_phi                 , "mu2_phi/F"                 );
  float          mu2_chIso               ; regTree->Branch("mu2_chIso"                 , &mu2_chIso               , "mu2_chIso/F"               ); 
  float          mu2_nhIso               ; regTree->Branch("mu2_nhIso"                 , &mu2_nhIso               , "mu2_nhIso/F"               );
  float          mu2_phIso               ; regTree->Branch("mu2_phIso"                 , &mu2_phIso               , "mu2_phIso/F"               );
  float          mu2_dxy                 ; regTree->Branch("mu2_dxy"                   , &mu2_dxy                 , "mu2_dxy/F"                 );
  float          mu2_dz                  ; regTree->Branch("mu2_dz"                    , &mu2_dz                  , "mu2_dz/F"                  );
  float          mu2_validFraction       ; regTree->Branch("mu2_validFraction"         , &mu2_validFraction       , "mu2_validFraction/F"       ); 
  unsigned short mu2_nValidMuon          ; regTree->Branch("mu2_nValidMuon"            , &mu2_nValidMuon          , "mu2_nValidMuon/s"          ); 
  unsigned short mu2_nValidPixel         ; regTree->Branch("mu2_nValidPixel"           , &mu2_nValidPixel         , "mu2_nValidPixel/s"         ); 
  unsigned short mu2_trkLayersWithMmt    ; regTree->Branch("mu2_trkLayersWithMmt"      , &mu2_trkLayersWithMmt    , "mu2_trkLayersWithMmt/s"    ); 
  unsigned short mu2_pixLayersWithMmt    ; regTree->Branch("mu2_pixLayersWithMmt"      , &mu2_pixLayersWithMmt    , "mu2_pixLayersWithMmt/s"    ); 
  unsigned short mu2_nMatched            ; regTree->Branch("mu2_nMatched"              , &mu2_nMatched            , "mu2_nMatched/s"            ); 
  float          mu2_normChi2            ; regTree->Branch("mu2_normChi2"              , &mu2_normChi2            , "mu2_normChi2/F"            ); 
  unsigned short mu2_chi2LocalPosition   ; regTree->Branch("mu2_chi2LocalPosition"     , &mu2_chi2LocalPosition   , "mu2_chi2LocalPosition/s"   ); 
  unsigned short mu2_trkKink             ; regTree->Branch("mu2_trkKink"               , &mu2_trkKink             , "mu2_trkKink/s"             ); 
  float          mu2_segmentCompatibility; regTree->Branch("mu2_segmentCompatibility"  , &mu2_segmentCompatibility, "mu2_segmentCompatibility/F"); 
  float          sv1_pt                  ; regTree->Branch("sv1_pt"                    , &sv1_pt                  , "sv1_pt/F"                  ); 
  float          sv1_m                   ; regTree->Branch("sv1_m"                     , &sv1_m                   , "sv1_m/F"                   ); 
  float          sv1_3Dval               ; regTree->Branch("sv1_3Dval"                 , &sv1_3Dval               , "sv1_3Dval/F"               ); 
  float          sv1_3Derr               ; regTree->Branch("sv1_3Derr"                 , &sv1_3Derr               , "sv1_3Derr/F"               ); 
  unsigned short sv1_ntrk                ; regTree->Branch("sv1_ntrk"                  , &sv1_ntrk                , "sv1_ntrk/s"                ); 
  float          sv2_pt                  ; regTree->Branch("sv2_pt"                    , &sv2_pt                  , "sv2_pt/F"                  ); 
  float          sv2_m                   ; regTree->Branch("sv2_m"                     , &sv2_m                   , "sv2_m/F"                   ); 
  float          sv2_3Dval               ; regTree->Branch("sv2_3Dval"                 , &sv2_3Dval               , "sv2_3Dval/F"               ); 
  float          sv2_3Derr               ; regTree->Branch("sv2_3Derr"                 , &sv2_3Derr               , "sv2_3Derr/F"               ); 
  unsigned short sv2_ntrk                ; regTree->Branch("sv2_ntrk"                  , &sv2_ntrk                , "sv2_ntrk/s"                ); 
  unsigned char nEle, nMu;
  TLorentzVector sumNuV4, fatjetV4, nuV4, fatjetNuV4;
  Long64_t nEntries = tree->GetEntries();
  if(maxEntries>0) nEntries=TMath::Min(maxEntries,nEntries);
  Long64_t oneTenth = nEntries/10;
  for (Long64_t iEntry = 0; iEntry<nEntries; ++iEntry) {
    if(iEntry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",iEntry,nEntries); 
    event.getEntry(*tree, iEntry);
    
    for (unsigned nJ = 0; nJ<event.puppiAK8Jets.size(); nJ++){
      auto& fatjet = event.puppiAK8Jets[nJ];
      // Fatfatjet preselection
      if(fatjet.pt()<200) continue;
      // End preselection (That was short!)

      // Reset output tree branches
      nEle=0; nMu=0;
      sumNuV4.SetPtEtaPhiM(0,0,0,0);
      fatjet_pt               =  -1; 
      fatjet_eta              = -99; 
      fatjet_phi              =  -1; 
      fatjet_MSD              =  -1; 
      fatjet_Mnu              =  -1; 
      fatjet_nhf              =  -1;
      fatjet_chf              =  -1;
      fatjet_nef              =  -1;
      fatjet_cef              =  -1;
    
      sumNu_pt                =  -1;
      sumNu_eta               = -99;
      sumNu_phi               =  -1;
      sumNu_E                 =  -1;
      
      ele1_pt                 =  -1; 
      ele1_eta                = -99; 
      ele1_phi                =  -1; 
      ele1_chIso              =  -1; 
      ele1_nhIso              =  -1; 
      ele1_phIso              =  -1; 
      ele1_dxy                =  -7; 
      ele1_dz                 =  -7; 
      ele1_sieie              =  -1; 
      ele1_sipip              =  -1; 
      ele1_r9                 =  -1; 
      ele1_dEtaInSeed         =  -1; 
      ele1_dPhiIn             =  -1; 
      ele1_eseed              =  -1; 
      ele1_hOverE             =  -1; 
      ele1_ecalE              =  -1; 
      ele1_trackP             =  -1; 
      ele2_pt                 =  -1; 
      ele2_eta                = -99; 
      ele2_phi                =  -1; 
      ele2_chIso              =  -1; 
      ele2_nhIso              =  -1; 
      ele2_phIso              =  -1; 
      ele2_dxy                =  -7; 
      ele2_dz                 =  -7; 
      ele2_sieie              =  -1; 
      ele2_sipip              =  -1; 
      ele2_r9                 =  -1; 
      ele2_dEtaInSeed         =  -1; 
      ele2_dPhiIn             =  -1; 
      ele2_eseed              =  -1; 
      ele2_hOverE             =  -1; 
      ele2_ecalE              =  -1; 
      ele2_trackP             =  -1; 
      
      mu1_pt                  =  -1;
      mu1_eta                 = -99;
      mu1_phi                 =  -1;
      mu1_chIso               =  -1;
      mu1_nhIso               =  -1;
      mu1_phIso               =  -1;
      mu1_dxy                 =  -7;
      mu1_dz                  =  -7;
      mu1_validFraction       =  -1;
      mu1_nValidMuon          =  -1;
      mu1_nValidPixel         =  -1;
      mu1_trkLayersWithMmt    =  -1;
      mu1_pixLayersWithMmt    =  -1;
      mu1_nMatched            =  -1;
      mu1_normChi2            =  -1;
      mu1_chi2LocalPosition   =  -1;
      mu1_trkKink             =  -1;
      mu1_segmentCompatibility=  -1;
      mu2_pt                  =  -1;
      mu2_eta                 = -99;
      mu2_phi                 =  -1;
      mu2_chIso               =  -1;
      mu2_nhIso               =  -1;
      mu2_phIso               =  -1;
      mu2_dxy                 =  -7;
      mu2_dz                  =  -7;
      mu2_validFraction       =  -1;
      mu2_nValidMuon          =  -1;
      mu2_nValidPixel         =  -1;
      mu2_trkLayersWithMmt    =  -1;
      mu2_pixLayersWithMmt    =  -1;
      mu2_nMatched            =  -1;
      mu2_normChi2            =  -1;
      mu2_chi2LocalPosition   =  -1;
      mu2_trkKink             =  -1;
      mu2_segmentCompatibility=  -1;

      sv1_pt    = -1;
      sv1_m     = -1;
      sv1_3Dval = -1;
      sv1_3Derr = -1;
      sv1_ntrk  = -1;
      sv2_pt    = -1;
      sv2_m     = -1;
      sv2_3Dval = -1;
      sv2_3Derr = -1;
      sv2_ntrk  = -1;
      
      // Populate output tree info
      fatjet_pt               = fatjet.pt(); 
      fatjet_eta              = fatjet.eta();
      fatjet_phi              = fatjet.phi();
      fatjet_MSD              = fatjet.mSD;
      fatjet_nhf              = fatjet.nhf;
      fatjet_chf              = fatjet.chf; 
      fatjet_nef              = fatjet.nef; 
      fatjet_cef              = fatjet.cef; 

      // Final state neutrinos inside the fatjet
      // No parentage calculation for now
      std::vector<const panda::GenParticle*> validGenP;
      for (unsigned iG = 0; iG<event.genParticles.size(); iG++) {
        auto& genParticle = event.genParticles[iG];
        if (genParticle.finalState != 1) continue;
        unsigned absid = abs(genParticle.pdgid);
        if (absid!=12 && absid!=14 && absid!=16) continue;
        if (fatjet.dR2(genParticle)>0.64) continue;
        
        bool isDuplicate=false;
        for (auto* vgp : validGenP)
          if(genParticle.pdgid==vgp->pdgid &&
            fabs(genParticle.pt()/vgp->pt()-1)<0.10 &&
            sqrt(pow(genParticle.eta()-vgp->eta(),2)+pow(TVector2::Phi_mpi_pi(genParticle.phi()-vgp->phi()),2)) < 0.00001
          ) { isDuplicate=true; break; }
        if(isDuplicate) continue;
        else validGenP.push_back(&genParticle);
        nuV4.SetPtEtaPhiM(genParticle.pt(), genParticle.eta(), genParticle.phi(), 0);
        sumNuV4+=nuV4;
      }
      fatjetV4.SetPtEtaPhiM(fatjet_pt, fatjet_eta, fatjet_phi, fatjet_MSD);
      fatjetNuV4 = sumNuV4 + fatjetV4;
      fatjet_Mnu =  sumNuV4.M(); 
      
      sumNu_pt                = sumNuV4.Pt() ;
      sumNu_eta               = sumNuV4.Eta();
      sumNu_phi               = sumNuV4.Phi();
      sumNu_E                 = sumNuV4.E()  ;

      // Properties of the leading electrons
      for(unsigned iE = 0; iE<event.electrons.size(); iE++) {
        auto& electron  = event.electrons[iE];
        if(nEle>=2 || electron.pt()<10) break;
        if(fabs(electron.eta()) > 2.5) continue;
        if(fatjet.dR2(electron)>0.64) continue;
        if(nEle==0) {
          ele1_pt        =electron.pt()      ;
          ele1_eta       =electron.eta()     ;
          ele1_phi       =electron.phi()     ;
          ele1_chIso     =electron.chIso     ;
          ele1_nhIso     =electron.nhIso     ;
          ele1_phIso     =electron.phIso     ;
          ele1_dxy       =electron.dxy       ;
          ele1_dz        =electron.dz        ;
          ele1_sieie     =electron.sieie     ;
          ele1_sipip     =electron.sipip     ;
          ele1_r9        =electron.r9        ;
          ele1_dEtaInSeed=electron.dEtaInSeed;
          ele1_dPhiIn    =electron.dPhiIn    ;
          ele1_eseed     =electron.eseed     ;
          ele1_hOverE    =electron.hOverE    ;
          ele1_ecalE     =electron.ecalE     ;
          ele1_trackP    =electron.trackP    ;
        } else if(nEle==1) {
          ele2_pt        =electron.pt()      ;
          ele2_eta       =electron.eta()     ;
          ele2_phi       =electron.phi()     ;
          ele2_chIso     =electron.chIso     ;
          ele2_nhIso     =electron.nhIso     ;
          ele2_phIso     =electron.phIso     ;
          ele2_dxy       =electron.dxy       ;
          ele2_dz        =electron.dz        ;
          ele2_sieie     =electron.sieie     ;
          ele2_sipip     =electron.sipip     ;
          ele2_r9        =electron.r9        ;
          ele2_dEtaInSeed=electron.dEtaInSeed;
          ele2_dPhiIn    =electron.dPhiIn    ;
          ele2_eseed     =electron.eseed     ;
          ele2_hOverE    =electron.hOverE    ;
          ele2_ecalE     =electron.ecalE     ;
          ele2_trackP    =electron.trackP    ;
        }
        nEle++;
      }
      // Properties of the leading muons 
      for(unsigned iM = 0; iM<event.muons.size(); iM++){
        auto& muon  = event.muons[iM];
        if(nMu>=2 || muon.pt()<5) break;
        if(!muon.global || fabs(muon.eta())>2.4) continue;
        if(fatjet.dR2(muon)>0.64) continue;
        if(nMu==0) {
          mu1_pt                  =muon.pt()                ;
          mu1_eta                 =muon.eta()               ;
          mu1_phi                 =muon.phi()               ;
          mu1_chIso               =muon.chIso               ;
          mu1_nhIso               =muon.nhIso               ;
          mu1_phIso               =muon.phIso               ;
          mu1_dxy                 =muon.dxy                 ;
          mu1_dz                  =muon.dz                  ;
          mu1_validFraction       =muon.validFraction       ;
          mu1_nValidMuon          =muon.nValidMuon          ;
          mu1_nValidPixel         =muon.nValidPixel         ;
          mu1_trkLayersWithMmt    =muon.trkLayersWithMmt    ;
          mu1_pixLayersWithMmt    =muon.pixLayersWithMmt    ;
          mu1_nMatched            =muon.nMatched            ;
          mu1_normChi2            =muon.normChi2            ;
          mu1_chi2LocalPosition   =muon.chi2LocalPosition   ;
          mu1_trkKink             =muon.trkKink             ;
          mu1_segmentCompatibility=muon.segmentCompatibility;
        } else if(nMu==1) {
          mu2_pt                  =muon.pt()                ;
          mu2_eta                 =muon.eta()               ;
          mu2_phi                 =muon.phi()               ;
          mu2_chIso               =muon.chIso               ;
          mu2_nhIso               =muon.nhIso               ;
          mu2_phIso               =muon.phIso               ;
          mu2_dxy                 =muon.dxy                 ;
          mu2_dz                  =muon.dz                  ;
          mu2_validFraction       =muon.validFraction       ;
          mu2_nValidMuon          =muon.nValidMuon          ;
          mu2_nValidPixel         =muon.nValidPixel         ;
          mu2_trkLayersWithMmt    =muon.trkLayersWithMmt    ;
          mu2_pixLayersWithMmt    =muon.pixLayersWithMmt    ;
          mu2_nMatched            =muon.nMatched            ;
          mu2_normChi2            =muon.normChi2            ;
          mu2_chi2LocalPosition   =muon.chi2LocalPosition   ;
          mu2_trkKink             =muon.trkKink             ;
          mu2_segmentCompatibility=muon.segmentCompatibility;
        }
        nMu++;
      }

      // AK8 -> AK4 jets within -> Their secondary vertices
      vector<panda::SecondaryVertex*> sv_;
      for (unsigned nJ2 = 0; nJ<event.chsAK4Jets.size(); nJ2++){
        auto& jet = event.chsAK4Jets[nJ2];
        if(jet.pt()<20) break;
        if(fabs(jet.eta())>2.4) continue;
        if(jet.dR2(fatjet)>0.64) continue;
        if(!jet.secondaryVertex.isValid()) continue;
        panda::SecondaryVertex* sv = (panda::SecondaryVertex*)jet.secondaryVertex.get();
        sv_.push_back(sv);
      }
      sort(
        sv_.begin(),
        sv_.end(),
          [](panda::SecondaryVertex *x, panda::SecondaryVertex *y) -> bool { return x->pt()  > y->pt() ; }
      );
      if(sv_.size()>0) {
        sv1_pt    = sv_[0]->pt();
        sv1_m     = sv_[0]->m();
        sv1_3Dval = sv_[0]->vtx3DVal;
        sv1_3Derr = sv_[0]->vtx3DeVal;
        sv1_ntrk  = sv_[0]->ntrk;
      } if(sv_.size()>1) {
        sv2_pt    = sv_[1]->pt();
        sv2_m     = sv_[1]->m();
        sv2_3Dval = sv_[1]->vtx3DVal;
        sv2_3Derr = sv_[1]->vtx3DeVal;
        sv2_ntrk  = sv_[1]->ntrk;
      }
      // Fill output tree
      regTree->Fill();
    }

  }
  regTree->Write("regTree",TObject::kOverwrite);
  outputFile->Close();
  inputFile->Close();
 
}

