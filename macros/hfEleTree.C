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
//#include "RoccoR.cc"
//#include <TRandom3.h>

void hfEleTree(
  TString inputFileName,
  TString outputFileName,
  Long64_t maxEntries=-1,
  int debug=0
) {
  gSystem->Load("libPandaTreeObjects.so");
  TFile *inputFile=0;
  int retries=0;
  while(true) {
    inputFile = TFile::Open(inputFileName,"read");
    if(inputFile && inputFile->IsOpen()) break;
    retries++;
    if(retries>100) { throw std::runtime_error("Error opening input file"); return; }
  }
    
  TTree* events = (TTree*)inputFile->Get("events"); // get the tree object from the file
  panda::Event event; // create an Event object
  event.setStatus(*events, {"!*"});
  event.setAddress(*events, {
    "runNumber", "lumiNumber", "eventNumber",
    "electrons", "taus", "photons",
    "genParticles", "pfCandidates","tracks",
    "puppiAK8Jets", "chsAK4Jets",
    "weight"
    //"pfCandidates","vertices","weight","tracks",
  }); 

  TH1D* hSumW=(TH1D*)inputFile->FindObjectAny("hSumW");
  Long64_t sumMcWeights = hSumW->GetBinContent(1);
  
  TFile *outputFile = new TFile(outputFileName, "recreate");
  TTree *tree = new TTree("tree","tree");
  float weight             ; tree->Branch("weight"            , &weight             , "weight/F"            );
  float fatjet_pt          ; tree->Branch("fatjet_pt"         , &fatjet_pt          , "fatjet_pt/F"         );
  float dEta_eleFatjet     ; tree->Branch("dEta_eleFatjet"    , &dEta_eleFatjet     , "dEta_eleFatjet/F"    );
  float dPhi_eleFatjet     ; tree->Branch("dPhi_eleFatjet"    , &dPhi_eleFatjet     , "dPhi_eleFatjet/F"    );
  float ele_pt             ; tree->Branch("ele_pt"            , &ele_pt             , "ele_pt/F"            );
  float ele_eta            ; tree->Branch("ele_eta"           , &ele_eta            , "ele_eta/F"           );
  float ele_phi            ; tree->Branch("ele_phi"           , &ele_phi            , "ele_phi/F"           );
  float ele_chIso          ; tree->Branch("ele_chIso"         , &ele_chIso          , "ele_chIso/F"         ); 
  float ele_nhIso          ; tree->Branch("ele_nhIso"         , &ele_nhIso          , "ele_nhIso/F"         );
  float ele_phIso          ; tree->Branch("ele_phIso"         , &ele_phIso          , "ele_phIso/F"         );
  float ele_ecalIso        ; tree->Branch("ele_ecalIso"       , &ele_ecalIso        , "ele_ecalIso/F"       );
  float ele_hcalIso        ; tree->Branch("ele_hcalIso"       , &ele_hcalIso        , "ele_hcalIso/F"       );
  float ele_trackIso       ; tree->Branch("ele_trackIso"      , &ele_trackIso       , "ele_trackIso/F"      );
  float ele_combIso        ; tree->Branch("ele_combIso"       , &ele_combIso        , "ele_combIso/F"       );
  float ele_dxy            ; tree->Branch("ele_dxy"           , &ele_dxy            , "ele_dxy/F"           );
  float ele_dz             ; tree->Branch("ele_dz"            , &ele_dz             , "ele_dz/F"            );
  float ele_sieie          ; tree->Branch("ele_sieie"         , &ele_sieie          , "ele_sieie/F"         ); 
  float ele_sipip          ; tree->Branch("ele_sipip"         , &ele_sipip          , "ele_sipip/F"         ); 
  float ele_r9             ; tree->Branch("ele_r9"            , &ele_r9             , "ele_r9/F"            );
  float ele_dEtaInSeed     ; tree->Branch("ele_dEtaInSeed"    , &ele_dEtaInSeed     , "ele_dEtaInSeed/F"    );
  float ele_dPhiIn         ; tree->Branch("ele_dPhiIn"        , &ele_dPhiIn         , "ele_dPhiIn/F"        );
  float ele_hOverE         ; tree->Branch("ele_hOverE"        , &ele_hOverE         , "ele_hOverE/F"        );
  float ele_ooEmooP        ; tree->Branch("ele_ooEmooP"       , &ele_ooEmooP        , "ele_ooEmooP/F"       );
  int   ele_nMissingHits   ; tree->Branch("ele_nMissingHits"  , &ele_nMissingHits   , "ele_nMissingHits/I"  );
  int   ele_conversionVeto ; tree->Branch("ele_conversionVeto", &ele_conversionVeto , "ele_conversionVeto/I");
  int   ele_tripleCharge   ; tree->Branch("ele_tripleCharge"  , &ele_tripleCharge   , "ele_tripleCharge/I"  );
  int   ele_loose          ; tree->Branch("ele_loose"         , &ele_loose          , "ele_loose/I"         );
  int   ele_tight          ; tree->Branch("ele_tight"         , &ele_tight          , "ele_tight/I"         );
  int   ele_matchedGHFE    ; tree->Branch("ele_matchedGHFE"   , &ele_matchedGHFE    , "ele_matchedGHFE/I"   );
  float ele_ptBal_GHFE     ; tree->Branch("ele_ptBal_GHFE"    , &ele_ptBal_GHFE     , "ele_ptBal_GHFE/F"    );
  float ele_deltaR_GHFE    ; tree->Branch("ele_deltaR_GHFE"   , &ele_deltaR_GHFE    , "ele_deltaR_GHFE/F"   );
  TTree *genTree = new TTree("genTree","genTree");
  float genEle_pt   ; genTree->Branch("genEle_pt"   , &genEle_pt    , "genEle_pt/F"   );
  float genEle_eta  ; genTree->Branch("genEle_eta"  , &genEle_eta   , "genEle_eta/F"  );
  float genEle_phi  ; genTree->Branch("genEle_phi"  , &genEle_phi   , "genEle_phi/F"  );
  float genNu_pt    ; genTree->Branch("genNu_pt"    , &genNu_pt     , "genNu_pt/F"    );
  float genNu_eta   ; genTree->Branch("genNu_eta"   , &genNu_eta    , "genNu_eta/F"   );
  float genNu_phi   ; genTree->Branch("genNu_phi"   , &genNu_phi    , "genNu_phi/F"   );
  int   foundRecoEle; genTree->Branch("foundRecoEle", &foundRecoEle , "foundRecoEle/I");
  int   foundPFCand ; genTree->Branch("foundPFCand" , &foundPFCand  , "foundPFCand/I" );
  genTree->Branch("weight", &weight, "weight/F");
  
  TLorentzVector sumNuV4, fatjetV4, nuV4, fatjetNuV4;
  Long64_t nEntries = events->GetEntries();
  if(maxEntries>0) nEntries=TMath::Min(maxEntries,nEntries);
  Long64_t oneTenth = nEntries/10;
  for (Long64_t iEntry = 0; iEntry<nEntries; ++iEntry) {
    if(!debug&&iEntry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",iEntry,nEntries); 
    // Reset output tree branches
    event.getEntry(*events, iEntry);
    
    // Look for a 250 GeV reco fatjet passing the loose ID
    panda::FatJet *theFatJet=0;
    for (unsigned nJ = 0; nJ<event.puppiAK8Jets.size(); nJ++){
      auto& fatjet = event.puppiAK8Jets[nJ];
      // Fatjet preselection
      if(fatjet.pt()<250) break;
      if(fabs(fatjet.eta())>2.4) continue;
      if(!fatjet.loose) continue;
      theFatJet = &fatjet;
    }
    if(!theFatJet) continue;
    if(debug) printf("######## Reading entry %lld/%lld ########################################################\n",iEntry,nEntries); 
    if(debug) printf("  Fatjet (pt,eta,phi)=(%.1f,%.2f,%.2f)\n",theFatJet->pt(),theFatJet->eta(),theFatJet->phi());
    fatjet_pt = theFatJet->pt();
    std::vector<const panda::GenParticle*> validGenP, genHFElectrons, genNuSisters;
    std::vector<panda::Electron*> genElectronRecoMatch;
    if(debug>=3) printf("  Looping over %u gen particles\n",event.genParticles.size());
    for (unsigned iG = 0; iG<event.genParticles.size(); iG++) {
      auto& genParticle = event.genParticles[iG];
      unsigned absid = abs(genParticle.pdgid);
      if (absid!=11) continue;
      if (genParticle.finalState != 1) continue;
      if(genParticle.pt()<10.) continue;
      if(fabs(genParticle.eta())>2.5) continue;
      if(theFatJet->dR2(genParticle)>0.64) continue;
      short parentId=0; if(genParticle.parent.isValid()) parentId = int(genParticle.parent.get()->pdgid);
      if(debug>=3) printf("    Gen electron at index %u (pt,eta,phi)=(%.2f,%.2f,%.2f), id %d, finalState %d, parent %d \n", iG, genParticle.pt(), genParticle.eta(), genParticle.phi(), genParticle.pdgid, genParticle.finalState, parentId);
      
      bool isDuplicate=false;
      for (auto* vgp : validGenP)
        if(genParticle.pdgid==vgp->pdgid &&
          fabs(genParticle.pt()/vgp->pt()-1)<0.10 &&
          sqrt(pow(genParticle.eta()-vgp->eta(),2)+pow(TVector2::Phi_mpi_pi(genParticle.phi()-vgp->phi()),2)) < 0.00001
        ) { isDuplicate=true; break; }
      if(isDuplicate) continue;
      validGenP.push_back(&genParticle);
      
      // Find a HF hadron
      // Allowed codes: 400-499 (charmed meson), 500-599 (B meson),
      // above 4000-5999 (Lambda_C^+ or something), and some weird D/B mesons over 10000
      // To do: go up the chain for taus
      unsigned short absParentIdMod10k = unsigned(abs(parentId))%10000;
      bool fromHFhadron = (absParentIdMod10k>=400 && absParentIdMod10k<600) || (absParentIdMod10k>=4000 && absParentIdMod10k<6000);
      if(!fromHFhadron) continue;
      
      genHFElectrons.push_back(&genParticle);
      // Find the neutrino sister
      bool foundSister=false; int target=(genParticle.pdgid==11? -12 : 12);
      for (unsigned iG2 = 0; iG2<event.genParticles.size(); iG2++) {
        auto& genParticle2 = event.genParticles[iG2];
        if(!genParticle2.parent.isValid()) continue;
        if(genParticle2.parent.get() != genParticle.parent.get()) continue;
        if(genParticle2.pdgid != target) continue;
        foundSister=true;
        genNuSisters.push_back(&genParticle2);
      }
      if(!foundSister) genNuSisters.push_back(nullptr);


    }
    genElectronRecoMatch.resize(genHFElectrons.size());
    for(unsigned iGHFE=0; iGHFE<genHFElectrons.size(); iGHFE++)
      genElectronRecoMatch[iGHFE]=nullptr;
    if(debug) printf("  Found %zu gen electrons from HF hadrons\n", genHFElectrons.size());
    
    weight = event.weight/sumMcWeights;
    
    // Properties of the leading electrons
    for(unsigned iE = 0; iE<event.electrons.size(); iE++) {
      auto& electron  = event.electrons[iE];
      float trackPt = electron.trackP / TMath::CosH(electron.eta());
      //if(trackPt<5.) continue;
      if(fabs(electron.eta()) > 2.5) continue;
      if(theFatJet->dR2(electron)>0.64) continue;
      bool matchedToHFHadron=false; float ptBalance_GHFE=3, dR_GHFE=3; // dummy distances
      if(genHFElectrons.size()>0 && debug>=3) 
        printf("  Attempting to match reco electron of (pt,eta,phi)=(%.2f,%.2f,%.2f) to %zu gen HF electrons\n",trackPt,electron.eta(),electron.phi(),genHFElectrons.size());
      for(unsigned iGHFE=0; iGHFE<genHFElectrons.size(); iGHFE++) {
        panda::GenParticle *genHFEle = (panda::GenParticle*)genHFElectrons[iGHFE];
        ptBalance_GHFE = fabs(trackPt/genHFEle->pt() - 1.);
        dR_GHFE = genHFEle->dR(electron);
        //if(ptBalance_GHFE < 0.5 && dR_GHFE < 0.01)  matchedToHFHadron=true;
        if(dR_GHFE < 0.1)  matchedToHFHadron=true;
        if(debug>=3) printf("    Gen HF electron %d has pt balance %.2f, dR %.2f => matched=%d\n",iGHFE,ptBalance_GHFE,dR_GHFE,matchedToHFHadron);
        if(matchedToHFHadron) {
          genElectronRecoMatch[iGHFE] = &electron;
          break;
        }
      }
      dEta_eleFatjet     = fabs(electron.eta()-theFatJet->eta());
      dPhi_eleFatjet     = fabs(TVector2::Phi_mpi_pi(electron.phi()-theFatJet->phi()));
      ele_pt             = trackPt              ;
      ele_eta            = electron.eta()       ;
      ele_phi            = electron.phi()       ;
      ele_chIso          = electron.chIso       ;
      ele_nhIso          = electron.nhIso       ;
      ele_phIso          = electron.phIso       ;
      ele_ecalIso        = electron.ecalIso     ;
      ele_hcalIso        = electron.hcalIso     ;
      ele_trackIso       = electron.trackIso    ;
      ele_combIso        = electron.combIso()   ;
      ele_dxy            = electron.dxy         ;
      ele_dz             = electron.dz          ;
      ele_sieie          = electron.sieie       ;
      ele_sipip          = electron.sipip       ;
      ele_r9             = electron.r9          ;
      ele_dEtaInSeed     = electron.dEtaInSeed  ;
      ele_dPhiIn         = electron.dPhiIn      ;
      ele_hOverE         = electron.hOverE      ;
      ele_ooEmooP        = fabs(1.-electron.eseed/max(.0001f,electron.trackP))/electron.ecalE;
      ele_nMissingHits   = electron.nMissingHits;
      ele_conversionVeto = electron.conversionVeto;
      ele_tripleCharge   = electron.tripleCharge;
      ele_loose          = electron.loose;
      ele_tight          = electron.tight;
      ele_matchedGHFE    = matchedToHFHadron    ;
      if(matchedToHFHadron) { 
        ele_ptBal_GHFE   = ptBalance_GHFE;
        ele_deltaR_GHFE  = dR_GHFE;
      } else {
        ele_ptBal_GHFE   = 3;
        ele_deltaR_GHFE  = 3;
      }
      // Fill output tree
      if(debug) printf("  Filling tree with a reco electron (gen HF match=%d)\n",matchedToHFHadron);
      tree->Fill();
    }
    for(unsigned iGHFE=0; iGHFE<genHFElectrons.size(); iGHFE++) {
      foundRecoEle = 0;
      foundPFCand = 0;
      if(genElectronRecoMatch[iGHFE]!=nullptr) foundRecoEle = 1;
      if(debug>=3) printf("  Looking for charged PF candidate to match to gen HF electron (pt,eta,phi)=(%.2f,%.2f,%.2f)\n",
        genHFElectrons[iGHFE]->pt(),genHFElectrons[iGHFE]->eta(),genHFElectrons[iGHFE]->phi());
      for (auto& pfCand : event.pfCandidates) {
        bool isRightType=false;
        switch( pfCand.ptype) {
          case panda::PFCand::hp:
          case panda::PFCand::ep:
          case panda::PFCand::mup:
          case panda::PFCand::hm:
          case panda::PFCand::em:
          case panda::PFCand::mum:
            isRightType=true; break;
          default: break;
        }
        if(!isRightType) continue;
        if (!pfCand.track.isValid()) continue;
        if(pfCand.pt()<5.) continue; // not optimal, should use pfRangeMax of the vertices to skip unpacking the soft pf candidates from each vertex
        if(pfCand.eta()<-2.5 || pfCand.eta()>2.5) continue;
        float ptBalance_GHFE = fabs(pfCand.pt()/genHFElectrons[iGHFE]->pt() - 1.);
        float dR_GHFE = genHFElectrons[iGHFE]->dR(pfCand);
        if(debug>=3) printf("    charged PF candidate: ptype=%d, (pt,eta,phi)=(%.2f,%.2f,%.2f)\n", pfCand.ptype, pfCand.pt(), pfCand.eta(), pfCand.phi());
        if(ptBalance_GHFE<0.5 && dR_GHFE<0.01) {
          if(debug>=3) printf("      above candidate matches with pt balance %.2f, dR %.3f\n", ptBalance_GHFE,dR_GHFE);
          foundPFCand=1;
          break; 
        }
      }
      // Look for photon nearby?
      if(foundPFCand) {
        if(debug>=3) printf("  Looking for a photon too\n");
          for(auto& photon: event.photons) {


          }
      }
      genEle_pt   = genHFElectrons[iGHFE]->pt();
      genEle_eta  = genHFElectrons[iGHFE]->eta();
      genEle_phi  = genHFElectrons[iGHFE]->phi();
      genNu_pt    = -1;
      genNu_eta   = -7;
      genNu_phi   = -7;
      if(genNuSisters[iGHFE]!=nullptr) {
        genNu_pt    = genNuSisters[iGHFE]->pt();
        genNu_eta   = genNuSisters[iGHFE]->eta();
        genNu_phi   = genNuSisters[iGHFE]->phi();
      }
      genTree->Fill();
    }

  }
  tree->Write("tree",TObject::kOverwrite);
  genTree->Write("genTree",TObject::kOverwrite);
  outputFile->Close();
  inputFile->Close();
 
}

