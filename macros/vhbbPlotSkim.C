#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unistd.h>
enum selectionType { 
  kWHLightFlavorCR, kWHHeavyFlavorCR, kWH2TopCR, kWHSR, 
  nSelectionTypes
};
enum sampleType {
  kData       , // 0 
  kQCD        , // 1
  kVZ         , // 2
  kWW         , // 3
  kTT         , // 4
  kTop        , // 5
  kWjets      , // 6
  kZjets      , // 7
  kVH         , // 7
  nSampleTypes 
};
enum plotCategory {
  kPlotData , // 0  
  kPlotQCD  , // 1  
  kPlotVZbb , // 2  
  kPlotVVLF , // 3  
  kPlotTT   , // 4  
  kPlotTop  , // 5  
  kPlotWbb  , // 6  
  kPlotWb   , // 7  
  kPlotWLF  , // 8  
  kPlotZbb  , // 9  
  kPlotZb   , //10  
  kPlotZLF  , //11   
  kPlotVH   , //12   
  nPlotCategories
};
// This function loads the ith entry of the branch, only if it has not already been loaded
int bLoad(TBranch *branch, Long64_t ientry) {
  int bytesRead=0;
  Long64_t readEntry = branch->GetReadEntry();
  if(readEntry != ientry) bytesRead = branch->GetEntry(ientry);
  return bytesRead;
}

void vhbbPlotSkim(
  TString inputFileName="", 
  TString outputFileName="",
  sampleType sample=nSampleTypes,
  selectionType selection=nSelectionTypes,
  bool debug=false,
  Long64_t maxEntries=0
) {
  const double theLumi=35900.;
  const bool usePandaHbb = true;
  const float bDiscrLoose  = 0.5426; 
  const float bDiscrMedium = 0.8484; 
  const float bDiscrTight  = 0.9535; 
  
  // Handle File IO
  assert(inputFileName!=outputFileName);
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  assert(inputFile && inputFile->IsOpen());
  TTree *events = (TTree*)inputFile->Get("events");
  assert(events);

  Int_t           nJet;
  Int_t           nJot;
  Int_t           nLooseLep;
  Int_t           nLooseElectron;
  Int_t           nLooseMuon;
  Int_t           nTightLep;
  Int_t           nTightElectron;
  Int_t           nTightMuon;
  Float_t         muonPt[4];
  Float_t         muonEta[4];
  Float_t         muonPhi[4];
  Float_t         muonSfLoose[4];
  Float_t         muonSfTight[4];
  Float_t         muonSfUnc[4];
  Float_t         muonSfReco[4];
  Int_t           muonSelBit[4];
  Int_t           muonPdgId[4];
  Float_t         electronPt[4];
  Float_t         electronEta[4];
  Float_t         electronPhi[4];
  Float_t         electronSfLoose[4];
  Float_t         electronSfTight[4];
  Float_t         electronSfUnc[4];
  Float_t         electronSfReco[4];
  Int_t           electronSelBit[4];
  Float_t         jetPt[15];
  Float_t         jetEta[15];
  Float_t         jetPhi[15];
  Float_t         jetE[15];
  Float_t         jetCSV[15];
  Float_t         jetIso[15];
  Float_t         jetQGL[15];
  Float_t         jetLeadingLepPt[15];
  Float_t         jetLeadingTrkPt[15];
  Float_t         jetEMFrac[15];
  Float_t         jetHadFrac[15];
  Int_t           jetNLep[15];
  Float_t         jetGenPt[15];
  Int_t           jetGenFlavor[15];
  Float_t         fj1sjPt[2];
  Float_t         fj1sjPhi[2];
  Float_t         fj1sjEta[2];
  Float_t         fj1sjM[2];
  Float_t         fj1sjCSV[2];
  Float_t         fj1sjQGL[2];
  Int_t           fj1Nbs;
  Int_t           fj1gbb;
  Float_t         hbbpt;
  Float_t         hbbeta;
  Float_t         hbbphi;
  Float_t         hbbm;
  Int_t           hbbjtidx[2];
  Float_t         scale[6];
  Float_t         fj1ECFN_1_1_05;
  Float_t         fj1ECFN_2_1_05;
  Float_t         fj1ECFN_3_1_05;
  Float_t         fj1ECFN_1_2_05;
  Float_t         fj1ECFN_2_2_05;
  Float_t         fj1ECFN_3_2_05;
  Float_t         fj1ECFN_1_3_05;
  Float_t         fj1ECFN_2_3_05;
  Float_t         fj1ECFN_3_3_05;
  Float_t         fj1ECFN_1_4_05;
  Float_t         fj1ECFN_2_4_05;
  Float_t         fj1ECFN_3_4_05;
  Float_t         fj1ECFN_1_1_10;
  Float_t         fj1ECFN_2_1_10;
  Float_t         fj1ECFN_3_1_10;
  Float_t         fj1ECFN_1_2_10;
  Float_t         fj1ECFN_2_2_10;
  Float_t         fj1ECFN_3_2_10;
  Float_t         fj1ECFN_1_3_10;
  Float_t         fj1ECFN_2_3_10;
  Float_t         fj1ECFN_3_3_10;
  Float_t         fj1ECFN_1_4_10;
  Float_t         fj1ECFN_2_4_10;
  Float_t         fj1ECFN_3_4_10;
  Float_t         fj1ECFN_1_1_20;
  Float_t         fj1ECFN_2_1_20;
  Float_t         fj1ECFN_3_1_20;
  Float_t         fj1ECFN_1_2_20;
  Float_t         fj1ECFN_2_2_20;
  Float_t         fj1ECFN_3_2_20;
  Float_t         fj1ECFN_1_3_20;
  Float_t         fj1ECFN_2_3_20;
  Float_t         fj1ECFN_3_3_20;
  Float_t         fj1ECFN_1_4_20;
  Float_t         fj1ECFN_2_4_20;
  Float_t         fj1ECFN_3_4_20;
  Float_t         fj1ECFN_1_1_40;
  Float_t         fj1ECFN_2_1_40;
  Float_t         fj1ECFN_3_1_40;
  Float_t         fj1ECFN_1_2_40;
  Float_t         fj1ECFN_2_2_40;
  Float_t         fj1ECFN_3_2_40;
  Float_t         fj1ECFN_1_3_40;
  Float_t         fj1ECFN_2_3_40;
  Float_t         fj1ECFN_3_3_40;
  Float_t         fj1ECFN_1_4_40;
  Float_t         fj1ECFN_2_4_40;
  Float_t         fj1ECFN_3_4_40;
  Float_t         sf_btag0;
  Float_t         sf_btag1;
  Float_t         sf_btag2;
  Float_t         sf_btagGT0;
  Float_t         sf_sjbtag0;
  Float_t         sf_sjbtag1;
  Float_t         sf_sjbtag2;
  Float_t         sf_sjbtagGT0;
  Float_t         sf_btag0BUp;
  Float_t         sf_btag1BUp;
  Float_t         sf_btag2BUp;
  Float_t         sf_btagGT0BUp;
  Float_t         sf_sjbtag0BUp;
  Float_t         sf_sjbtag1BUp;
  Float_t         sf_sjbtag2BUp;
  Float_t         sf_sjbtagGT0BUp;
  Float_t         sf_btag0BDown;
  Float_t         sf_btag1BDown;
  Float_t         sf_btag2BDown;
  Float_t         sf_btagGT0BDown;
  Float_t         sf_sjbtag0BDown;
  Float_t         sf_sjbtag1BDown;
  Float_t         sf_sjbtag2BDown;
  Float_t         sf_sjbtagGT0BDown;
  Float_t         sf_btag0MUp;
  Float_t         sf_btag1MUp;
  Float_t         sf_btag2MUp;
  Float_t         sf_btagGT0MUp;
  Float_t         sf_sjbtag0MUp;
  Float_t         sf_sjbtag1MUp;
  Float_t         sf_sjbtag2MUp;
  Float_t         sf_sjbtagGT0MUp;
  Float_t         sf_btag0MDown;
  Float_t         sf_btag1MDown;
  Float_t         sf_btag2MDown;
  Float_t         sf_btagGT0MDown;
  Float_t         sf_sjbtag0MDown;
  Float_t         sf_sjbtag1MDown;
  Float_t         sf_sjbtag2MDown;
  Float_t         sf_sjbtagGT0MDown;
  Float_t         muonD0[4];
  Float_t         muonDZ[4];
  Float_t         muonSfMedium[4];
  Int_t           muonIsSoftMuon[4];
  Int_t           muonIsGlobalMuon[4];
  Int_t           muonIsTrackerMuon[4];
  Int_t           muonNValidMuon[4];
  Int_t           muonNValidPixel[4];
  Int_t           muonTrkLayersWithMmt[4];
  Int_t           muonPixLayersWithMmt[4];
  Int_t           muonNMatched[4];
  Int_t           muonChi2LocalPosition[4];
  Int_t           muonTrkKink[4];
  Float_t         muonValidFraction[4];
  Float_t         muonNormChi2[4];
  Float_t         muonSegmentCompatibility[4];
  Int_t           electronPdgId[4];
  Float_t         electronD0[4];
  Float_t         electronDZ[4];
  Float_t         electronSfMedium[4];
  Float_t         electronChIsoPh[4];
  Float_t         electronNhIsoPh[4];
  Float_t         electronPhIsoPh[4];
  Float_t         electronEcalIso[4];
  Float_t         electronHcalIso[4];
  Float_t         electronTrackIso[4];
  Float_t         electronIsoPUOffset[4];
  Float_t         electronSieie[4];
  Float_t         electronSipip[4];
  Float_t         electronDEtaInSeed[4];
  Float_t         electronDPhiIn[4];
  Float_t         electronEseed[4];
  Float_t         electronHOverE[4];
  Float_t         electronEcalE[4];
  Float_t         electronTrackP[4];
  Int_t           electronNMissingHits[4];
  Int_t           electronTripleCharge[4];
  Float_t         sf_zz;
  Float_t         sf_zzUnc;
  Float_t         sf_wz;
  Float_t         sf_zh;
  Float_t         sf_zhUp;
  Float_t         sf_zhDown;
  Float_t         genLep1Pt;
  Float_t         genLep1Eta;
  Float_t         genLep1Phi;
  Int_t           genLep1PdgId;
  Float_t         genLep2Pt;
  Float_t         genLep2Eta;
  Float_t         genLep2Phi;
  Int_t           genLep2PdgId;
  Int_t           looseGenLep1PdgId;
  Int_t           looseGenLep2PdgId;
  Int_t           looseGenLep3PdgId;
  Int_t           looseGenLep4PdgId;
  Int_t           whichRecoil;
  Float_t         genJet1Pt;
  Float_t         genJet2Pt;
  Float_t         genJet1Eta;
  Float_t         genJet2Eta;
  Float_t         genMjj;
  Int_t           badECALFilter;
  Int_t           jetNMBtags;
  Float_t         pfmetRaw;
  Int_t           nAK8jet;
  Int_t           nB;
  Float_t         fj1MSDScaleUp_sj;
  Float_t         fj1MSDScaleDown_sj;
  Float_t         fj1MSDSmeared_sj;
  Float_t         fj1MSDSmearedUp_sj;
  Float_t         fj1MSDSmearedDown_sj;
  Float_t         fj1PtScaleUp_sj;
  Float_t         fj1PtScaleDown_sj;
  Float_t         fj1PtSmeared_sj;
  Float_t         fj1PtSmearedUp_sj;
  Float_t         fj1PtSmearedDown_sj;
  Int_t           isGS;
  Float_t         fj1SubMaxCSV;
  Int_t           runNumber;
  Int_t           lumiNumber;
  ULong64_t       eventNumber;
  Int_t           npv;
  Int_t           pu;
  Float_t         mcWeight;
  Int_t           trigger;
  Int_t           metFilter;
  Int_t           egmFilter;
  Float_t         filter_maxRecoil;
  Float_t         filter_whichRecoil;
  Float_t         sf_ewkV;
  Float_t         sf_qcdV;
  Float_t         sf_ewkV2j;
  Float_t         sf_qcdV2j;
  Float_t         sf_qcdTT;
  Float_t         sf_lepID;
  Float_t         sf_lepIso;
  Float_t         sf_lepTrack;
  Float_t         sf_pho;
  Float_t         sf_eleTrig;
  Float_t         sf_phoTrig;
  Float_t         sf_metTrig;
  Float_t         sf_pu;
  Float_t         sf_npv;
  Float_t         sf_tt;
  Float_t         sf_tt_ext;
  Float_t         sf_tt_bound;
  Float_t         sf_tt8TeV;
  Float_t         sf_tt8TeV_ext;
  Float_t         sf_tt8TeV_bound;
  Float_t         sf_phoPurity;
  Float_t         pfmet;
  Float_t         pfmetphi;
  Float_t         pfmetnomu;
  Float_t         puppimet;
  Float_t         puppimetphi;
  Float_t         calomet;
  Float_t         calometphi;
  Float_t         pfcalobalance;
  Float_t         sumET;
  Float_t         trkmet;
  Float_t         puppiUWmag;
  Float_t         puppiUWphi;
  Float_t         puppiUZmag;
  Float_t         puppiUZphi;
  Float_t         puppiUAmag;
  Float_t         puppiUAphi;
  Float_t         puppiUperp;
  Float_t         puppiUpara;
  Float_t         puppiUmag;
  Float_t         puppiUphi;
  Float_t         pfUWmag;
  Float_t         pfUWphi;
  Float_t         pfUZmag;
  Float_t         pfUZphi;
  Float_t         pfUAmag;
  Float_t         pfUAphi;
  Float_t         pfUperp;
  Float_t         pfUpara;
  Float_t         pfUmag;
  Float_t         pfUphi;
  Float_t         dphipfmet;
  Float_t         dphipuppimet;
  Float_t         dphipuppiUW;
  Float_t         dphipuppiUZ;
  Float_t         dphipuppiUA;
  Float_t         dphipfUW;
  Float_t         dphipfUZ;
  Float_t         dphipfUA;
  Float_t         dphipuppiU;
  Float_t         dphipfU;
  Float_t         trueGenBosonPt;
  Float_t         genBosonPt;
  Float_t         genBosonEta;
  Float_t         genBosonMass;
  Float_t         genBosonPhi;
  Float_t         genWPlusPt;
  Float_t         genWMinusPt;
  Float_t         genWPlusEta;
  Float_t         genWMinusEta;
  Float_t         genTopPt;
  Int_t           genTopIsHad;
  Float_t         genTopEta;
  Float_t         genAntiTopPt;
  Int_t           genAntiTopIsHad;
  Float_t         genAntiTopEta;
  Float_t         genTTPt;
  Float_t         genTTEta;
  Int_t           nIsoJet;
  Int_t           jet1Flav;
  Float_t         jet1Phi;
  Float_t         jet1Pt;
  Float_t         jet1GenPt;
  Float_t         jet1Eta;
  Float_t         jet1CSV;
  Int_t           jet1IsTight;
  Int_t           jet2Flav;
  Float_t         jet2Phi;
  Float_t         jet2Pt;
  Float_t         jet2GenPt;
  Float_t         jet2Eta;
  Float_t         jet2CSV;
  Float_t         isojet1Pt;
  Float_t         isojet1CSV;
  Int_t           isojet1Flav;
  Float_t         isojet2Pt;
  Float_t         isojet2CSV;
  Int_t           isojet2Flav;
  Int_t           jetNBtags;
  Int_t           isojetNBtags;
  Int_t           nFatjet;
  Float_t         fj1Tau32;
  Float_t         fj1Tau21;
  Float_t         fj1Tau32SD;
  Float_t         fj1Tau21SD;
  Float_t         fj1MSD;
  Float_t         fj1MSDScaleUp;
  Float_t         fj1MSDScaleDown;
  Float_t         fj1MSDSmeared;
  Float_t         fj1MSDSmearedUp;
  Float_t         fj1MSDSmearedDown;
  Float_t         fj1MSD_corr;
  Float_t         fj1Pt;
  Float_t         fj1PtScaleUp;
  Float_t         fj1PtScaleDown;
  Float_t         fj1PtSmeared;
  Float_t         fj1PtSmearedUp;
  Float_t         fj1PtSmearedDown;
  Float_t         fj1Phi;
  Float_t         fj1Eta;
  Float_t         fj1M;
  Float_t         fj1MaxCSV;
  Float_t         fj1MinCSV;
  Float_t         fj1DoubleCSV;
  Float_t         fj1GenPt;
  Float_t         fj1GenSize;
  Int_t           fj1IsMatched;
  Float_t         fj1GenWPt;
  Float_t         fj1GenWSize;
  Int_t           fj1IsWMatched;
  Int_t           fj1HighestPtGen;
  Float_t         fj1HighestPtGenPt;
  Int_t           fj1IsTight;
  Int_t           fj1IsLoose;
  Float_t         fj1RawPt;
  Int_t           fj1NHF;
  Float_t         fj1HTTMass;
  Float_t         fj1HTTFRec;
  Int_t           fj1IsClean;
  Int_t           nHF;
  Int_t           nLoosePhoton;
  Int_t           nTightPhoton;
  Int_t           loosePho1IsTight;
  Float_t         loosePho1Pt;
  Float_t         loosePho1Eta;
  Float_t         loosePho1Phi;
  Float_t         diLepMass;
  Int_t           nTau;
  Float_t         mT;
  Float_t         scaleUp;
  Float_t         scaleDown;
  Float_t         pdfUp;
  Float_t         pdfDown;
  Float_t         normalizedWeight;

   // Book all branches explicitly as TBranches 
  TBranch *b_nJet                     = events->GetBranch("nJet")                     ; b_nJet                     ->SetAddress(&nJet                    );
  TBranch *b_nJot                     = events->GetBranch("nJot")                     ; b_nJot                     ->SetAddress(&nJot                    );
  TBranch *b_nLooseLep                = events->GetBranch("nLooseLep")                ; b_nLooseLep                ->SetAddress(&nLooseLep               );
  TBranch *b_nLooseElectron           = events->GetBranch("nLooseElectron")           ; b_nLooseElectron           ->SetAddress(&nLooseElectron          );
  TBranch *b_nLooseMuon               = events->GetBranch("nLooseMuon")               ; b_nLooseMuon               ->SetAddress(&nLooseMuon              );
  TBranch *b_nTightLep                = events->GetBranch("nTightLep")                ; b_nTightLep                ->SetAddress(&nTightLep               );
  TBranch *b_nTightElectron           = events->GetBranch("nTightElectron")           ; b_nTightElectron           ->SetAddress(&nTightElectron          );
  TBranch *b_nTightMuon               = events->GetBranch("nTightMuon")               ; b_nTightMuon               ->SetAddress(&nTightMuon              );
  TBranch *b_muonPt                   = events->GetBranch("muonPt")                   ; b_muonPt                   ->SetAddress(muonPt                   );
  TBranch *b_muonEta                  = events->GetBranch("muonEta")                  ; b_muonEta                  ->SetAddress(muonEta                  );
  TBranch *b_muonPhi                  = events->GetBranch("muonPhi")                  ; b_muonPhi                  ->SetAddress(muonPhi                  );
  TBranch *b_muonSfLoose              = events->GetBranch("muonSfLoose")              ; b_muonSfLoose              ->SetAddress(muonSfLoose              );
  TBranch *b_muonSfTight              = events->GetBranch("muonSfTight")              ; b_muonSfTight              ->SetAddress(muonSfTight              );
  TBranch *b_muonSfUnc                = events->GetBranch("muonSfUnc")                ; b_muonSfUnc                ->SetAddress(muonSfUnc                );
  TBranch *b_muonSfReco               = events->GetBranch("muonSfReco")               ; b_muonSfReco               ->SetAddress(muonSfReco               );
  TBranch *b_muonSelBit               = events->GetBranch("muonSelBit")               ; b_muonSelBit               ->SetAddress(muonSelBit               );
  TBranch *b_muonPdgId                = events->GetBranch("muonPdgId")                ; b_muonPdgId                ->SetAddress(muonPdgId                );
  TBranch *b_electronPt               = events->GetBranch("electronPt")               ; b_electronPt               ->SetAddress(electronPt               );
  TBranch *b_electronEta              = events->GetBranch("electronEta")              ; b_electronEta              ->SetAddress(electronEta              );
  TBranch *b_electronPhi              = events->GetBranch("electronPhi")              ; b_electronPhi              ->SetAddress(electronPhi              );
  TBranch *b_electronSfLoose          = events->GetBranch("electronSfLoose")          ; b_electronSfLoose          ->SetAddress(electronSfLoose          );
  TBranch *b_electronSfTight          = events->GetBranch("electronSfTight")          ; b_electronSfTight          ->SetAddress(electronSfTight          );
  TBranch *b_electronSfUnc            = events->GetBranch("electronSfUnc")            ; b_electronSfUnc            ->SetAddress(electronSfUnc            );
  TBranch *b_electronSfReco           = events->GetBranch("electronSfReco")           ; b_electronSfReco           ->SetAddress(electronSfReco           );
  TBranch *b_electronSelBit           = events->GetBranch("electronSelBit")           ; b_electronSelBit           ->SetAddress(electronSelBit           );
  TBranch *b_jetPt                    = events->GetBranch("jetPt")                    ; b_jetPt                    ->SetAddress(jetPt                    );
  TBranch *b_jetEta                   = events->GetBranch("jetEta")                   ; b_jetEta                   ->SetAddress(jetEta                   );
  TBranch *b_jetPhi                   = events->GetBranch("jetPhi")                   ; b_jetPhi                   ->SetAddress(jetPhi                   );
  TBranch *b_jetE                     = events->GetBranch("jetE")                     ; b_jetE                     ->SetAddress(jetE                     );
  TBranch *b_jetCSV                   = events->GetBranch("jetCSV")                   ; b_jetCSV                   ->SetAddress(jetCSV                   );
  TBranch *b_jetIso                   = events->GetBranch("jetIso")                   ; b_jetIso                   ->SetAddress(jetIso                   );
  TBranch *b_jetQGL                   = events->GetBranch("jetQGL")                   ; b_jetQGL                   ->SetAddress(jetQGL                   );
  TBranch *b_jetLeadingLepPt          = events->GetBranch("jetLeadingLepPt")          ; b_jetLeadingLepPt          ->SetAddress(jetLeadingLepPt          );
  TBranch *b_jetLeadingTrkPt          = events->GetBranch("jetLeadingTrkPt")          ; b_jetLeadingTrkPt          ->SetAddress(jetLeadingTrkPt          );
  TBranch *b_jetEMFrac                = events->GetBranch("jetEMFrac")                ; b_jetEMFrac                ->SetAddress(jetEMFrac                );
  TBranch *b_jetHadFrac               = events->GetBranch("jetHadFrac")               ; b_jetHadFrac               ->SetAddress(jetHadFrac               );
  TBranch *b_jetNLep                  = events->GetBranch("jetNLep")                  ; b_jetNLep                  ->SetAddress(jetNLep                  );
  TBranch *b_jetGenPt                 = events->GetBranch("jetGenPt")                 ; b_jetGenPt                 ->SetAddress(jetGenPt                 );
  TBranch *b_jetGenFlavor             = events->GetBranch("jetGenFlavor")             ; b_jetGenFlavor             ->SetAddress(jetGenFlavor             );
  TBranch *b_fj1sjPt                  = events->GetBranch("fj1sjPt")                  ; b_fj1sjPt                  ->SetAddress(fj1sjPt                  );
  TBranch *b_fj1sjPhi                 = events->GetBranch("fj1sjPhi")                 ; b_fj1sjPhi                 ->SetAddress(fj1sjPhi                 );
  TBranch *b_fj1sjEta                 = events->GetBranch("fj1sjEta")                 ; b_fj1sjEta                 ->SetAddress(fj1sjEta                 );
  TBranch *b_fj1sjM                   = events->GetBranch("fj1sjM")                   ; b_fj1sjM                   ->SetAddress(fj1sjM                   );
  TBranch *b_fj1sjCSV                 = events->GetBranch("fj1sjCSV")                 ; b_fj1sjCSV                 ->SetAddress(fj1sjCSV                 );
  TBranch *b_fj1sjQGL                 = events->GetBranch("fj1sjQGL")                 ; b_fj1sjQGL                 ->SetAddress(fj1sjQGL                 );
  TBranch *b_fj1Nbs                   = events->GetBranch("fj1Nbs")                   ; b_fj1Nbs                   ->SetAddress(&fj1Nbs                  );
  TBranch *b_fj1gbb                   = events->GetBranch("fj1gbb")                   ; b_fj1gbb                   ->SetAddress(&fj1gbb                  );
  TBranch *b_hbbpt                    = events->GetBranch("hbbpt")                    ; b_hbbpt                    ->SetAddress(&hbbpt                   );
  TBranch *b_hbbeta                   = events->GetBranch("hbbeta")                   ; b_hbbeta                   ->SetAddress(&hbbeta                  );
  TBranch *b_hbbphi                   = events->GetBranch("hbbphi")                   ; b_hbbphi                   ->SetAddress(&hbbphi                  );
  TBranch *b_hbbm                     = events->GetBranch("hbbm")                     ; b_hbbm                     ->SetAddress(&hbbm                    );
  TBranch *b_hbbjtidx                 = events->GetBranch("hbbjtidx")                 ; b_hbbjtidx                 ->SetAddress(hbbjtidx                 );
  TBranch *b_scale                    = events->GetBranch("scale")                    ; b_scale                    ->SetAddress(scale                    );
  TBranch *b_fj1ECFN_1_1_05           = events->GetBranch("fj1ECFN_1_1_05")           ; b_fj1ECFN_1_1_05           ->SetAddress(&fj1ECFN_1_1_05          );
  TBranch *b_fj1ECFN_2_1_05           = events->GetBranch("fj1ECFN_2_1_05")           ; b_fj1ECFN_2_1_05           ->SetAddress(&fj1ECFN_2_1_05          );
  TBranch *b_fj1ECFN_3_1_05           = events->GetBranch("fj1ECFN_3_1_05")           ; b_fj1ECFN_3_1_05           ->SetAddress(&fj1ECFN_3_1_05          );
  TBranch *b_fj1ECFN_1_2_05           = events->GetBranch("fj1ECFN_1_2_05")           ; b_fj1ECFN_1_2_05           ->SetAddress(&fj1ECFN_1_2_05          );
  TBranch *b_fj1ECFN_2_2_05           = events->GetBranch("fj1ECFN_2_2_05")           ; b_fj1ECFN_2_2_05           ->SetAddress(&fj1ECFN_2_2_05          );
  TBranch *b_fj1ECFN_3_2_05           = events->GetBranch("fj1ECFN_3_2_05")           ; b_fj1ECFN_3_2_05           ->SetAddress(&fj1ECFN_3_2_05          );
  TBranch *b_fj1ECFN_1_3_05           = events->GetBranch("fj1ECFN_1_3_05")           ; b_fj1ECFN_1_3_05           ->SetAddress(&fj1ECFN_1_3_05          );
  TBranch *b_fj1ECFN_2_3_05           = events->GetBranch("fj1ECFN_2_3_05")           ; b_fj1ECFN_2_3_05           ->SetAddress(&fj1ECFN_2_3_05          );
  TBranch *b_fj1ECFN_3_3_05           = events->GetBranch("fj1ECFN_3_3_05")           ; b_fj1ECFN_3_3_05           ->SetAddress(&fj1ECFN_3_3_05          );
  TBranch *b_fj1ECFN_1_4_05           = events->GetBranch("fj1ECFN_1_4_05")           ; b_fj1ECFN_1_4_05           ->SetAddress(&fj1ECFN_1_4_05          );
  TBranch *b_fj1ECFN_2_4_05           = events->GetBranch("fj1ECFN_2_4_05")           ; b_fj1ECFN_2_4_05           ->SetAddress(&fj1ECFN_2_4_05          );
  TBranch *b_fj1ECFN_3_4_05           = events->GetBranch("fj1ECFN_3_4_05")           ; b_fj1ECFN_3_4_05           ->SetAddress(&fj1ECFN_3_4_05          );
  TBranch *b_fj1ECFN_1_1_10           = events->GetBranch("fj1ECFN_1_1_10")           ; b_fj1ECFN_1_1_10           ->SetAddress(&fj1ECFN_1_1_10          );
  TBranch *b_fj1ECFN_2_1_10           = events->GetBranch("fj1ECFN_2_1_10")           ; b_fj1ECFN_2_1_10           ->SetAddress(&fj1ECFN_2_1_10          );
  TBranch *b_fj1ECFN_3_1_10           = events->GetBranch("fj1ECFN_3_1_10")           ; b_fj1ECFN_3_1_10           ->SetAddress(&fj1ECFN_3_1_10          );
  TBranch *b_fj1ECFN_1_2_10           = events->GetBranch("fj1ECFN_1_2_10")           ; b_fj1ECFN_1_2_10           ->SetAddress(&fj1ECFN_1_2_10          );
  TBranch *b_fj1ECFN_2_2_10           = events->GetBranch("fj1ECFN_2_2_10")           ; b_fj1ECFN_2_2_10           ->SetAddress(&fj1ECFN_2_2_10          );
  TBranch *b_fj1ECFN_3_2_10           = events->GetBranch("fj1ECFN_3_2_10")           ; b_fj1ECFN_3_2_10           ->SetAddress(&fj1ECFN_3_2_10          );
  TBranch *b_fj1ECFN_1_3_10           = events->GetBranch("fj1ECFN_1_3_10")           ; b_fj1ECFN_1_3_10           ->SetAddress(&fj1ECFN_1_3_10          );
  TBranch *b_fj1ECFN_2_3_10           = events->GetBranch("fj1ECFN_2_3_10")           ; b_fj1ECFN_2_3_10           ->SetAddress(&fj1ECFN_2_3_10          );
  TBranch *b_fj1ECFN_3_3_10           = events->GetBranch("fj1ECFN_3_3_10")           ; b_fj1ECFN_3_3_10           ->SetAddress(&fj1ECFN_3_3_10          );
  TBranch *b_fj1ECFN_1_4_10           = events->GetBranch("fj1ECFN_1_4_10")           ; b_fj1ECFN_1_4_10           ->SetAddress(&fj1ECFN_1_4_10          );
  TBranch *b_fj1ECFN_2_4_10           = events->GetBranch("fj1ECFN_2_4_10")           ; b_fj1ECFN_2_4_10           ->SetAddress(&fj1ECFN_2_4_10          );
  TBranch *b_fj1ECFN_3_4_10           = events->GetBranch("fj1ECFN_3_4_10")           ; b_fj1ECFN_3_4_10           ->SetAddress(&fj1ECFN_3_4_10          );
  TBranch *b_fj1ECFN_1_1_20           = events->GetBranch("fj1ECFN_1_1_20")           ; b_fj1ECFN_1_1_20           ->SetAddress(&fj1ECFN_1_1_20          );
  TBranch *b_fj1ECFN_2_1_20           = events->GetBranch("fj1ECFN_2_1_20")           ; b_fj1ECFN_2_1_20           ->SetAddress(&fj1ECFN_2_1_20          );
  TBranch *b_fj1ECFN_3_1_20           = events->GetBranch("fj1ECFN_3_1_20")           ; b_fj1ECFN_3_1_20           ->SetAddress(&fj1ECFN_3_1_20          );
  TBranch *b_fj1ECFN_1_2_20           = events->GetBranch("fj1ECFN_1_2_20")           ; b_fj1ECFN_1_2_20           ->SetAddress(&fj1ECFN_1_2_20          );
  TBranch *b_fj1ECFN_2_2_20           = events->GetBranch("fj1ECFN_2_2_20")           ; b_fj1ECFN_2_2_20           ->SetAddress(&fj1ECFN_2_2_20          );
  TBranch *b_fj1ECFN_3_2_20           = events->GetBranch("fj1ECFN_3_2_20")           ; b_fj1ECFN_3_2_20           ->SetAddress(&fj1ECFN_3_2_20          );
  TBranch *b_fj1ECFN_1_3_20           = events->GetBranch("fj1ECFN_1_3_20")           ; b_fj1ECFN_1_3_20           ->SetAddress(&fj1ECFN_1_3_20          );
  TBranch *b_fj1ECFN_2_3_20           = events->GetBranch("fj1ECFN_2_3_20")           ; b_fj1ECFN_2_3_20           ->SetAddress(&fj1ECFN_2_3_20          );
  TBranch *b_fj1ECFN_3_3_20           = events->GetBranch("fj1ECFN_3_3_20")           ; b_fj1ECFN_3_3_20           ->SetAddress(&fj1ECFN_3_3_20          );
  TBranch *b_fj1ECFN_1_4_20           = events->GetBranch("fj1ECFN_1_4_20")           ; b_fj1ECFN_1_4_20           ->SetAddress(&fj1ECFN_1_4_20          );
  TBranch *b_fj1ECFN_2_4_20           = events->GetBranch("fj1ECFN_2_4_20")           ; b_fj1ECFN_2_4_20           ->SetAddress(&fj1ECFN_2_4_20          );
  TBranch *b_fj1ECFN_3_4_20           = events->GetBranch("fj1ECFN_3_4_20")           ; b_fj1ECFN_3_4_20           ->SetAddress(&fj1ECFN_3_4_20          );
  TBranch *b_fj1ECFN_1_1_40           = events->GetBranch("fj1ECFN_1_1_40")           ; b_fj1ECFN_1_1_40           ->SetAddress(&fj1ECFN_1_1_40          );
  TBranch *b_fj1ECFN_2_1_40           = events->GetBranch("fj1ECFN_2_1_40")           ; b_fj1ECFN_2_1_40           ->SetAddress(&fj1ECFN_2_1_40          );
  TBranch *b_fj1ECFN_3_1_40           = events->GetBranch("fj1ECFN_3_1_40")           ; b_fj1ECFN_3_1_40           ->SetAddress(&fj1ECFN_3_1_40          );
  TBranch *b_fj1ECFN_1_2_40           = events->GetBranch("fj1ECFN_1_2_40")           ; b_fj1ECFN_1_2_40           ->SetAddress(&fj1ECFN_1_2_40          );
  TBranch *b_fj1ECFN_2_2_40           = events->GetBranch("fj1ECFN_2_2_40")           ; b_fj1ECFN_2_2_40           ->SetAddress(&fj1ECFN_2_2_40          );
  TBranch *b_fj1ECFN_3_2_40           = events->GetBranch("fj1ECFN_3_2_40")           ; b_fj1ECFN_3_2_40           ->SetAddress(&fj1ECFN_3_2_40          );
  TBranch *b_fj1ECFN_1_3_40           = events->GetBranch("fj1ECFN_1_3_40")           ; b_fj1ECFN_1_3_40           ->SetAddress(&fj1ECFN_1_3_40          );
  TBranch *b_fj1ECFN_2_3_40           = events->GetBranch("fj1ECFN_2_3_40")           ; b_fj1ECFN_2_3_40           ->SetAddress(&fj1ECFN_2_3_40          );
  TBranch *b_fj1ECFN_3_3_40           = events->GetBranch("fj1ECFN_3_3_40")           ; b_fj1ECFN_3_3_40           ->SetAddress(&fj1ECFN_3_3_40          );
  TBranch *b_fj1ECFN_1_4_40           = events->GetBranch("fj1ECFN_1_4_40")           ; b_fj1ECFN_1_4_40           ->SetAddress(&fj1ECFN_1_4_40          );
  TBranch *b_fj1ECFN_2_4_40           = events->GetBranch("fj1ECFN_2_4_40")           ; b_fj1ECFN_2_4_40           ->SetAddress(&fj1ECFN_2_4_40          );
  TBranch *b_fj1ECFN_3_4_40           = events->GetBranch("fj1ECFN_3_4_40")           ; b_fj1ECFN_3_4_40           ->SetAddress(&fj1ECFN_3_4_40          );
  TBranch *b_sf_btag0                 = events->GetBranch("sf_btag0")                 ; b_sf_btag0                 ->SetAddress(&sf_btag0                );
  TBranch *b_sf_btag1                 = events->GetBranch("sf_btag1")                 ; b_sf_btag1                 ->SetAddress(&sf_btag1                );
  TBranch *b_sf_btag2                 = events->GetBranch("sf_btag2")                 ; b_sf_btag2                 ->SetAddress(&sf_btag2                );
  TBranch *b_sf_btagGT0               = events->GetBranch("sf_btagGT0")               ; b_sf_btagGT0               ->SetAddress(&sf_btagGT0              );
  TBranch *b_sf_sjbtag0               = events->GetBranch("sf_sjbtag0")               ; b_sf_sjbtag0               ->SetAddress(&sf_sjbtag0              );
  TBranch *b_sf_sjbtag1               = events->GetBranch("sf_sjbtag1")               ; b_sf_sjbtag1               ->SetAddress(&sf_sjbtag1              );
  TBranch *b_sf_sjbtag2               = events->GetBranch("sf_sjbtag2")               ; b_sf_sjbtag2               ->SetAddress(&sf_sjbtag2              );
  TBranch *b_sf_sjbtagGT0             = events->GetBranch("sf_sjbtagGT0")             ; b_sf_sjbtagGT0             ->SetAddress(&sf_sjbtagGT0            );
  TBranch *b_sf_btag0BUp              = events->GetBranch("sf_btag0BUp")              ; b_sf_btag0BUp              ->SetAddress(&sf_btag0BUp             );
  TBranch *b_sf_btag1BUp              = events->GetBranch("sf_btag1BUp")              ; b_sf_btag1BUp              ->SetAddress(&sf_btag1BUp             );
  TBranch *b_sf_btag2BUp              = events->GetBranch("sf_btag2BUp")              ; b_sf_btag2BUp              ->SetAddress(&sf_btag2BUp             );
  TBranch *b_sf_btagGT0BUp            = events->GetBranch("sf_btagGT0BUp")            ; b_sf_btagGT0BUp            ->SetAddress(&sf_btagGT0BUp           );
  TBranch *b_sf_sjbtag0BUp            = events->GetBranch("sf_sjbtag0BUp")            ; b_sf_sjbtag0BUp            ->SetAddress(&sf_sjbtag0BUp           );
  TBranch *b_sf_sjbtag1BUp            = events->GetBranch("sf_sjbtag1BUp")            ; b_sf_sjbtag1BUp            ->SetAddress(&sf_sjbtag1BUp           );
  TBranch *b_sf_sjbtag2BUp            = events->GetBranch("sf_sjbtag2BUp")            ; b_sf_sjbtag2BUp            ->SetAddress(&sf_sjbtag2BUp           );
  TBranch *b_sf_sjbtagGT0BUp          = events->GetBranch("sf_sjbtagGT0BUp")          ; b_sf_sjbtagGT0BUp          ->SetAddress(&sf_sjbtagGT0BUp         );
  TBranch *b_sf_btag0BDown            = events->GetBranch("sf_btag0BDown")            ; b_sf_btag0BDown            ->SetAddress(&sf_btag0BDown           );
  TBranch *b_sf_btag1BDown            = events->GetBranch("sf_btag1BDown")            ; b_sf_btag1BDown            ->SetAddress(&sf_btag1BDown           );
  TBranch *b_sf_btag2BDown            = events->GetBranch("sf_btag2BDown")            ; b_sf_btag2BDown            ->SetAddress(&sf_btag2BDown           );
  TBranch *b_sf_btagGT0BDown          = events->GetBranch("sf_btagGT0BDown")          ; b_sf_btagGT0BDown          ->SetAddress(&sf_btagGT0BDown         );
  TBranch *b_sf_sjbtag0BDown          = events->GetBranch("sf_sjbtag0BDown")          ; b_sf_sjbtag0BDown          ->SetAddress(&sf_sjbtag0BDown         );
  TBranch *b_sf_sjbtag1BDown          = events->GetBranch("sf_sjbtag1BDown")          ; b_sf_sjbtag1BDown          ->SetAddress(&sf_sjbtag1BDown         );
  TBranch *b_sf_sjbtag2BDown          = events->GetBranch("sf_sjbtag2BDown")          ; b_sf_sjbtag2BDown          ->SetAddress(&sf_sjbtag2BDown         );
  TBranch *b_sf_sjbtagGT0BDown        = events->GetBranch("sf_sjbtagGT0BDown")        ; b_sf_sjbtagGT0BDown        ->SetAddress(&sf_sjbtagGT0BDown       );
  TBranch *b_sf_btag0MUp              = events->GetBranch("sf_btag0MUp")              ; b_sf_btag0MUp              ->SetAddress(&sf_btag0MUp             );
  TBranch *b_sf_btag1MUp              = events->GetBranch("sf_btag1MUp")              ; b_sf_btag1MUp              ->SetAddress(&sf_btag1MUp             );
  TBranch *b_sf_btag2MUp              = events->GetBranch("sf_btag2MUp")              ; b_sf_btag2MUp              ->SetAddress(&sf_btag2MUp             );
  TBranch *b_sf_btagGT0MUp            = events->GetBranch("sf_btagGT0MUp")            ; b_sf_btagGT0MUp            ->SetAddress(&sf_btagGT0MUp           );
  TBranch *b_sf_sjbtag0MUp            = events->GetBranch("sf_sjbtag0MUp")            ; b_sf_sjbtag0MUp            ->SetAddress(&sf_sjbtag0MUp           );
  TBranch *b_sf_sjbtag1MUp            = events->GetBranch("sf_sjbtag1MUp")            ; b_sf_sjbtag1MUp            ->SetAddress(&sf_sjbtag1MUp           );
  TBranch *b_sf_sjbtag2MUp            = events->GetBranch("sf_sjbtag2MUp")            ; b_sf_sjbtag2MUp            ->SetAddress(&sf_sjbtag2MUp           );
  TBranch *b_sf_sjbtagGT0MUp          = events->GetBranch("sf_sjbtagGT0MUp")          ; b_sf_sjbtagGT0MUp          ->SetAddress(&sf_sjbtagGT0MUp         );
  TBranch *b_sf_btag0MDown            = events->GetBranch("sf_btag0MDown")            ; b_sf_btag0MDown            ->SetAddress(&sf_btag0MDown           );
  TBranch *b_sf_btag1MDown            = events->GetBranch("sf_btag1MDown")            ; b_sf_btag1MDown            ->SetAddress(&sf_btag1MDown           );
  TBranch *b_sf_btag2MDown            = events->GetBranch("sf_btag2MDown")            ; b_sf_btag2MDown            ->SetAddress(&sf_btag2MDown           );
  TBranch *b_sf_btagGT0MDown          = events->GetBranch("sf_btagGT0MDown")          ; b_sf_btagGT0MDown          ->SetAddress(&sf_btagGT0MDown         );
  TBranch *b_sf_sjbtag0MDown          = events->GetBranch("sf_sjbtag0MDown")          ; b_sf_sjbtag0MDown          ->SetAddress(&sf_sjbtag0MDown         );
  TBranch *b_sf_sjbtag1MDown          = events->GetBranch("sf_sjbtag1MDown")          ; b_sf_sjbtag1MDown          ->SetAddress(&sf_sjbtag1MDown         );
  TBranch *b_sf_sjbtag2MDown          = events->GetBranch("sf_sjbtag2MDown")          ; b_sf_sjbtag2MDown          ->SetAddress(&sf_sjbtag2MDown         );
  TBranch *b_sf_sjbtagGT0MDown        = events->GetBranch("sf_sjbtagGT0MDown")        ; b_sf_sjbtagGT0MDown        ->SetAddress(&sf_sjbtagGT0MDown       );
  TBranch *b_muonD0                   = events->GetBranch("muonD0")                   ; b_muonD0                   ->SetAddress(muonD0                   );
  TBranch *b_muonDZ                   = events->GetBranch("muonDZ")                   ; b_muonDZ                   ->SetAddress(muonDZ                   );
  TBranch *b_muonSfMedium             = events->GetBranch("muonSfMedium")             ; b_muonSfMedium             ->SetAddress(muonSfMedium             );
  TBranch *b_muonIsSoftMuon           = events->GetBranch("muonIsSoftMuon")           ; b_muonIsSoftMuon           ->SetAddress(muonIsSoftMuon           );
  TBranch *b_muonIsGlobalMuon         = events->GetBranch("muonIsGlobalMuon")         ; b_muonIsGlobalMuon         ->SetAddress(muonIsGlobalMuon         );
  TBranch *b_muonIsTrackerMuon        = events->GetBranch("muonIsTrackerMuon")        ; b_muonIsTrackerMuon        ->SetAddress(muonIsTrackerMuon        );
  TBranch *b_muonNValidMuon           = events->GetBranch("muonNValidMuon")           ; b_muonNValidMuon           ->SetAddress(muonNValidMuon           );
  TBranch *b_muonNValidPixel          = events->GetBranch("muonNValidPixel")          ; b_muonNValidPixel          ->SetAddress(muonNValidPixel          );
  TBranch *b_muonTrkLayersWithMmt     = events->GetBranch("muonTrkLayersWithMmt")     ; b_muonTrkLayersWithMmt     ->SetAddress(muonTrkLayersWithMmt     );
  TBranch *b_muonPixLayersWithMmt     = events->GetBranch("muonPixLayersWithMmt")     ; b_muonPixLayersWithMmt     ->SetAddress(muonPixLayersWithMmt     );
  TBranch *b_muonNMatched             = events->GetBranch("muonNMatched")             ; b_muonNMatched             ->SetAddress(muonNMatched             );
  TBranch *b_muonChi2LocalPosition    = events->GetBranch("muonChi2LocalPosition")    ; b_muonChi2LocalPosition    ->SetAddress(muonChi2LocalPosition    );
  TBranch *b_muonTrkKink              = events->GetBranch("muonTrkKink")              ; b_muonTrkKink              ->SetAddress(muonTrkKink              );
  TBranch *b_muonValidFraction        = events->GetBranch("muonValidFraction")        ; b_muonValidFraction        ->SetAddress(muonValidFraction        );
  TBranch *b_muonNormChi2             = events->GetBranch("muonNormChi2")             ; b_muonNormChi2             ->SetAddress(muonNormChi2             );
  TBranch *b_muonSegmentCompatibility = events->GetBranch("muonSegmentCompatibility") ; b_muonSegmentCompatibility ->SetAddress(muonSegmentCompatibility );
  TBranch *b_electronPdgId            = events->GetBranch("electronPdgId")            ; b_electronPdgId            ->SetAddress(electronPdgId            );
  TBranch *b_electronD0               = events->GetBranch("electronD0")               ; b_electronD0               ->SetAddress(electronD0               );
  TBranch *b_electronDZ               = events->GetBranch("electronDZ")               ; b_electronDZ               ->SetAddress(electronDZ               );
  TBranch *b_electronSfMedium         = events->GetBranch("electronSfMedium")         ; b_electronSfMedium         ->SetAddress(electronSfMedium         );
  TBranch *b_electronChIsoPh          = events->GetBranch("electronChIsoPh")          ; b_electronChIsoPh          ->SetAddress(electronChIsoPh          );
  TBranch *b_electronNhIsoPh          = events->GetBranch("electronNhIsoPh")          ; b_electronNhIsoPh          ->SetAddress(electronNhIsoPh          );
  TBranch *b_electronPhIsoPh          = events->GetBranch("electronPhIsoPh")          ; b_electronPhIsoPh          ->SetAddress(electronPhIsoPh          );
  TBranch *b_electronEcalIso          = events->GetBranch("electronEcalIso")          ; b_electronEcalIso          ->SetAddress(electronEcalIso          );
  TBranch *b_electronHcalIso          = events->GetBranch("electronHcalIso")          ; b_electronHcalIso          ->SetAddress(electronHcalIso          );
  TBranch *b_electronTrackIso         = events->GetBranch("electronTrackIso")         ; b_electronTrackIso         ->SetAddress(electronTrackIso         );
  TBranch *b_electronIsoPUOffset      = events->GetBranch("electronIsoPUOffset")      ; b_electronIsoPUOffset      ->SetAddress(electronIsoPUOffset      );
  TBranch *b_electronSieie            = events->GetBranch("electronSieie")            ; b_electronSieie            ->SetAddress(electronSieie            );
  TBranch *b_electronSipip            = events->GetBranch("electronSipip")            ; b_electronSipip            ->SetAddress(electronSipip            );
  TBranch *b_electronDEtaInSeed       = events->GetBranch("electronDEtaInSeed")       ; b_electronDEtaInSeed       ->SetAddress(electronDEtaInSeed       );
  TBranch *b_electronDPhiIn           = events->GetBranch("electronDPhiIn")           ; b_electronDPhiIn           ->SetAddress(electronDPhiIn           );
  TBranch *b_electronEseed            = events->GetBranch("electronEseed")            ; b_electronEseed            ->SetAddress(electronEseed            );
  TBranch *b_electronHOverE           = events->GetBranch("electronHOverE")           ; b_electronHOverE           ->SetAddress(electronHOverE           );
  TBranch *b_electronEcalE            = events->GetBranch("electronEcalE")            ; b_electronEcalE            ->SetAddress(electronEcalE            );
  TBranch *b_electronTrackP           = events->GetBranch("electronTrackP")           ; b_electronTrackP           ->SetAddress(electronTrackP           );
  TBranch *b_electronNMissingHits     = events->GetBranch("electronNMissingHits")     ; b_electronNMissingHits     ->SetAddress(electronNMissingHits     );
  TBranch *b_electronTripleCharge     = events->GetBranch("electronTripleCharge")     ; b_electronTripleCharge     ->SetAddress(electronTripleCharge     );
  TBranch *b_sf_zz                    = events->GetBranch("sf_zz")                    ; b_sf_zz                    ->SetAddress(&sf_zz                   );
  TBranch *b_sf_zzUnc                 = events->GetBranch("sf_zzUnc")                 ; b_sf_zzUnc                 ->SetAddress(&sf_zzUnc                );
  TBranch *b_sf_wz                    = events->GetBranch("sf_wz")                    ; b_sf_wz                    ->SetAddress(&sf_wz                   );
  TBranch *b_sf_zh                    = events->GetBranch("sf_zh")                    ; b_sf_zh                    ->SetAddress(&sf_zh                   );
  TBranch *b_sf_zhUp                  = events->GetBranch("sf_zhUp")                  ; b_sf_zhUp                  ->SetAddress(&sf_zhUp                 );
  TBranch *b_sf_zhDown                = events->GetBranch("sf_zhDown")                ; b_sf_zhDown                ->SetAddress(&sf_zhDown               );
  TBranch *b_genLep1Pt                = events->GetBranch("genLep1Pt")                ; b_genLep1Pt                ->SetAddress(&genLep1Pt               );
  TBranch *b_genLep1Eta               = events->GetBranch("genLep1Eta")               ; b_genLep1Eta               ->SetAddress(&genLep1Eta              );
  TBranch *b_genLep1Phi               = events->GetBranch("genLep1Phi")               ; b_genLep1Phi               ->SetAddress(&genLep1Phi              );
  TBranch *b_genLep1PdgId             = events->GetBranch("genLep1PdgId")             ; b_genLep1PdgId             ->SetAddress(&genLep1PdgId            );
  TBranch *b_genLep2Pt                = events->GetBranch("genLep2Pt")                ; b_genLep2Pt                ->SetAddress(&genLep2Pt               );
  TBranch *b_genLep2Eta               = events->GetBranch("genLep2Eta")               ; b_genLep2Eta               ->SetAddress(&genLep2Eta              );
  TBranch *b_genLep2Phi               = events->GetBranch("genLep2Phi")               ; b_genLep2Phi               ->SetAddress(&genLep2Phi              );
  TBranch *b_genLep2PdgId             = events->GetBranch("genLep2PdgId")             ; b_genLep2PdgId             ->SetAddress(&genLep2PdgId            );
  TBranch *b_looseGenLep1PdgId        = events->GetBranch("looseGenLep1PdgId")        ; b_looseGenLep1PdgId        ->SetAddress(&looseGenLep1PdgId       );
  TBranch *b_looseGenLep2PdgId        = events->GetBranch("looseGenLep2PdgId")        ; b_looseGenLep2PdgId        ->SetAddress(&looseGenLep2PdgId       );
  TBranch *b_looseGenLep3PdgId        = events->GetBranch("looseGenLep3PdgId")        ; b_looseGenLep3PdgId        ->SetAddress(&looseGenLep3PdgId       );
  TBranch *b_looseGenLep4PdgId        = events->GetBranch("looseGenLep4PdgId")        ; b_looseGenLep4PdgId        ->SetAddress(&looseGenLep4PdgId       );
  TBranch *b_whichRecoil              = events->GetBranch("whichRecoil")              ; b_whichRecoil              ->SetAddress(&whichRecoil             );
  TBranch *b_genJet1Pt                = events->GetBranch("genJet1Pt")                ; b_genJet1Pt                ->SetAddress(&genJet1Pt               );
  TBranch *b_genJet2Pt                = events->GetBranch("genJet2Pt")                ; b_genJet2Pt                ->SetAddress(&genJet2Pt               );
  TBranch *b_genJet1Eta               = events->GetBranch("genJet1Eta")               ; b_genJet1Eta               ->SetAddress(&genJet1Eta              );
  TBranch *b_genJet2Eta               = events->GetBranch("genJet2Eta")               ; b_genJet2Eta               ->SetAddress(&genJet2Eta              );
  TBranch *b_genMjj                   = events->GetBranch("genMjj")                   ; b_genMjj                   ->SetAddress(&genMjj                  );
  TBranch *b_badECALFilter            = events->GetBranch("badECALFilter")            ; b_badECALFilter            ->SetAddress(&badECALFilter           );
  TBranch *b_jetNMBtags               = events->GetBranch("jetNMBtags")               ; b_jetNMBtags               ->SetAddress(&jetNMBtags              );
  TBranch *b_pfmetRaw                 = events->GetBranch("pfmetRaw")                 ; b_pfmetRaw                 ->SetAddress(&pfmetRaw                );
  TBranch *b_nAK8jet                  = events->GetBranch("nAK8jet")                  ; b_nAK8jet                  ->SetAddress(&nAK8jet                 );
  TBranch *b_nB                       = events->GetBranch("nB")                       ; b_nB                       ->SetAddress(&nB                      );
  TBranch *b_fj1MSDScaleUp_sj         = events->GetBranch("fj1MSDScaleUp_sj")         ; b_fj1MSDScaleUp_sj         ->SetAddress(&fj1MSDScaleUp_sj        );
  TBranch *b_fj1MSDScaleDown_sj       = events->GetBranch("fj1MSDScaleDown_sj")       ; b_fj1MSDScaleDown_sj       ->SetAddress(&fj1MSDScaleDown_sj      );
  TBranch *b_fj1MSDSmeared_sj         = events->GetBranch("fj1MSDSmeared_sj")         ; b_fj1MSDSmeared_sj         ->SetAddress(&fj1MSDSmeared_sj        );
  TBranch *b_fj1MSDSmearedUp_sj       = events->GetBranch("fj1MSDSmearedUp_sj")       ; b_fj1MSDSmearedUp_sj       ->SetAddress(&fj1MSDSmearedUp_sj      );
  TBranch *b_fj1MSDSmearedDown_sj     = events->GetBranch("fj1MSDSmearedDown_sj")     ; b_fj1MSDSmearedDown_sj     ->SetAddress(&fj1MSDSmearedDown_sj    );
  TBranch *b_fj1PtScaleUp_sj          = events->GetBranch("fj1PtScaleUp_sj")          ; b_fj1PtScaleUp_sj          ->SetAddress(&fj1PtScaleUp_sj         );
  TBranch *b_fj1PtScaleDown_sj        = events->GetBranch("fj1PtScaleDown_sj")        ; b_fj1PtScaleDown_sj        ->SetAddress(&fj1PtScaleDown_sj       );
  TBranch *b_fj1PtSmeared_sj          = events->GetBranch("fj1PtSmeared_sj")          ; b_fj1PtSmeared_sj          ->SetAddress(&fj1PtSmeared_sj         );
  TBranch *b_fj1PtSmearedUp_sj        = events->GetBranch("fj1PtSmearedUp_sj")        ; b_fj1PtSmearedUp_sj        ->SetAddress(&fj1PtSmearedUp_sj       );
  TBranch *b_fj1PtSmearedDown_sj      = events->GetBranch("fj1PtSmearedDown_sj")      ; b_fj1PtSmearedDown_sj      ->SetAddress(&fj1PtSmearedDown_sj     );
  TBranch *b_isGS                     = events->GetBranch("isGS")                     ; b_isGS                     ->SetAddress(&isGS                    );
  TBranch *b_fj1SubMaxCSV             = events->GetBranch("fj1SubMaxCSV")             ; b_fj1SubMaxCSV             ->SetAddress(&fj1SubMaxCSV            );
  TBranch *b_runNumber                = events->GetBranch("runNumber")                ; b_runNumber                ->SetAddress(&runNumber               );
  TBranch *b_lumiNumber               = events->GetBranch("lumiNumber")               ; b_lumiNumber               ->SetAddress(&lumiNumber              );
  TBranch *b_eventNumber              = events->GetBranch("eventNumber")              ; b_eventNumber              ->SetAddress(&eventNumber             );
  TBranch *b_npv                      = events->GetBranch("npv")                      ; b_npv                      ->SetAddress(&npv                     );
  TBranch *b_pu                       = events->GetBranch("pu")                       ; b_pu                       ->SetAddress(&pu                      );
  TBranch *b_mcWeight                 = events->GetBranch("mcWeight")                 ; b_mcWeight                 ->SetAddress(&mcWeight                );
  TBranch *b_trigger                  = events->GetBranch("trigger")                  ; b_trigger                  ->SetAddress(&trigger                 );
  TBranch *b_metFilter                = events->GetBranch("metFilter")                ; b_metFilter                ->SetAddress(&metFilter               );
  TBranch *b_egmFilter                = events->GetBranch("egmFilter")                ; b_egmFilter                ->SetAddress(&egmFilter               );
  TBranch *b_filter_maxRecoil         = events->GetBranch("filter_maxRecoil")         ; b_filter_maxRecoil         ->SetAddress(&filter_maxRecoil        );
  TBranch *b_filter_whichRecoil       = events->GetBranch("filter_whichRecoil")       ; b_filter_whichRecoil       ->SetAddress(&filter_whichRecoil      );
  TBranch *b_sf_ewkV                  = events->GetBranch("sf_ewkV")                  ; b_sf_ewkV                  ->SetAddress(&sf_ewkV                 );
  TBranch *b_sf_qcdV                  = events->GetBranch("sf_qcdV")                  ; b_sf_qcdV                  ->SetAddress(&sf_qcdV                 );
  TBranch *b_sf_ewkV2j                = events->GetBranch("sf_ewkV2j")                ; b_sf_ewkV2j                ->SetAddress(&sf_ewkV2j               );
  TBranch *b_sf_qcdV2j                = events->GetBranch("sf_qcdV2j")                ; b_sf_qcdV2j                ->SetAddress(&sf_qcdV2j               );
  TBranch *b_sf_qcdTT                 = events->GetBranch("sf_qcdTT")                 ; b_sf_qcdTT                 ->SetAddress(&sf_qcdTT                );
  TBranch *b_sf_lepID                 = events->GetBranch("sf_lepID")                 ; b_sf_lepID                 ->SetAddress(&sf_lepID                );
  TBranch *b_sf_lepIso                = events->GetBranch("sf_lepIso")                ; b_sf_lepIso                ->SetAddress(&sf_lepIso               );
  TBranch *b_sf_lepTrack              = events->GetBranch("sf_lepTrack")              ; b_sf_lepTrack              ->SetAddress(&sf_lepTrack             );
  TBranch *b_sf_pho                   = events->GetBranch("sf_pho")                   ; b_sf_pho                   ->SetAddress(&sf_pho                  );
  TBranch *b_sf_eleTrig               = events->GetBranch("sf_eleTrig")               ; b_sf_eleTrig               ->SetAddress(&sf_eleTrig              );
  TBranch *b_sf_phoTrig               = events->GetBranch("sf_phoTrig")               ; b_sf_phoTrig               ->SetAddress(&sf_phoTrig              );
  TBranch *b_sf_metTrig               = events->GetBranch("sf_metTrig")               ; b_sf_metTrig               ->SetAddress(&sf_metTrig              );
  TBranch *b_sf_pu                    = events->GetBranch("sf_pu")                    ; b_sf_pu                    ->SetAddress(&sf_pu                   );
  TBranch *b_sf_npv                   = events->GetBranch("sf_npv")                   ; b_sf_npv                   ->SetAddress(&sf_npv                  );
  TBranch *b_sf_tt                    = events->GetBranch("sf_tt")                    ; b_sf_tt                    ->SetAddress(&sf_tt                   );
  TBranch *b_sf_tt_ext                = events->GetBranch("sf_tt_ext")                ; b_sf_tt_ext                ->SetAddress(&sf_tt_ext               );
  TBranch *b_sf_tt_bound              = events->GetBranch("sf_tt_bound")              ; b_sf_tt_bound              ->SetAddress(&sf_tt_bound             );
  TBranch *b_sf_tt8TeV                = events->GetBranch("sf_tt8TeV")                ; b_sf_tt8TeV                ->SetAddress(&sf_tt8TeV               );
  TBranch *b_sf_tt8TeV_ext            = events->GetBranch("sf_tt8TeV_ext")            ; b_sf_tt8TeV_ext            ->SetAddress(&sf_tt8TeV_ext           );
  TBranch *b_sf_tt8TeV_bound          = events->GetBranch("sf_tt8TeV_bound")          ; b_sf_tt8TeV_bound          ->SetAddress(&sf_tt8TeV_bound         );
  TBranch *b_sf_phoPurity             = events->GetBranch("sf_phoPurity")             ; b_sf_phoPurity             ->SetAddress(&sf_phoPurity            );
  TBranch *b_pfmet                    = events->GetBranch("pfmet")                    ; b_pfmet                    ->SetAddress(&pfmet                   );
  TBranch *b_pfmetphi                 = events->GetBranch("pfmetphi")                 ; b_pfmetphi                 ->SetAddress(&pfmetphi                );
  TBranch *b_pfmetnomu                = events->GetBranch("pfmetnomu")                ; b_pfmetnomu                ->SetAddress(&pfmetnomu               );
  TBranch *b_puppimet                 = events->GetBranch("puppimet")                 ; b_puppimet                 ->SetAddress(&puppimet                );
  TBranch *b_puppimetphi              = events->GetBranch("puppimetphi")              ; b_puppimetphi              ->SetAddress(&puppimetphi             );
  TBranch *b_calomet                  = events->GetBranch("calomet")                  ; b_calomet                  ->SetAddress(&calomet                 );
  TBranch *b_calometphi               = events->GetBranch("calometphi")               ; b_calometphi               ->SetAddress(&calometphi              );
  TBranch *b_pfcalobalance            = events->GetBranch("pfcalobalance")            ; b_pfcalobalance            ->SetAddress(&pfcalobalance           );
  TBranch *b_sumET                    = events->GetBranch("sumET")                    ; b_sumET                    ->SetAddress(&sumET                   );
  TBranch *b_trkmet                   = events->GetBranch("trkmet")                   ; b_trkmet                   ->SetAddress(&trkmet                  );
  TBranch *b_puppiUWmag               = events->GetBranch("puppiUWmag")               ; b_puppiUWmag               ->SetAddress(&puppiUWmag              );
  TBranch *b_puppiUWphi               = events->GetBranch("puppiUWphi")               ; b_puppiUWphi               ->SetAddress(&puppiUWphi              );
  TBranch *b_puppiUZmag               = events->GetBranch("puppiUZmag")               ; b_puppiUZmag               ->SetAddress(&puppiUZmag              );
  TBranch *b_puppiUZphi               = events->GetBranch("puppiUZphi")               ; b_puppiUZphi               ->SetAddress(&puppiUZphi              );
  TBranch *b_puppiUAmag               = events->GetBranch("puppiUAmag")               ; b_puppiUAmag               ->SetAddress(&puppiUAmag              );
  TBranch *b_puppiUAphi               = events->GetBranch("puppiUAphi")               ; b_puppiUAphi               ->SetAddress(&puppiUAphi              );
  TBranch *b_puppiUperp               = events->GetBranch("puppiUperp")               ; b_puppiUperp               ->SetAddress(&puppiUperp              );
  TBranch *b_puppiUpara               = events->GetBranch("puppiUpara")               ; b_puppiUpara               ->SetAddress(&puppiUpara              );
  TBranch *b_puppiUmag                = events->GetBranch("puppiUmag")                ; b_puppiUmag                ->SetAddress(&puppiUmag               );
  TBranch *b_puppiUphi                = events->GetBranch("puppiUphi")                ; b_puppiUphi                ->SetAddress(&puppiUphi               );
  TBranch *b_pfUWmag                  = events->GetBranch("pfUWmag")                  ; b_pfUWmag                  ->SetAddress(&pfUWmag                 );
  TBranch *b_pfUWphi                  = events->GetBranch("pfUWphi")                  ; b_pfUWphi                  ->SetAddress(&pfUWphi                 );
  TBranch *b_pfUZmag                  = events->GetBranch("pfUZmag")                  ; b_pfUZmag                  ->SetAddress(&pfUZmag                 );
  TBranch *b_pfUZphi                  = events->GetBranch("pfUZphi")                  ; b_pfUZphi                  ->SetAddress(&pfUZphi                 );
  TBranch *b_pfUAmag                  = events->GetBranch("pfUAmag")                  ; b_pfUAmag                  ->SetAddress(&pfUAmag                 );
  TBranch *b_pfUAphi                  = events->GetBranch("pfUAphi")                  ; b_pfUAphi                  ->SetAddress(&pfUAphi                 );
  TBranch *b_pfUperp                  = events->GetBranch("pfUperp")                  ; b_pfUperp                  ->SetAddress(&pfUperp                 );
  TBranch *b_pfUpara                  = events->GetBranch("pfUpara")                  ; b_pfUpara                  ->SetAddress(&pfUpara                 );
  TBranch *b_pfUmag                   = events->GetBranch("pfUmag")                   ; b_pfUmag                   ->SetAddress(&pfUmag                  );
  TBranch *b_pfUphi                   = events->GetBranch("pfUphi")                   ; b_pfUphi                   ->SetAddress(&pfUphi                  );
  TBranch *b_dphipfmet                = events->GetBranch("dphipfmet")                ; b_dphipfmet                ->SetAddress(&dphipfmet               );
  TBranch *b_dphipuppimet             = events->GetBranch("dphipuppimet")             ; b_dphipuppimet             ->SetAddress(&dphipuppimet            );
  TBranch *b_dphipuppiUW              = events->GetBranch("dphipuppiUW")              ; b_dphipuppiUW              ->SetAddress(&dphipuppiUW             );
  TBranch *b_dphipuppiUZ              = events->GetBranch("dphipuppiUZ")              ; b_dphipuppiUZ              ->SetAddress(&dphipuppiUZ             );
  TBranch *b_dphipuppiUA              = events->GetBranch("dphipuppiUA")              ; b_dphipuppiUA              ->SetAddress(&dphipuppiUA             );
  TBranch *b_dphipfUW                 = events->GetBranch("dphipfUW")                 ; b_dphipfUW                 ->SetAddress(&dphipfUW                );
  TBranch *b_dphipfUZ                 = events->GetBranch("dphipfUZ")                 ; b_dphipfUZ                 ->SetAddress(&dphipfUZ                );
  TBranch *b_dphipfUA                 = events->GetBranch("dphipfUA")                 ; b_dphipfUA                 ->SetAddress(&dphipfUA                );
  TBranch *b_dphipuppiU               = events->GetBranch("dphipuppiU")               ; b_dphipuppiU               ->SetAddress(&dphipuppiU              );
  TBranch *b_dphipfU                  = events->GetBranch("dphipfU")                  ; b_dphipfU                  ->SetAddress(&dphipfU                 );
  TBranch *b_trueGenBosonPt           = events->GetBranch("trueGenBosonPt")           ; b_trueGenBosonPt           ->SetAddress(&trueGenBosonPt          );
  TBranch *b_genBosonPt               = events->GetBranch("genBosonPt")               ; b_genBosonPt               ->SetAddress(&genBosonPt              );
  TBranch *b_genBosonEta              = events->GetBranch("genBosonEta")              ; b_genBosonEta              ->SetAddress(&genBosonEta             );
  TBranch *b_genBosonMass             = events->GetBranch("genBosonMass")             ; b_genBosonMass             ->SetAddress(&genBosonMass            );
  TBranch *b_genBosonPhi              = events->GetBranch("genBosonPhi")              ; b_genBosonPhi              ->SetAddress(&genBosonPhi             );
  TBranch *b_genWPlusPt               = events->GetBranch("genWPlusPt")               ; b_genWPlusPt               ->SetAddress(&genWPlusPt              );
  TBranch *b_genWMinusPt              = events->GetBranch("genWMinusPt")              ; b_genWMinusPt              ->SetAddress(&genWMinusPt             );
  TBranch *b_genWPlusEta              = events->GetBranch("genWPlusEta")              ; b_genWPlusEta              ->SetAddress(&genWPlusEta             );
  TBranch *b_genWMinusEta             = events->GetBranch("genWMinusEta")             ; b_genWMinusEta             ->SetAddress(&genWMinusEta            );
  TBranch *b_genTopPt                 = events->GetBranch("genTopPt")                 ; b_genTopPt                 ->SetAddress(&genTopPt                );
  TBranch *b_genTopIsHad              = events->GetBranch("genTopIsHad")              ; b_genTopIsHad              ->SetAddress(&genTopIsHad             );
  TBranch *b_genTopEta                = events->GetBranch("genTopEta")                ; b_genTopEta                ->SetAddress(&genTopEta               );
  TBranch *b_genAntiTopPt             = events->GetBranch("genAntiTopPt")             ; b_genAntiTopPt             ->SetAddress(&genAntiTopPt            );
  TBranch *b_genAntiTopIsHad          = events->GetBranch("genAntiTopIsHad")          ; b_genAntiTopIsHad          ->SetAddress(&genAntiTopIsHad         );
  TBranch *b_genAntiTopEta            = events->GetBranch("genAntiTopEta")            ; b_genAntiTopEta            ->SetAddress(&genAntiTopEta           );
  TBranch *b_genTTPt                  = events->GetBranch("genTTPt")                  ; b_genTTPt                  ->SetAddress(&genTTPt                 );
  TBranch *b_genTTEta                 = events->GetBranch("genTTEta")                 ; b_genTTEta                 ->SetAddress(&genTTEta                );
  TBranch *b_nIsoJet                  = events->GetBranch("nIsoJet")                  ; b_nIsoJet                  ->SetAddress(&nIsoJet                 );
  TBranch *b_jet1Flav                 = events->GetBranch("jet1Flav")                 ; b_jet1Flav                 ->SetAddress(&jet1Flav                );
  TBranch *b_jet1Phi                  = events->GetBranch("jet1Phi")                  ; b_jet1Phi                  ->SetAddress(&jet1Phi                 );
  TBranch *b_jet1Pt                   = events->GetBranch("jet1Pt")                   ; b_jet1Pt                   ->SetAddress(&jet1Pt                  );
  TBranch *b_jet1GenPt                = events->GetBranch("jet1GenPt")                ; b_jet1GenPt                ->SetAddress(&jet1GenPt               );
  TBranch *b_jet1Eta                  = events->GetBranch("jet1Eta")                  ; b_jet1Eta                  ->SetAddress(&jet1Eta                 );
  TBranch *b_jet1CSV                  = events->GetBranch("jet1CSV")                  ; b_jet1CSV                  ->SetAddress(&jet1CSV                 );
  TBranch *b_jet1IsTight              = events->GetBranch("jet1IsTight")              ; b_jet1IsTight              ->SetAddress(&jet1IsTight             );
  TBranch *b_jet2Flav                 = events->GetBranch("jet2Flav")                 ; b_jet2Flav                 ->SetAddress(&jet2Flav                );
  TBranch *b_jet2Phi                  = events->GetBranch("jet2Phi")                  ; b_jet2Phi                  ->SetAddress(&jet2Phi                 );
  TBranch *b_jet2Pt                   = events->GetBranch("jet2Pt")                   ; b_jet2Pt                   ->SetAddress(&jet2Pt                  );
  TBranch *b_jet2GenPt                = events->GetBranch("jet2GenPt")                ; b_jet2GenPt                ->SetAddress(&jet2GenPt               );
  TBranch *b_jet2Eta                  = events->GetBranch("jet2Eta")                  ; b_jet2Eta                  ->SetAddress(&jet2Eta                 );
  TBranch *b_jet2CSV                  = events->GetBranch("jet2CSV")                  ; b_jet2CSV                  ->SetAddress(&jet2CSV                 );
  TBranch *b_isojet1Pt                = events->GetBranch("isojet1Pt")                ; b_isojet1Pt                ->SetAddress(&isojet1Pt               );
  TBranch *b_isojet1CSV               = events->GetBranch("isojet1CSV")               ; b_isojet1CSV               ->SetAddress(&isojet1CSV              );
  TBranch *b_isojet1Flav              = events->GetBranch("isojet1Flav")              ; b_isojet1Flav              ->SetAddress(&isojet1Flav             );
  TBranch *b_isojet2Pt                = events->GetBranch("isojet2Pt")                ; b_isojet2Pt                ->SetAddress(&isojet2Pt               );
  TBranch *b_isojet2CSV               = events->GetBranch("isojet2CSV")               ; b_isojet2CSV               ->SetAddress(&isojet2CSV              );
  TBranch *b_isojet2Flav              = events->GetBranch("isojet2Flav")              ; b_isojet2Flav              ->SetAddress(&isojet2Flav             );
  TBranch *b_jetNBtags                = events->GetBranch("jetNBtags")                ; b_jetNBtags                ->SetAddress(&jetNBtags               );
  TBranch *b_isojetNBtags             = events->GetBranch("isojetNBtags")             ; b_isojetNBtags             ->SetAddress(&isojetNBtags            );
  TBranch *b_nFatjet                  = events->GetBranch("nFatjet")                  ; b_nFatjet                  ->SetAddress(&nFatjet                 );
  TBranch *b_fj1Tau32                 = events->GetBranch("fj1Tau32")                 ; b_fj1Tau32                 ->SetAddress(&fj1Tau32                );
  TBranch *b_fj1Tau21                 = events->GetBranch("fj1Tau21")                 ; b_fj1Tau21                 ->SetAddress(&fj1Tau21                );
  TBranch *b_fj1Tau32SD               = events->GetBranch("fj1Tau32SD")               ; b_fj1Tau32SD               ->SetAddress(&fj1Tau32SD              );
  TBranch *b_fj1Tau21SD               = events->GetBranch("fj1Tau21SD")               ; b_fj1Tau21SD               ->SetAddress(&fj1Tau21SD              );
  TBranch *b_fj1MSD                   = events->GetBranch("fj1MSD")                   ; b_fj1MSD                   ->SetAddress(&fj1MSD                  );
  TBranch *b_fj1MSDScaleUp            = events->GetBranch("fj1MSDScaleUp")            ; b_fj1MSDScaleUp            ->SetAddress(&fj1MSDScaleUp           );
  TBranch *b_fj1MSDScaleDown          = events->GetBranch("fj1MSDScaleDown")          ; b_fj1MSDScaleDown          ->SetAddress(&fj1MSDScaleDown         );
  TBranch *b_fj1MSDSmeared            = events->GetBranch("fj1MSDSmeared")            ; b_fj1MSDSmeared            ->SetAddress(&fj1MSDSmeared           );
  TBranch *b_fj1MSDSmearedUp          = events->GetBranch("fj1MSDSmearedUp")          ; b_fj1MSDSmearedUp          ->SetAddress(&fj1MSDSmearedUp         );
  TBranch *b_fj1MSDSmearedDown        = events->GetBranch("fj1MSDSmearedDown")        ; b_fj1MSDSmearedDown        ->SetAddress(&fj1MSDSmearedDown       );
  TBranch *b_fj1MSD_corr              = events->GetBranch("fj1MSD_corr")              ; b_fj1MSD_corr              ->SetAddress(&fj1MSD_corr             );
  TBranch *b_fj1Pt                    = events->GetBranch("fj1Pt")                    ; b_fj1Pt                    ->SetAddress(&fj1Pt                   );
  TBranch *b_fj1PtScaleUp             = events->GetBranch("fj1PtScaleUp")             ; b_fj1PtScaleUp             ->SetAddress(&fj1PtScaleUp            );
  TBranch *b_fj1PtScaleDown           = events->GetBranch("fj1PtScaleDown")           ; b_fj1PtScaleDown           ->SetAddress(&fj1PtScaleDown          );
  TBranch *b_fj1PtSmeared             = events->GetBranch("fj1PtSmeared")             ; b_fj1PtSmeared             ->SetAddress(&fj1PtSmeared            );
  TBranch *b_fj1PtSmearedUp           = events->GetBranch("fj1PtSmearedUp")           ; b_fj1PtSmearedUp           ->SetAddress(&fj1PtSmearedUp          );
  TBranch *b_fj1PtSmearedDown         = events->GetBranch("fj1PtSmearedDown")         ; b_fj1PtSmearedDown         ->SetAddress(&fj1PtSmearedDown        );
  TBranch *b_fj1Phi                   = events->GetBranch("fj1Phi")                   ; b_fj1Phi                   ->SetAddress(&fj1Phi                  );
  TBranch *b_fj1Eta                   = events->GetBranch("fj1Eta")                   ; b_fj1Eta                   ->SetAddress(&fj1Eta                  );
  TBranch *b_fj1M                     = events->GetBranch("fj1M")                     ; b_fj1M                     ->SetAddress(&fj1M                    );
  TBranch *b_fj1MaxCSV                = events->GetBranch("fj1MaxCSV")                ; b_fj1MaxCSV                ->SetAddress(&fj1MaxCSV               );
  TBranch *b_fj1MinCSV                = events->GetBranch("fj1MinCSV")                ; b_fj1MinCSV                ->SetAddress(&fj1MinCSV               );
  TBranch *b_fj1DoubleCSV             = events->GetBranch("fj1DoubleCSV")             ; b_fj1DoubleCSV             ->SetAddress(&fj1DoubleCSV            );
  TBranch *b_fj1GenPt                 = events->GetBranch("fj1GenPt")                 ; b_fj1GenPt                 ->SetAddress(&fj1GenPt                );
  TBranch *b_fj1GenSize               = events->GetBranch("fj1GenSize")               ; b_fj1GenSize               ->SetAddress(&fj1GenSize              );
  TBranch *b_fj1IsMatched             = events->GetBranch("fj1IsMatched")             ; b_fj1IsMatched             ->SetAddress(&fj1IsMatched            );
  TBranch *b_fj1GenWPt                = events->GetBranch("fj1GenWPt")                ; b_fj1GenWPt                ->SetAddress(&fj1GenWPt               );
  TBranch *b_fj1GenWSize              = events->GetBranch("fj1GenWSize")              ; b_fj1GenWSize              ->SetAddress(&fj1GenWSize             );
  TBranch *b_fj1IsWMatched            = events->GetBranch("fj1IsWMatched")            ; b_fj1IsWMatched            ->SetAddress(&fj1IsWMatched           );
  TBranch *b_fj1HighestPtGen          = events->GetBranch("fj1HighestPtGen")          ; b_fj1HighestPtGen          ->SetAddress(&fj1HighestPtGen         );
  TBranch *b_fj1HighestPtGenPt        = events->GetBranch("fj1HighestPtGenPt")        ; b_fj1HighestPtGenPt        ->SetAddress(&fj1HighestPtGenPt       );
  TBranch *b_fj1IsTight               = events->GetBranch("fj1IsTight")               ; b_fj1IsTight               ->SetAddress(&fj1IsTight              );
  TBranch *b_fj1IsLoose               = events->GetBranch("fj1IsLoose")               ; b_fj1IsLoose               ->SetAddress(&fj1IsLoose              );
  TBranch *b_fj1RawPt                 = events->GetBranch("fj1RawPt")                 ; b_fj1RawPt                 ->SetAddress(&fj1RawPt                );
  TBranch *b_fj1NHF                   = events->GetBranch("fj1NHF")                   ; b_fj1NHF                   ->SetAddress(&fj1NHF                  );
  TBranch *b_fj1HTTMass               = events->GetBranch("fj1HTTMass")               ; b_fj1HTTMass               ->SetAddress(&fj1HTTMass              );
  TBranch *b_fj1HTTFRec               = events->GetBranch("fj1HTTFRec")               ; b_fj1HTTFRec               ->SetAddress(&fj1HTTFRec              );
  TBranch *b_fj1IsClean               = events->GetBranch("fj1IsClean")               ; b_fj1IsClean               ->SetAddress(&fj1IsClean              );
  TBranch *b_nHF                      = events->GetBranch("nHF")                      ; b_nHF                      ->SetAddress(&nHF                     );
  TBranch *b_nLoosePhoton             = events->GetBranch("nLoosePhoton")             ; b_nLoosePhoton             ->SetAddress(&nLoosePhoton            );
  TBranch *b_nTightPhoton             = events->GetBranch("nTightPhoton")             ; b_nTightPhoton             ->SetAddress(&nTightPhoton            );
  TBranch *b_loosePho1IsTight         = events->GetBranch("loosePho1IsTight")         ; b_loosePho1IsTight         ->SetAddress(&loosePho1IsTight        );
  TBranch *b_loosePho1Pt              = events->GetBranch("loosePho1Pt")              ; b_loosePho1Pt              ->SetAddress(&loosePho1Pt             );
  TBranch *b_loosePho1Eta             = events->GetBranch("loosePho1Eta")             ; b_loosePho1Eta             ->SetAddress(&loosePho1Eta            );
  TBranch *b_loosePho1Phi             = events->GetBranch("loosePho1Phi")             ; b_loosePho1Phi             ->SetAddress(&loosePho1Phi            );
  TBranch *b_diLepMass                = events->GetBranch("diLepMass")                ; b_diLepMass                ->SetAddress(&diLepMass               );
  TBranch *b_nTau                     = events->GetBranch("nTau")                     ; b_nTau                     ->SetAddress(&nTau                    );
  TBranch *b_mT                       = events->GetBranch("mT")                       ; b_mT                       ->SetAddress(&mT                      );
  TBranch *b_scaleUp                  = events->GetBranch("scaleUp")                  ; b_scaleUp                  ->SetAddress(&scaleUp                 );
  TBranch *b_scaleDown                = events->GetBranch("scaleDown")                ; b_scaleDown                ->SetAddress(&scaleDown               );
  TBranch *b_pdfUp                    = events->GetBranch("pdfUp")                    ; b_pdfUp                    ->SetAddress(&pdfUp                   );
  TBranch *b_pdfDown                  = events->GetBranch("pdfDown")                  ; b_pdfDown                  ->SetAddress(&pdfDown                 );
  TBranch *b_normalizedWeight         = events->GetBranch("normalizedWeight")         ; b_normalizedWeight         ->SetAddress(&normalizedWeight        );

  TFile *outputFile = TFile::Open(outputFileName,"RECREATE");
  TTree *plotTree = new TTree("plotTree","Tree for making plots");
  int typeLepSel, theCategory;
  float lepton1Pt, lepton1Eta, lepton1Phi;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1Eta, hbbJet1Phi;
  float hbbJet2Pt, hbbJet2Eta, hbbJet2Phi;
  float vectorBosonPt, hbbDijetPt, hbbDijetMass;
  float bDiscrMin, bDiscrMax;
  float deltaPhiMetLep1, deltaPhiVH;
  float weight;
  
  plotTree->Branch("theCategory"     , &theCategory     );
  plotTree->Branch("nJet"            , &nJet            );
  plotTree->Branch("weight"          , &weight          );
  plotTree->Branch("pfmet"           , &pfmet           );
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
    typeLepSel=-1; // 0: mixed e-mu, 1: all mu, 2: all e
    
    // Vectors to keep around in memory
    TVector2 vectorBosonV2;
    TLorentzVector hbbJet1P4, hbbJet2P4, hbbDijetP4;
    
    // Selection
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
      // Lepton multiplicity
      nBytesRead+=bLoad(b_nLooseLep,ientry);
      if(nLooseLep!=1) continue; //N_al = 0
      if(debug) printf("Passed lepton multiplicity\n");

      // Jet multiplicity
      nBytesRead+=bLoad(b_nJet,ientry);
      if     (nJet<2) continue;
      else if(selection==kWHHeavyFlavorCR && nJet!=2) continue;
      else if(selection==kWH2TopCR && nJet<3) continue;
      else if(selection==kWHSR && nJet>3) continue;
      if(debug) printf("Passed jet multiplicity\n");

      // Lepton ID and isolation
      nBytesRead+=bLoad(b_nTightMuon,ientry);
      if(nTightMuon==1) typeLepSel=1; else {
        nBytesRead+=bLoad(b_nTightElectron,ientry);
        if(nTightElectron==1) typeLepSel=2;
      } if(typeLepSel<0) continue;
      if(debug) printf("Passed lepton ID/iso multiplicity\n");

      // Lepton kinematics
      if     (typeLepSel==1) {
        nBytesRead+=bLoad(b_muonPt,ientry);
        nBytesRead+=bLoad(b_muonEta,ientry);
        nBytesRead+=bLoad(b_muonPhi,ientry);
        lepton1Pt=muonPt[0]; lepton1Eta=muonEta[0]; lepton1Phi=muonPhi[0];
      } else if(typeLepSel==2) {
        nBytesRead+=bLoad(b_electronPt,ientry);
        nBytesRead+=bLoad(b_electronEta,ientry);
        nBytesRead+=bLoad(b_electronPhi,ientry);
        lepton1Pt=electronPt[0]; lepton1Eta=electronEta[0]; lepton1Phi=electronPhi[0];
      } else continue;
      if(debug) printf("Passed lepton kinematics\n");
      
      // Met
      nBytesRead+=bLoad(b_pfmet,ientry);
      if(pfmet<30) continue;
      if(debug) printf("passed MET\n");

      // W reconstruction
      nBytesRead+=bLoad(b_pfmetphi,ientry);
      { 
        TVector2 lepton1V2, pfmetV2;
        lepton1V2.SetMagPhi(lepton1Pt,lepton1Phi);
        pfmetV2.SetMagPhi(pfmet,pfmetphi);
        vectorBosonV2 = lepton1V2+pfmetV2; 
        deltaPhiMetLep1 = pfmetV2.DeltaPhi(lepton1V2);
      } vectorBosonPt=vectorBosonV2.Mod();
      if(vectorBosonPt<50.) continue;
      if(debug) printf("passed W reco\n");

      // Jet kinematics
      nBytesRead+=bLoad(b_hbbjtidx,ientry); // indices of Higgs daughter jets
      nBytesRead+=bLoad(b_jetPt,ientry);
      if(debug) printf("hbb jet1 pt %.2f, jet2 pt %.2f\n", jetPt[hbbjtidx[0]], jetPt[hbbjtidx[1]]);
      if(jetPt[hbbjtidx[0]]<25 || jetPt[hbbjtidx[1]]<25) continue;
      if(usePandaHbb) {
        nBytesRead+=bLoad(b_hbbpt,ientry);
        nBytesRead+=bLoad(b_hbbeta,ientry);
        nBytesRead+=bLoad(b_hbbphi,ientry);
        nBytesRead+=bLoad(b_hbbm,ientry);
        nBytesRead+=bLoad(b_jetEta,ientry);
        nBytesRead+=bLoad(b_jetPhi,ientry);
        hbbDijetPt=hbbpt;
        hbbDijetMass=hbbm;
      } else { 
        nBytesRead+=bLoad(b_jetE,ientry);
        hbbJet1P4.SetPtEtaPhiE(jetPt[hbbjtidx[0]],jetEta[hbbjtidx[0]],jetPhi[hbbjtidx[0]],jetE[hbbjtidx[0]]);
        hbbJet2P4.SetPtEtaPhiE(jetPt[hbbjtidx[1]],jetEta[hbbjtidx[1]],jetPhi[hbbjtidx[1]],jetE[hbbjtidx[1]]);
        hbbDijetP4 = hbbJet1P4+hbbJet2P4;
        hbbDijetPt=hbbDijetP4.Pt(); 
        hbbDijetMass=hbbDijetP4.M();
      }
      hbbJet1Pt=jetPt[hbbjtidx[0]]; hbbJet1Eta=jetEta[hbbjtidx[0]]; hbbJet1Phi=jetPhi[hbbjtidx[0]];
      hbbJet2Pt=jetPt[hbbjtidx[1]]; hbbJet1Eta=jetEta[hbbjtidx[1]]; hbbJet1Phi=jetPhi[hbbjtidx[1]];
      if(debug) printf("hbbDijetPt %.2f, hbbDijetMass %.2f\n", hbbDijetPt, hbbDijetMass); 
      if(hbbDijetPt<50.) continue;
      bool passDijetMass=false;
      if(hbbDijetMass<0 || hbbDijetMass>250.) continue;
      if     (selection==kWHSR           ) passDijetMass = (hbbDijetMass >= 90) && (hbbDijetMass <  150);
      else if(selection==kWHHeavyFlavorCR) passDijetMass = (hbbDijetMass <  90) || (hbbDijetMass >= 150);
      else if(selection==kWHLightFlavorCR) passDijetMass = true;
      else if(selection==kWH2TopCR       ) passDijetMass = true;
      if(!passDijetMass) continue;
      if(debug) printf("passed jet kinematics\n");

      // Jet B-tagging
      nBytesRead+=bLoad(b_jetCSV,ientry);
      bDiscrMax=jetCSV[hbbjtidx[0]];
      bDiscrMin=jetCSV[hbbjtidx[1]];
      bool passBTag=false;
      if     (selection==kWHSR           ) passBTag = (bDiscrMax >= bDiscrTight) && (bDiscrMin >= bDiscrLoose );
      else if(selection==kWHHeavyFlavorCR) passBTag = (bDiscrMax >= bDiscrTight) && (bDiscrMin >= bDiscrLoose );
      else if(selection==kWHLightFlavorCR) passBTag = (bDiscrMax >= bDiscrLoose) && (bDiscrMax <  bDiscrMedium);
      else if(selection==kWH2TopCR       ) passBTag = (bDiscrMax >= bDiscrTight) && (bDiscrMin >= bDiscrLoose );
      if(!passBTag) {
        if(debug) printf("failed btag, bDiscrMax=%.4f, bDiscrMin=%.4f\n",bDiscrMax,bDiscrMin);
        continue;
      }
      if(debug) printf("passed B tagging\n");

      // Other stuff
      nBytesRead+=bLoad(b_hbbpt,ientry);
      nBytesRead+=bLoad(b_hbbphi,ientry);
      { 
        TVector2 hbbDijetV2; hbbDijetV2.SetMagPhi(hbbpt, hbbphi);
        deltaPhiVH = vectorBosonV2.DeltaPhi(hbbDijetV2);
      }
    }

    // Weighting
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
      nBytesRead+=bLoad(b_normalizedWeight,ientry);
	  nBytesRead+=bLoad(b_sf_pu,ientry);
	  if(sample==kWjets || sample==kZjets) { nBytesRead+=bLoad(b_sf_qcdV,ientry); nBytesRead+=bLoad(b_sf_ewkV,ientry); } 
      else { sf_qcdV=1; sf_ewkV=1; }
	  //nBytesRead+=bLoad(sf_btag1,ientry);
 	  if(sample==kTT) nBytesRead+=bLoad(b_sf_tt,ientry); else sf_tt=1;
      weight = normalizedWeight * sf_pu * sf_ewkV * sf_qcdV * sf_tt;
      if (typeLepSel==1) {
        nBytesRead+=bLoad(b_muonSfReco,ientry);
        nBytesRead+=bLoad(b_muonSfTight,ientry);
        // need muon trigger efficiency
        weight *= muonSfReco[0] * muonSfTight[0];
      } else if(typeLepSel==2) {
        nBytesRead+=bLoad(b_electronSfReco,ientry);
        nBytesRead+=bLoad(b_electronSfTight,ientry);
        nBytesRead+=bLoad(b_sf_eleTrig,ientry);
        weight *= sf_eleTrig * electronSfReco[0] * electronSfTight[0];
      }
      if(sample==kVZ) {
        nBytesRead+=bLoad(b_sf_wz,ientry);
        nBytesRead+=bLoad(b_sf_zz,ientry);
        //nBytesRead+=bLoad(b_sf_zzUnc,ientry);
        weight *= sf_wz * sf_zz;
      } else if(sample==kVH) {
        nBytesRead+=bLoad(b_sf_zh,ientry);
        //nBytesRead+=bLoad(b_sf_zhUp,ientry);
        //nBytesRead+=bLoad(b_sf_zhDown,ientry);
        weight *= sf_zh;
      }
    }

    // Category Assignment
    int nBQuark=0;
    if(sample==kVZ || sample==kWjets || kZjets) {
      nBytesRead+=bLoad(b_jetGenFlavor,ientry);
      for(int iJ=0; iJ<nJet; iJ++) if(jetGenFlavor[iJ]==5) nBQuark++;
    }
    switch(sample) {
      case kData:
        theCategory=kPlotData; break;
      case kQCD:
        theCategory=kPlotQCD;  break;
      case kVZ:
        if(nBQuark<2) theCategory=kPlotVVLF;
        else          theCategory=kPlotVVLF;
        break;
      case kWW:
        theCategory=kPlotVVLF; break;
      case kTT:
        theCategory=kPlotTT; break;
      case kTop:
        theCategory=kPlotTop; break;
      case kWjets:
        if(nBQuark>=2) theCategory=kPlotWbb;
        else if(nBQuark==1) theCategory=kPlotWb;
        else theCategory=kPlotWLF; // light flavor
        break;
      case kZjets:
        if(nBQuark>=2) theCategory=kPlotZbb;
        else if(nBQuark==1) theCategory=kPlotZb;
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




}
