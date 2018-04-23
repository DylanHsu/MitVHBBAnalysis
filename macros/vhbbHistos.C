//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PandaAnalysis/Utilities/src/CSVHelper.cc"
#include <Compression.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TTree.h>
#include <TMath.h>
#include "TMVA/Reader.h"
#include "parseTmvaXmlFile.C"

#include <cassert>
#include <iostream>
#include <fstream>

#include "vhbbPlot.h"

const vector<TString> wptCorrFilenames = {
  "MitVHBBAnalysis/wptCorrections_InclusiveResolved_Linear.root",
  "MitVHBBAnalysis/wptCorrections_Resolved_Linear.root",
  "MitVHBBAnalysis/wptCorrections_Boosted_Linear.root",
};
const bool useIgnorantVHFSFs=true;

using namespace vhbbPlot;
void vhbbHistos(
  TString plotTreeFileName,
  TString dataCardDir,
  vhbbPlot::selectionType selection,
  int theLepSel=-1, // -1: any; 0, e-mu; 1, mu only; 2, ele only
  bool isBlinded=true,
  int MVAVarType=3, //1=Higgs pT; 2=Multiclass BDT Output; 3=Single-Class BDT;
  bool lowMass=true, // only for WH HF CR
  int wptCorrType=-1,
  TString bdtWeights=""
) {
  if(dataCardDir!="") system(Form("mkdir -p MitVHBBAnalysis/datacards/%s",dataCardDir.Data()));
  bool useWptCorr=(wptCorrType>=0); // && wptCorrTypewptCorrFilenames.size());
  TString wptCorrFilename;
  //if(useWptCorr) wptCorrFilename = wptCorrFilenames[wptCorrType];
  if(useWptCorr) wptCorrFilename = "MitVHBBAnalysis/wptCorrections_InclusiveResolved_Linear.root";
    
  char shapeType[64]; TString bdtOutputName="BDT output";
  vector<float> MVAbins; TString MVAVarName;
  if(MVAVarType==1) {
    MVAVarName="Higgs p_{T} classifier [GeV]";
    if(selection>=kWHLightFlavorCR && selection<=kWHPresel) {
      MVAbins={100,120,140,160,180,200,250,300,350};
      MVAVarName="H(bb) pT";
      sprintf(shapeType,"ptShape");
    } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
      MVAbins={250,300,350,400,450,500,550,600};
      MVAVarName="H(bb) pT";
      sprintf(shapeType,"ptShape");
    } else throw std::runtime_error("bad selection argument");
  } else if(MVAVarType==2) {
    if(selection==kWHLightFlavorCR) {
      MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667};
      //MVAbins={-1,-.8,-.6,-.4,-.2,0,.2,.4,.6}; 
      MVAVarName="Subleading H(bb) CMVA";
      sprintf(shapeType,"lesserCMVAShape");
    } else if(selection==kWHHeavyFlavorCR || selection==kWH2TopCR) {
      MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000, 0.7333, 0.8667, 1.0000};
      MVAVarName="Subleading H(bb) CMVA";
      sprintf(shapeType,"lesserCMVAShape");
    } else if(selection>=kWHLightFlavorFJCR && selection<=kWHTT1bFJCR) {
      MVAbins={50,70,90,110,130,150,170,190,210,230,250};
      MVAVarName="Fatjet soft drop mass [GeV]";
      sprintf(shapeType,"softDropMassShape");
    } else if(selection==kWHSR || selection==kWHFJSR) {
      MVAbins={0,0.20,0.40,0.60,0.70,0.80,0.85,0.9,0.95,1.00};
      MVAVarName="Multiclass BDT Output";
      sprintf(shapeType,"multiClassBDTShape");
    }
    else throw std::runtime_error("bad selection argument");
    if(selection==kWHLightFlavorCR || selection==kWHLightFlavorFJCR) bdtOutputName="BDT output for W+LF likelihood";
    else if(selection==kWHHeavyFlavorCR || selection==kWHHeavyFlavorFJCR) bdtOutputName="BDT output for W+b(b) likelihood";
    else if(selection==kWH2TopCR || selection==kWHTT2bFJCR || selection==kWHTT1bFJCR) bdtOutputName="BDT output for top/t#bar{t} likelihood";
    else bdtOutputName="BDT output for WH(bb) likelihood";
  } else if(MVAVarType==3) {
    if(selection==kWHLightFlavorCR) {
      //MVAbins={-1,-.8,-.6,-.4,-.2,0,.2,.4,.6}; 
      MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000};
      MVAVarName="Subleading H(bb) CMVA";
      sprintf(shapeType,"lesserCMVAShape");
    } else if(selection==kWHHeavyFlavorCR || selection==kWH2TopCR) {
      //MVAbins={-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1}; 
      MVAbins={-1.0000, -0.8667, -0.7333, -0.6000, -0.4667, -0.3333, -0.2000, -0.0667, 0.0667, 0.2000, 0.3333, 0.4667, 0.6000, 0.7333, 0.8667, 1.0000};
      MVAVarName="Subleading H(bb) CMVA";
      sprintf(shapeType,"lesserCMVAShape");
    } else if(selection==kWHSR) {
      MVAbins={-1,-0.5,0, 0.20,0.40,0.60,0.70,0.80,0.85,0.9,0.94,1.00};
      //MVAbins={-1,-0.5,0,0.30,0.5,0.70,0.85,1.00};
      MVAVarName="BDT Output";
      sprintf(shapeType,"singleClassBDTShape"); 
    } else if(selection>=kWHLightFlavorFJCR && selection<kWHFJSR) {
      if(selection==kWHHeavyFlavorFJCR) MVAbins={40,50,60,70,80};
      else                              MVAbins={40,50,60,70,80,100,150,200};
      MVAVarName="Fatjet soft drop mass [GeV]";
      sprintf(shapeType,"softDropMassShape");
    } else if(selection==kWHFJSR) {
      //MVAbins={-1,-.5,0,0.40,0.75,1.00};
      MVAbins={-1,0,.2,.4,0.6,0.7,1};
      MVAVarName="BDT Output";
      sprintf(shapeType,"singleClassBDTShape"); 
    }
    else throw std::runtime_error("bad selection argument");
  } else if(MVAVarType==4) { // cut and count sanity check
    if(selection==kWHLightFlavorCR) {
      MVAbins={-1,1};
      MVAVarName="Subleading H(bb) CMVA";
      sprintf(shapeType,"lesserCMVAShape");
    } else if(selection==kWHHeavyFlavorCR || selection==kWH2TopCR) {
      MVAbins={-1,1};
      MVAVarName="Subleading H(bb) CMVA";
      sprintf(shapeType,"lesserCMVAShape");
    } else if(selection==kWHSR) {
      MVAbins={-1,1};
      MVAVarName="BDT Output";
      sprintf(shapeType,"singleClassBDTShape"); 
    } else if(selection>=kWHLightFlavorFJCR && selection<kWHFJSR) {
      if(selection==kWHHeavyFlavorFJCR) MVAbins={40,80};
      else                              MVAbins={40,200};
      MVAVarName="Fatjet soft drop mass [GeV]";
      sprintf(shapeType,"softDropMassShape");
    } else if(selection==kWHFJSR) {
      MVAbins={-1,1};
      MVAVarName="BDT Output";
      sprintf(shapeType,"singleClassBDTShape"); 
    }
    else throw std::runtime_error("bad selection argument");
  } else throw std::runtime_error("bad MVAVarType");
  
  // Load Files
  gSystem->Load("libCondFormatsJetMETObjects.so");
  system("mkdir -p MitVHBBAnalysis/datacards");
  system("mkdir -p MitVHBBAnalysis/plots");
  TFile *plotTreeFile = TFile::Open(plotTreeFileName, "READ");
  assert(plotTreeFile && plotTreeFile->IsOpen());
  //string theJecUncertainties = "PandaAnalysis/data/jec/23Sep2016V4/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt";
  //JetCorrectionUncertainty *jcuTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters(theJecUncertainties, "Total")));
  CSVHelper *cmvaReweighter = new CSVHelper("PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root"   , "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root"   , 5);
  TFile *wptCorrFile=0; TF1 *theWptCorr=0; TH1F *theWptCorrHist=0;
  if(useWptCorr && wptCorrType==0) {
    wptCorrFile = TFile::Open(wptCorrFilename, "READ");
    assert(wptCorrFile && wptCorrFile->IsOpen());
    theWptCorr = (TF1*)wptCorrFile->Get("wptCorr_nominal");
    assert(theWptCorr);
    theWptCorrHist = (TH1F*)wptCorrFile->Get("histWithStatUnc_wptCorr_nominal");
    assert(theWptCorrHist);
  }
  // Initialize the tree
  TTree *plotTree = (TTree*)plotTreeFile->Get("plotTree"); assert(plotTree);
  int runNumber = -1;
  int lumiNumber = -1;
  ULong64_t eventNumber = -1;
  unsigned char typeLepSel, theCategory;
  unsigned selectionBits, selectionBits_jesUp, selectionBits_jesDown, nMinusOneBits;
  int nJet, nSoft2, nSoft5, nSoft10, nIsojet, nFatjet;
  unsigned char isojetNBtags;
  float sumEtSoft1;
  float pfmet, pfmetphi, pfmetsig, pfmetUp, pfmetDown;
  float lepton1Pt, lepton1Eta, lepton1Phi, lepton1RelIso;
  int lepton1Flav, lepton1Charge;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1Eta, hbbJet1Phi, hbbJet1PtUp, hbbJet1PtDown;
  float hbbJet2Pt, hbbJet2Eta, hbbJet2Phi, hbbJet2PtUp, hbbJet2PtDown;
  float topWBosonPt, topWBosonEta, topWBosonPhi, topWBosonCosThetaCS;
  float topWBosonPt_jesUp, topWBosonPt_jesDown, topWBosonPhi_jesUp, topWBosonPhi_jesDown;
  float hbbDijetPt, hbbDijetPtUp, hbbDijetPtDown;
  float hbbDijetMass, hbbDijetMassUp, hbbDijetMassDown;
  float hbbCosThetaJJ, hbbCosThetaCSJ1;
  float bDiscrMin, bDiscrMax;
  float deltaPhiLep1Met; 
  float deltaPhiVH, deltaPhiVH_jesUp, deltaPhiVH_jesDown; 
  float topMassLep1Met, topMassLep1Met_jesUp, topMassLep1Met_jesDown;
  float weight;
  float weight_VHCorrUp, weight_VHCorrDown;
  float weight_pdfUp, weight_pdfDown;
  float weight_QCDr1f2, weight_QCDr1f5, weight_QCDr2f1, weight_QCDr2f2, weight_QCDr5f1, weight_QCDr5f5, weight_lepSFUp;
  float weight_NLOQCDfUp,weight_NLOQCDfDown,weight_NLOQCDrUp,weight_NLOQCDrDown;
  float weight_cmvaJESUp     [5][3] , weight_cmvaJESDown     [5][3] , 
        weight_cmvaLFUp      [5][3] , weight_cmvaLFDown      [5][3] , 
        weight_cmvaHFUp      [5][3] , weight_cmvaHFDown      [5][3] , 
        weight_cmvaHFStats1Up[5][3] , weight_cmvaHFStats1Down[5][3] , 
        weight_cmvaHFStats2Up[5][3] , weight_cmvaHFStats2Down[5][3] , 
        weight_cmvaLFStats1Up[5][3] , weight_cmvaLFStats1Down[5][3] , 
        weight_cmvaLFStats2Up[5][3] , weight_cmvaLFStats2Down[5][3] , 
        weight_cmvaCErr1Up   [5][3] , weight_cmvaCErr1Down   [5][3] , 
        weight_cmvaCErr2Up   [5][3] , weight_cmvaCErr2Down   [5][3] ; 
  float mT,mT_jesUp,mT_jesDown;
  float fj1Tau32;
  float fj1Tau21;
  float fj1Tau32SD;
  float fj1Tau21SD;
  float fj1MSD, fj1MSD_corr,fj1MSD_corr_jesUp,fj1MSD_corr_jesDown;
  float fj1MSDScaleUp;
  float fj1MSDScaleDown;
  float fj1MSDSmeared;
  float fj1MSDSmearedUp;
  float fj1MSDSmearedDown;
  float fj1Pt;
  float fj1PtScaleUp;
  float fj1PtScaleDown;
  float fj1PtSmeared;
  float fj1PtSmearedUp;
  float fj1PtSmearedDown;
  float fj1Phi;
  float fj1Eta, adjustedFatjetEta;
  float fj1M;
  float fj1MaxCSV;
  float fj1MinCSV;
  float fj1DoubleCSV;
  float fj1HTTMass;
  float fj1HTTFRec;
  float fj1SDEFrac100;
  int fj1HighestPtGen;
  std::map<TString, float> fj1ECFN, mvaPsi;
  fj1ECFN["1_1_05"]=-3393; 
  fj1ECFN["2_1_05"]=-3393; 
  fj1ECFN["3_1_05"]=-3393; 
  fj1ECFN["1_2_05"]=-3393; 
  fj1ECFN["2_2_05"]=-3393; 
  fj1ECFN["3_2_05"]=-3393; 
  fj1ECFN["1_3_05"]=-3393; 
  fj1ECFN["2_3_05"]=-3393; 
  fj1ECFN["3_3_05"]=-3393; 
  fj1ECFN["1_4_05"]=-3393; 
  fj1ECFN["2_4_05"]=-3393; 
  fj1ECFN["3_4_05"]=-3393; 
  fj1ECFN["1_1_10"]=-3393; 
  fj1ECFN["2_1_10"]=-3393; 
  fj1ECFN["3_1_10"]=-3393; 
  fj1ECFN["1_2_10"]=-3393; 
  fj1ECFN["2_2_10"]=-3393; 
  fj1ECFN["3_2_10"]=-3393; 
  fj1ECFN["1_3_10"]=-3393; 
  fj1ECFN["2_3_10"]=-3393; 
  fj1ECFN["3_3_10"]=-3393; 
  fj1ECFN["1_4_10"]=-3393; 
  fj1ECFN["2_4_10"]=-3393; 
  fj1ECFN["3_4_10"]=-3393; 
  fj1ECFN["1_1_20"]=-3393; 
  fj1ECFN["2_1_20"]=-3393; 
  fj1ECFN["3_1_20"]=-3393; 
  fj1ECFN["1_2_20"]=-3393; 
  fj1ECFN["2_2_20"]=-3393; 
  fj1ECFN["3_2_20"]=-3393; 
  fj1ECFN["1_3_20"]=-3393; 
  fj1ECFN["2_3_20"]=-3393; 
  fj1ECFN["3_3_20"]=-3393; 
  fj1ECFN["1_4_20"]=-3393; 
  fj1ECFN["2_4_20"]=-3393; 
  fj1ECFN["3_4_20"]=-3393; 
  fj1ECFN["1_1_40"]=-3393; 
  fj1ECFN["2_1_40"]=-3393; 
  fj1ECFN["3_1_40"]=-3393; 
  fj1ECFN["1_2_40"]=-3393; 
  fj1ECFN["2_2_40"]=-3393; 
  fj1ECFN["3_2_40"]=-3393; 
  fj1ECFN["1_3_40"]=-3393; 
  fj1ECFN["2_3_40"]=-3393; 
  fj1ECFN["3_3_40"]=-3393; 
  fj1ECFN["1_4_40"]=-3393; 
  fj1ECFN["2_4_40"]=-3393; 
  fj1ECFN["3_4_40"]=-3393; 
  mvaPsi["021004010502"] = -3393; 
  mvaPsi["012004010502"] = -3393; 
  mvaPsi["021003011002"] = -3393; 
  mvaPsi["022004011002"] = -3393; 
  mvaPsi["020503010502"] = -3393; 
  mvaPsi["022003012002"] = -3393; 
  mvaPsi["012004020503"] = -3393; 
  mvaPsi["021004020503"] = -3393; 
  mvaPsi["032003012002"] = -3393; 
  mvaPsi["012003011002"] = -3393; 
  mvaPsi["011003010502"] = -3393; 
  mvaPsi["031003011002"] = -3393; 
  mvaPsi["031003012002"] = -3393; 
  mvaPsi["030503011002"] = -3393; 
  mvaPsi["022004012002"] = -3393; 
  mvaPsi["021004011002"] = -3393; 
  mvaPsi["012004011002"] = -3393; 
  mvaPsi["022004021003"] = -3393; 
  mvaPsi["012003010502"] = -3393; 
  mvaPsi["022003011002"] = -3393; 
  mvaPsi["022004031003"] = -3393; 
  mvaPsi["012004030503"] = -3393; 
  mvaPsi["021004030503"] = -3393; 
  mvaPsi["020504010502"] = -3393; 
  mvaPsi["022004030503"] = -3393; 
  mvaPsi["011004010502"] = -3393; 
  mvaPsi["030503010502"] = -3393; 
  mvaPsi["021003010502"] = -3393; 
  mvaPsi["012004021004"] = -3393; 
  mvaPsi["021004012004"] = -3393; 
  mvaPsi["011003020503"] = -3393; 
  mvaPsi["020503011003"] = -3393; 
  mvaPsi["012003030503"] = -3393; 
  mvaPsi["030503012003"] = -3393; 
  mvaPsi["010503010502"] = -3393; 
  mvaPsi["020503011002"] = -3393; 
  mvaPsi["011003011002"] = -3393; 
  mvaPsi["020504011002"] = -3393; 
  mvaPsi["012004021003"] = -3393; 
  mvaPsi["012004012002"] = -3393; 
  mvaPsi["020504020503"] = -3393; 
  mvaPsi["011004011002"] = -3393; 
  mvaPsi["021004012002"] = -3393; 
  mvaPsi["021003012002"] = -3393; 
  mvaPsi["011004020503"] = -3393; 
  mvaPsi["010504010502"] = -3393; 
  mvaPsi["021004021003"] = -3393; 
  mvaPsi["012003012002"] = -3393; 
  mvaPsi["022004022003"] = -3393; 
  // Sid's "clean" monotop set
  mvaPsi["012002011002"] = -3393;
  mvaPsi["014003022003"] = -3393;
  mvaPsi["031003014003"] = -3393;
  mvaPsi["031003022003"] = -3393;
  mvaPsi["032003034003"] = -3393;
  mvaPsi["012004011003"] = -3393;
  mvaPsi["014004012003"] = -3393;
  mvaPsi["020504010503"] = -3393;
  mvaPsi["021004011003"] = -3393;
  mvaPsi["021004020503"] = -3393;
  mvaPsi["022004012003"] = -3393;


  plotTree->SetBranchAddress("runNumber"       , &runNumber           );
  plotTree->SetBranchAddress("lumiNumber"      , &lumiNumber          );
  plotTree->SetBranchAddress("eventNumber"     , &eventNumber         );
  plotTree->SetBranchAddress("selectionBits"   , &selectionBits   );
  plotTree->SetBranchAddress("selectionBits_jesUp"     , &selectionBits_jesUp   );
  plotTree->SetBranchAddress("selectionBits_jesDown"   , &selectionBits_jesDown );
  plotTree->SetBranchAddress("nMinusOneBits"   , &nMinusOneBits   );
  plotTree->SetBranchAddress("theCategory"     , &theCategory     );
  plotTree->SetBranchAddress("pfmet"           , &pfmet           );
  plotTree->SetBranchAddress("pfmetUp"         , &pfmetUp         );
  plotTree->SetBranchAddress("pfmetDown"       , &pfmetDown       );
  plotTree->SetBranchAddress("pfmetsig"        , &pfmetsig        );
  plotTree->SetBranchAddress("pfmetphi"        , &pfmetphi        );
  plotTree->SetBranchAddress("nFatjet"         , &nFatjet         );
  plotTree->SetBranchAddress("weight"          , &weight          );
  plotTree->SetBranchAddress("weight_VHCorrUp"    , &weight_VHCorrUp     );
  plotTree->SetBranchAddress("weight_VHCorrDown"  , &weight_VHCorrDown   );
  //plotTree->SetBranchAddress("weight_pdfUp"       , &weight_pdfUp        );
  //plotTree->SetBranchAddress("weight_pdfDown"     , &weight_pdfDown      );
  plotTree->SetBranchAddress("weight_QCDr1f2"     , &weight_QCDr1f2      );
  plotTree->SetBranchAddress("weight_QCDr1f5"     , &weight_QCDr1f5      );
  plotTree->SetBranchAddress("weight_QCDr2f1"     , &weight_QCDr2f1      );
  plotTree->SetBranchAddress("weight_QCDr2f2"     , &weight_QCDr2f2      );
  plotTree->SetBranchAddress("weight_QCDr5f1"     , &weight_QCDr5f1      );
  plotTree->SetBranchAddress("weight_QCDr5f5"     , &weight_QCDr5f5      );
  plotTree->SetBranchAddress("weight_NLOQCDfUp"   , &weight_NLOQCDfUp    );
  plotTree->SetBranchAddress("weight_NLOQCDfDown" , &weight_NLOQCDfDown  );
  plotTree->SetBranchAddress("weight_NLOQCDrUp"   , &weight_NLOQCDrUp    );
  plotTree->SetBranchAddress("weight_NLOQCDrDown" , &weight_NLOQCDrDown  );
  plotTree->SetBranchAddress("weight_cmvaJESUp"       , weight_cmvaJESUp       );
  plotTree->SetBranchAddress("weight_cmvaLFUp"        , weight_cmvaLFUp        );
  plotTree->SetBranchAddress("weight_cmvaHFUp"        , weight_cmvaHFUp        );
  plotTree->SetBranchAddress("weight_cmvaHFStats1Up"  , weight_cmvaHFStats1Up  );
  plotTree->SetBranchAddress("weight_cmvaHFStats2Up"  , weight_cmvaHFStats2Up  );
  plotTree->SetBranchAddress("weight_cmvaLFStats1Up"  , weight_cmvaLFStats1Up  );
  plotTree->SetBranchAddress("weight_cmvaLFStats2Up"  , weight_cmvaLFStats2Up  );
  plotTree->SetBranchAddress("weight_cmvaCErr1Up"     , weight_cmvaCErr1Up     );
  plotTree->SetBranchAddress("weight_cmvaCErr2Up"     , weight_cmvaCErr2Up     );
  plotTree->SetBranchAddress("weight_cmvaJESDown"     , weight_cmvaJESDown     );
  plotTree->SetBranchAddress("weight_cmvaLFDown"      , weight_cmvaLFDown      );
  plotTree->SetBranchAddress("weight_cmvaHFDown"      , weight_cmvaHFDown      );
  plotTree->SetBranchAddress("weight_cmvaHFStats1Down", weight_cmvaHFStats1Down);
  plotTree->SetBranchAddress("weight_cmvaHFStats2Down", weight_cmvaHFStats2Down);
  plotTree->SetBranchAddress("weight_cmvaLFStats1Down", weight_cmvaLFStats1Down);
  plotTree->SetBranchAddress("weight_cmvaLFStats2Down", weight_cmvaLFStats2Down);
  plotTree->SetBranchAddress("weight_cmvaCErr1Down"   , weight_cmvaCErr1Down   );
  plotTree->SetBranchAddress("weight_cmvaCErr2Down"   , weight_cmvaCErr2Down   );
  plotTree->SetBranchAddress("weight_lepSFUp"     , &weight_lepSFUp      );
  if(selection>=kWHLightFlavorCR && selection<=kWHPresel) {
    plotTree->SetBranchAddress("nJet"                  , &nJet                  );
    plotTree->SetBranchAddress("typeLepSel"            , &typeLepSel            );
    plotTree->SetBranchAddress("hbbDijetPt"            , &hbbDijetPt            );
    plotTree->SetBranchAddress("hbbDijetPtUp"          , &hbbDijetPtUp          );
    plotTree->SetBranchAddress("hbbDijetPtDown"        , &hbbDijetPtDown        );
    plotTree->SetBranchAddress("hbbDijetMass"          , &hbbDijetMass          );
    plotTree->SetBranchAddress("hbbDijetMassUp"        , &hbbDijetMassUp        );
    plotTree->SetBranchAddress("hbbDijetMassDown"      , &hbbDijetMassDown      );
    plotTree->SetBranchAddress("hbbCosThetaJJ"         , &hbbCosThetaJJ         );
    plotTree->SetBranchAddress("hbbCosThetaCSJ1"       , &hbbCosThetaCSJ1       );
    plotTree->SetBranchAddress("topMassLep1Met"        , &topMassLep1Met        );
    plotTree->SetBranchAddress("topMassLep1Met_jesUp"  , &topMassLep1Met_jesUp  );
    plotTree->SetBranchAddress("topMassLep1Met_jesDown", &topMassLep1Met_jesDown);
    plotTree->SetBranchAddress("bDiscrMin"             , &bDiscrMin             );
    plotTree->SetBranchAddress("bDiscrMax"             , &bDiscrMax             );
    plotTree->SetBranchAddress("hbbJet1Pt"             , &hbbJet1Pt             );
    plotTree->SetBranchAddress("hbbJet1PtUp"           , &hbbJet1PtUp           );
    plotTree->SetBranchAddress("hbbJet1PtDown"         , &hbbJet1PtDown         );
    plotTree->SetBranchAddress("hbbJet1Eta"            , &hbbJet1Eta            );
    plotTree->SetBranchAddress("hbbJet1Phi"            , &hbbJet1Phi            );
    plotTree->SetBranchAddress("hbbJet2Pt"             , &hbbJet2Pt             );
    plotTree->SetBranchAddress("hbbJet2PtUp"           , &hbbJet2PtUp           );
    plotTree->SetBranchAddress("hbbJet2PtDown"         , &hbbJet2PtDown         );
    plotTree->SetBranchAddress("hbbJet2Eta"            , &hbbJet2Eta            );
    plotTree->SetBranchAddress("hbbJet2Phi"            , &hbbJet2Phi            );
    plotTree->SetBranchAddress("lepton1Pt"             , &lepton1Pt             );
    plotTree->SetBranchAddress("lepton1Eta"            , &lepton1Eta            );
    plotTree->SetBranchAddress("lepton1Phi"            , &lepton1Phi            );
    plotTree->SetBranchAddress("lepton1Flav"           , &lepton1Flav           );
    plotTree->SetBranchAddress("lepton1Charge"         , &lepton1Charge         );
    plotTree->SetBranchAddress("topWBosonCosThetaCS"   , &topWBosonCosThetaCS   );
    plotTree->SetBranchAddress("topWBosonPt"           , &topWBosonPt           );
    plotTree->SetBranchAddress("topWBosonPt_jesUp"     , &topWBosonPt_jesUp     );
    plotTree->SetBranchAddress("topWBosonPt_jesDown"   , &topWBosonPt_jesDown   );
    plotTree->SetBranchAddress("topWBosonEta"          , &topWBosonEta          );
    plotTree->SetBranchAddress("topWBosonPhi"          , &topWBosonPhi          );
    plotTree->SetBranchAddress("topWBosonPhi_jesUp"    , &topWBosonPhi_jesUp    );
    plotTree->SetBranchAddress("topWBosonPhi_jesDown"  , &topWBosonPhi_jesDown  );
    plotTree->SetBranchAddress("mT"                    , &mT                    );
    plotTree->SetBranchAddress("mT_jesUp"              , &mT_jesUp              );
    plotTree->SetBranchAddress("mT_jesDown"            , &mT_jesDown            );
    plotTree->SetBranchAddress("deltaPhiLep1Met"       , &deltaPhiLep1Met       );
    plotTree->SetBranchAddress("deltaPhiVH"            , &deltaPhiVH            );
    plotTree->SetBranchAddress("deltaPhiVH_jesUp"      , &deltaPhiVH_jesUp      );
    plotTree->SetBranchAddress("deltaPhiVH_jesDown"    , &deltaPhiVH_jesDown    );
    plotTree->SetBranchAddress("nSoft2"                , &nSoft2                );
    plotTree->SetBranchAddress("nSoft5"                , &nSoft5                );
    plotTree->SetBranchAddress("nSoft10"               , &nSoft10               );
    plotTree->SetBranchAddress("sumEtSoft1"            , &sumEtSoft1            );
    /*
    plotTree->SetBranchAddress("nIsojet"               , &nIsojet               );
    plotTree->SetBranchAddress("fj1Tau32"              , &fj1Tau32              );
    plotTree->SetBranchAddress("fj1Tau21"              , &fj1Tau21              );
    plotTree->SetBranchAddress("fj1Tau32SD"            , &fj1Tau32SD            ); 
    plotTree->SetBranchAddress("fj1Tau21SD"            , &fj1Tau21SD            ); 
    plotTree->SetBranchAddress("fj1MSD"                , &fj1MSD                );
    plotTree->SetBranchAddress("fj1MSDScaleUp"         , &fj1MSDScaleUp         );    
    plotTree->SetBranchAddress("fj1MSDScaleDown"       , &fj1MSDScaleDown       );      
    plotTree->SetBranchAddress("fj1MSDSmeared"         , &fj1MSDSmeared         );    
    plotTree->SetBranchAddress("fj1MSDSmearedUp"       , &fj1MSDSmearedUp       );      
    plotTree->SetBranchAddress("fj1MSDSmearedDown"     , &fj1MSDSmearedDown     );        
    plotTree->SetBranchAddress("fj1Pt"                 , &fj1Pt                 );
    plotTree->SetBranchAddress("fj1PtScaleUp"          , &fj1PtScaleUp          );   
    plotTree->SetBranchAddress("fj1PtScaleDown"        , &fj1PtScaleDown        );     
    plotTree->SetBranchAddress("fj1PtSmeared"          , &fj1PtSmeared          );   
    plotTree->SetBranchAddress("fj1PtSmearedUp"        , &fj1PtSmearedUp        );     
    plotTree->SetBranchAddress("fj1PtSmearedDown"      , &fj1PtSmearedDown      );       
    plotTree->SetBranchAddress("fj1Phi"                , &fj1Phi                );
    plotTree->SetBranchAddress("fj1Eta"                , &fj1Eta                );
    plotTree->SetBranchAddress("fj1M"                  , &fj1M                  );
    plotTree->SetBranchAddress("fj1MaxCSV"             , &fj1MaxCSV             );
    plotTree->SetBranchAddress("fj1MinCSV"             , &fj1MinCSV             );
    plotTree->SetBranchAddress("fj1DoubleCSV"          , &fj1DoubleCSV          );   
    plotTree->SetBranchAddress("fj1HTTMass"            , &fj1HTTMass            ); 
    plotTree->SetBranchAddress("fj1HTTFRec"            , &fj1HTTFRec            ); 
    plotTree->SetBranchAddress("fj1SDEFrac100"         , &fj1SDEFrac100         );    
    plotTree->SetBranchAddress("fj1ECFN_1_1_05", &fj1ECFN["1_1_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_05", &fj1ECFN["2_1_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_05", &fj1ECFN["3_1_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_05", &fj1ECFN["1_2_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_05", &fj1ECFN["2_2_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_05", &fj1ECFN["3_2_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_05", &fj1ECFN["1_3_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_05", &fj1ECFN["2_3_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_05", &fj1ECFN["3_3_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_05", &fj1ECFN["1_4_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_05", &fj1ECFN["2_4_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_05", &fj1ECFN["3_4_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_1_10", &fj1ECFN["1_1_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_10", &fj1ECFN["2_1_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_10", &fj1ECFN["3_1_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_10", &fj1ECFN["1_2_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_10", &fj1ECFN["2_2_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_10", &fj1ECFN["3_2_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_10", &fj1ECFN["1_3_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_10", &fj1ECFN["2_3_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_10", &fj1ECFN["3_3_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_10", &fj1ECFN["1_4_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_10", &fj1ECFN["2_4_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_10", &fj1ECFN["3_4_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_1_20", &fj1ECFN["1_1_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_20", &fj1ECFN["2_1_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_20", &fj1ECFN["3_1_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_20", &fj1ECFN["1_2_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_20", &fj1ECFN["2_2_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_20", &fj1ECFN["3_2_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_20", &fj1ECFN["1_3_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_20", &fj1ECFN["2_3_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_20", &fj1ECFN["3_3_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_20", &fj1ECFN["1_4_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_20", &fj1ECFN["2_4_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_20", &fj1ECFN["3_4_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_1_40", &fj1ECFN["1_1_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_40", &fj1ECFN["2_1_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_40", &fj1ECFN["3_1_40"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_40", &fj1ECFN["1_2_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_40", &fj1ECFN["2_2_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_40", &fj1ECFN["3_2_40"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_40", &fj1ECFN["1_3_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_40", &fj1ECFN["2_3_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_40", &fj1ECFN["3_3_40"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_40", &fj1ECFN["1_4_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_40", &fj1ECFN["2_4_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_40", &fj1ECFN["3_4_40"]);
    */
  } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
    plotTree->SetBranchAddress("topWBosonPt"           , &topWBosonPt           );
    plotTree->SetBranchAddress("topWBosonPt_jesUp"     , &topWBosonPt_jesUp     );
    plotTree->SetBranchAddress("topWBosonPt_jesDown"   , &topWBosonPt_jesDown   );
    plotTree->SetBranchAddress("topWBosonPhi"          , &topWBosonPhi          );
    plotTree->SetBranchAddress("topWBosonPhi_jesUp"    , &topWBosonPhi_jesUp    );
    plotTree->SetBranchAddress("topWBosonPhi_jesDown"  , &topWBosonPhi_jesDown  );
    plotTree->SetBranchAddress("deltaPhiLep1Met"       , &deltaPhiLep1Met       );
    plotTree->SetBranchAddress("deltaPhiVH"            , &deltaPhiVH            );
    plotTree->SetBranchAddress("deltaPhiVH_jesUp"      , &deltaPhiVH_jesUp      );
    plotTree->SetBranchAddress("deltaPhiVH_jesDown"    , &deltaPhiVH_jesDown    );
    plotTree->SetBranchAddress("typeLepSel"            , &typeLepSel            );
    plotTree->SetBranchAddress("lepton1Pt"             , &lepton1Pt             );
    plotTree->SetBranchAddress("lepton1Eta"            , &lepton1Eta            );
    plotTree->SetBranchAddress("lepton1Phi"            , &lepton1Phi            );
    plotTree->SetBranchAddress("lepton1RelIso"         , &lepton1RelIso         );
    plotTree->SetBranchAddress("lepton1Flav"           , &lepton1Flav           );
    plotTree->SetBranchAddress("lepton1Charge"         , &lepton1Charge         );
    plotTree->SetBranchAddress("nIsojet"               , &nIsojet               );
    plotTree->SetBranchAddress("isojetNBtags"          , &isojetNBtags          );
    plotTree->SetBranchAddress("mT"                    , &mT                    );
    plotTree->SetBranchAddress("mT_jesUp"              , &mT_jesUp              );
    plotTree->SetBranchAddress("mT_jesDown"            , &mT_jesDown            );
    plotTree->SetBranchAddress("fj1Tau32"              , &fj1Tau32              );
    plotTree->SetBranchAddress("fj1Tau21"              , &fj1Tau21              );
    plotTree->SetBranchAddress("fj1Tau32SD"            , &fj1Tau32SD            ); 
    plotTree->SetBranchAddress("fj1Tau21SD"            , &fj1Tau21SD            ); 
    plotTree->SetBranchAddress("fj1MSD"                , &fj1MSD                );
    plotTree->SetBranchAddress("fj1MSD_corr"           , &fj1MSD_corr           );
    plotTree->SetBranchAddress("fj1MSD_corr_jesUp"     , &fj1MSD_corr_jesUp     );
    plotTree->SetBranchAddress("fj1MSD_corr_jesDown"   , &fj1MSD_corr_jesDown   );
    plotTree->SetBranchAddress("fj1MSDScaleUp"         , &fj1MSDScaleUp         );    
    plotTree->SetBranchAddress("fj1MSDScaleDown"       , &fj1MSDScaleDown       );      
    plotTree->SetBranchAddress("fj1MSDSmeared"         , &fj1MSDSmeared         );    
    plotTree->SetBranchAddress("fj1MSDSmearedUp"       , &fj1MSDSmearedUp       );      
    plotTree->SetBranchAddress("fj1MSDSmearedDown"     , &fj1MSDSmearedDown     );        
    plotTree->SetBranchAddress("fj1Pt"                 , &fj1Pt                 );
    plotTree->SetBranchAddress("fj1PtScaleUp"          , &fj1PtScaleUp          );   
    plotTree->SetBranchAddress("fj1PtScaleDown"        , &fj1PtScaleDown        );     
    plotTree->SetBranchAddress("fj1PtSmeared"          , &fj1PtSmeared          );   
    plotTree->SetBranchAddress("fj1PtSmearedUp"        , &fj1PtSmearedUp        );     
    plotTree->SetBranchAddress("fj1PtSmearedDown"      , &fj1PtSmearedDown      );       
    plotTree->SetBranchAddress("fj1Phi"                , &fj1Phi                );
    plotTree->SetBranchAddress("fj1Eta"                , &fj1Eta                );
    plotTree->SetBranchAddress("fj1M"                  , &fj1M                  );
    plotTree->SetBranchAddress("fj1MaxCSV"             , &fj1MaxCSV             );
    plotTree->SetBranchAddress("fj1MinCSV"             , &fj1MinCSV             );
    plotTree->SetBranchAddress("fj1DoubleCSV"          , &fj1DoubleCSV          );   
    plotTree->SetBranchAddress("fj1HTTMass"            , &fj1HTTMass            ); 
    plotTree->SetBranchAddress("fj1HTTFRec"            , &fj1HTTFRec            ); 
    plotTree->SetBranchAddress("fj1SDEFrac100"         , &fj1SDEFrac100         );    
    plotTree->SetBranchAddress("fj1HighestPtGen"       , &fj1HighestPtGen       );    
    plotTree->SetBranchAddress("fj1ECFN_1_1_05", &fj1ECFN["1_1_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_05", &fj1ECFN["2_1_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_05", &fj1ECFN["3_1_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_05", &fj1ECFN["1_2_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_05", &fj1ECFN["2_2_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_05", &fj1ECFN["3_2_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_05", &fj1ECFN["1_3_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_05", &fj1ECFN["2_3_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_05", &fj1ECFN["3_3_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_05", &fj1ECFN["1_4_05"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_05", &fj1ECFN["2_4_05"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_05", &fj1ECFN["3_4_05"]);
    plotTree->SetBranchAddress("fj1ECFN_1_1_10", &fj1ECFN["1_1_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_10", &fj1ECFN["2_1_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_10", &fj1ECFN["3_1_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_10", &fj1ECFN["1_2_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_10", &fj1ECFN["2_2_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_10", &fj1ECFN["3_2_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_10", &fj1ECFN["1_3_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_10", &fj1ECFN["2_3_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_10", &fj1ECFN["3_3_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_10", &fj1ECFN["1_4_10"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_10", &fj1ECFN["2_4_10"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_10", &fj1ECFN["3_4_10"]);
    plotTree->SetBranchAddress("fj1ECFN_1_1_20", &fj1ECFN["1_1_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_20", &fj1ECFN["2_1_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_20", &fj1ECFN["3_1_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_20", &fj1ECFN["1_2_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_20", &fj1ECFN["2_2_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_20", &fj1ECFN["3_2_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_20", &fj1ECFN["1_3_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_20", &fj1ECFN["2_3_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_20", &fj1ECFN["3_3_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_20", &fj1ECFN["1_4_20"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_20", &fj1ECFN["2_4_20"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_20", &fj1ECFN["3_4_20"]);
    plotTree->SetBranchAddress("fj1ECFN_1_1_40", &fj1ECFN["1_1_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_1_40", &fj1ECFN["2_1_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_1_40", &fj1ECFN["3_1_40"]);
    plotTree->SetBranchAddress("fj1ECFN_1_2_40", &fj1ECFN["1_2_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_2_40", &fj1ECFN["2_2_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_2_40", &fj1ECFN["3_2_40"]);
    plotTree->SetBranchAddress("fj1ECFN_1_3_40", &fj1ECFN["1_3_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_3_40", &fj1ECFN["2_3_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_3_40", &fj1ECFN["3_3_40"]);
    plotTree->SetBranchAddress("fj1ECFN_1_4_40", &fj1ECFN["1_4_40"]);
    plotTree->SetBranchAddress("fj1ECFN_2_4_40", &fj1ECFN["2_4_40"]);
    plotTree->SetBranchAddress("fj1ECFN_3_4_40", &fj1ECFN["3_4_40"]);
  }
  float balanceVH=-3393;
  float dPhiHbbJet1MET=-3393;
  float dPhiHbbJet2MET=-3393;
  float MVAVar=-3393, MVAVar_jesUp=-3393, MVAVar_jesDown=-3393;
  float nAddJet,nAddJet_jesUp,nAddJet_jesDown;
  float mvaNSoft2, mvaNSoft5, mvaNSoft10;
  float dPhil1W, dPhil1b1, dPhil1b2, dPhiWb1, dPhiWb2, dPhib1b2, dEtal1W, dEtal1b1, dEtal1b2, dEtaWb1, dEtaWb2, dEtab1b2;
  float dPhil1W_jesUp, dPhiWb1_jesUp, dPhiWb2_jesUp;
  float dPhil1W_jesDown, dPhiWb1_jesDown, dPhiWb2_jesDown;
  float dPhil1fj1, dPhiWfj1, dEtal1fj1;
  float mva_nIsojet;
  float mva_lepton1Charge, mva_lepton1Flav;

  // construct histograms
  int nPlots=512;
  vector<float> xmin(nPlots), xmax(nPlots);  vector<int> nbins(nPlots);
  vector<TString> histoNames(nPlots), histoTitles(nPlots);
  
  int pMVAVar=-1;
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) { int p=0;
    histoNames[p]="pTH"                    ; histoTitles[p]="Higgs daughter dijet p_{T} [GeV]"           ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="mH"                     ; histoTitles[p]="Higgs daughter dijet mass [GeV]"            ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="WpT"                    ; histoTitles[p]="W boson p_{T} [GeV]"                        ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="deltaPhiVH"             ; histoTitles[p]="#Delta#phi(H,W) [Rad]"                      ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="pTBalanceDijetW"        ; histoTitles[p]="Higgs daughter dijet p_{T} / W boson p_{T}" ; nbins[p]=  30; xmin[p]=     0; xmax[p]=     2; p++; 
    histoNames[p]="lepton1Pt"              ; histoTitles[p]="Lepton p_{T} [GeV]"                         ; nbins[p]=  23; xmin[p]=    20; xmax[p]=   250; p++; 
    histoNames[p]="lepton1Flav"            ; histoTitles[p]="Lepton flavor"                              ; nbins[p]=   3; xmin[p]=    11; xmax[p]=    14; p++; 
    histoNames[p]="lepton1Charge"          ; histoTitles[p]="Lepton charge"                              ; nbins[p]=   3; xmin[p]=    -1; xmax[p]=     2; p++; 
    histoNames[p]="pfmet"                  ; histoTitles[p]="p_{T}^{miss} [GeV]"                         ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="pfmetphi"               ; histoTitles[p]="#phi(p_{T}^{miss}) [Rad]"                   ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="pfmetsig"               ; histoTitles[p]="Significance of p_{T}^{miss}"               ; nbins[p]=  30; xmin[p]=     0; xmax[p]=    30; p++; 
    histoNames[p]="mTW"                    ; histoTitles[p]="W transverse mass [GeV]"                    ; nbins[p]=  20; xmin[p]=     0; xmax[p]=   200; p++; 
    histoNames[p]="Hbjet1Pt"               ; histoTitles[p]="Leading H daughter jet p_{T} [GeV]"         ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="Hbjet2Pt"               ; histoTitles[p]="Trailing H daughter jet p_{T} [GeV]"        ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="dPhiHbbJet1MET"         ; histoTitles[p]="#Delta#phi(H daughter jet 1, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhiHbbJet2MET"         ; histoTitles[p]="#Delta#phi(H daughter jet 2, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="bDiscrMin"              ; histoTitles[p]="Lesser CMVA of Higgs daughter jets"         ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="bDiscrMax"              ; histoTitles[p]="Greater CMVA of Higgs daughter jets"        ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="topMassLep1Met"         ; histoTitles[p]="Reco top mass [GeV]"                        ; nbins[p]=  22; xmin[p]=   60.; xmax[p]=  500.; p++; 
    histoNames[p]="nJet"                   ; histoTitles[p]="N AK4 CHS jets pT>20"                       ; nbins[p]=   8; xmin[p]=    0.; xmax[p]=    8.; p++; 
    histoNames[p]="sumEtSoft1"             ; histoTitles[p]="#sum E_{T}(soft 1)"                         ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=  300.; p++; 
    histoNames[p]="nSoft2"                 ; histoTitles[p]="N^{soft}_{2}"                               ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=   25.; p++; 
    histoNames[p]="nSoft5"                 ; histoTitles[p]="N^{soft}_{5}"                               ; nbins[p]=  12; xmin[p]=    0.; xmax[p]=   12.; p++; 
    histoNames[p]="nSoft10"                ; histoTitles[p]="N^{soft}_{10}"                              ; nbins[p]=   8; xmin[p]=    0.; xmax[p]=    8.; p++; 
    histoNames[p]="topWBosonCosThetaCS"    ; histoTitles[p]="W boson cos #theta^{CS}"                    ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="hbbCosThetaJJ"          ; histoTitles[p]="cos #theta(bb)"                             ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="hbbCosThetaCSJ1"        ; histoTitles[p]="cos #theta^{CS} harder b"                   ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="dPhil1W"                ; histoTitles[p]="#Delta#phi(lep,W)"                          ; nbins[p]=  32; xmin[p]=    0.; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhil1b1"               ; histoTitles[p]="#Delta#phi(lep,b1)"                         ; nbins[p]=  32; xmin[p]=    0.; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhil1b2"               ; histoTitles[p]="#Delta#phi(lep,b2)"                         ; nbins[p]=  32; xmin[p]=    0.; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhiWb1"                ; histoTitles[p]="#Delta#phi(W,b1)"                           ; nbins[p]=  32; xmin[p]=    0.; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhiWb2"                ; histoTitles[p]="#Delta#phi(W,b2)"                           ; nbins[p]=  32; xmin[p]=    0.; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhib1b2"               ; histoTitles[p]="#Delta#phi(b1,b2)"                          ; nbins[p]=  32; xmin[p]=    0.; xmax[p]= 3.142; p++; 
    histoNames[p]="dEtal1W"                ; histoTitles[p]="|#Delta#eta(lep,W)|"                        ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=    5.; p++; 
    histoNames[p]="dEtal1b1"               ; histoTitles[p]="|#Delta#eta(lep,b1)|"                       ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=    5.; p++; 
    histoNames[p]="dEtal1b2"               ; histoTitles[p]="|#Delta#eta(lep,b2)|"                       ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=    5.; p++; 
    histoNames[p]="dEtaWb1"                ; histoTitles[p]="|#Delta#eta(W,b1)|"                         ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=    5.; p++; 
    histoNames[p]="dEtaWb2"                ; histoTitles[p]="|#Delta#eta(W,b2)|"                         ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=    5.; p++; 
    histoNames[p]="dEtab1b2"               ; histoTitles[p]="|#Delta#eta(b1,b2)|"                        ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=    5.; p++; 
    
    histoNames[p]="nIsojet"                ; histoTitles[p]="N central AK4 jets"                         ; nbins[p]=   7; xmin[p]=    0.; xmax[p]=    7.; p++;
    histoNames[p]="fj1Pt"                  ; histoTitles[p]="FJ p_{T} [GeV]"                             ; nbins[p]=  20; xmin[p]=  200.; xmax[p]=  600.; p++;
    histoNames[p]="fj1Eta"                 ; histoTitles[p]="FJ #eta"                                    ; nbins[p]=  20; xmin[p]=  -2.5; xmax[p]=   2.5; p++;
    histoNames[p]="maxSubjetCSV"           ; histoTitles[p]="FJ greater subjet B-tag (CSV)"              ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="minSubjetCSV"           ; histoTitles[p]="FJ lesser subjet B-tag (CSV)"               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1DoubleCSV"           ; histoTitles[p]="FJ double B-tag"                            ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++;
    histoNames[p]="fj1HTTMass"             ; histoTitles[p]="FJ HTT mass [GeV]"                          ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=  200.; p++;
    histoNames[p]="fj1HTTFRec"             ; histoTitles[p]="FJ HTT f_{rec}"                             ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=   0.4; p++;
    histoNames[p]="fj1MSD_corr"            ; histoTitles[p]="FJ soft drop mass"                          ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=  300.; p++;
    histoNames[p]="fj1Tau32"               ; histoTitles[p]="FJ #tau_{32}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau32SD"             ; histoTitles[p]="FJ #tau_{32}^{SD}"                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau21"               ; histoTitles[p]="FJ #tau_{21}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau21SD"             ; histoTitles[p]="FJ #tau_{21}^{SD}"                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="dPhil1fj1"              ; histoTitles[p]="#Delta#phi(lep,FJ) [Rad]"                   ; nbins[p]=  20; xmin[p]= 1.571; xmax[p]= 3.142; p++;
    histoNames[p]="dPhiWfj1"               ; histoTitles[p]="#Delta#phi(W,FJ) [Rad]"                     ; nbins[p]=  20; xmin[p]= 1.571; xmax[p]= 3.142; p++;
    histoNames[p]="dEtal1fj1"              ; histoTitles[p]="|#Delta#eta(lep,FJ)|"                       ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    5.; p++;
    histoNames[p]="psi021004010502"        ; histoTitles[p]="#psi(2,1.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.2  ; p++;
    histoNames[p]="psi012004010502"        ; histoTitles[p]="#psi(1,2.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi021003011002"        ; histoTitles[p]="#psi(2,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.6  ; p++;
    histoNames[p]="psi022004011002"        ; histoTitles[p]="#psi(2,2.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi020503010502"        ; histoTitles[p]="#psi(2,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.6  ; p++;
    histoNames[p]="psi022003012002"        ; histoTitles[p]="#psi(2,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.4  ; p++;
    histoNames[p]="psi012004020503"        ; histoTitles[p]="#psi(1,2.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.6  ; p++;
    histoNames[p]="psi021004020503"        ; histoTitles[p]="#psi(2,1.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 1    ; p++;
    histoNames[p]="psi032003012002"        ; histoTitles[p]="#psi(3,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 4    ; p++;
    histoNames[p]="psi012003011002"        ; histoTitles[p]="#psi(1,2.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.4  ; p++;
    histoNames[p]="psi011003010502"        ; histoTitles[p]="#psi(1,1.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.5  ; p++;
    histoNames[p]="psi031003011002"        ; histoTitles[p]="#psi(3,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3.5  ; p++;
    histoNames[p]="psi031003012002"        ; histoTitles[p]="#psi(3,1.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.4  ; p++;
    histoNames[p]="psi030503011002"        ; histoTitles[p]="#psi(3,0.5,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.5  ; p++;
    histoNames[p]="psi022004012002"        ; histoTitles[p]="#psi(2,2.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.01 ; p++;
    histoNames[p]="psi021004011002"        ; histoTitles[p]="#psi(2,1.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.04 ; p++;
    histoNames[p]="psi012004011002"        ; histoTitles[p]="#psi(1,2.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.02 ; p++;
    histoNames[p]="psi022004021003"        ; histoTitles[p]="#psi(2,2.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.8  ; p++;
    histoNames[p]="psi012003010502"        ; histoTitles[p]="#psi(1,2.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3    ; p++;
    histoNames[p]="psi022003011002"        ; histoTitles[p]="#psi(2,2.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 7    ; p++;
    histoNames[p]="psi022004031003"        ; histoTitles[p]="#psi(2,2.0,4,3,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.05 ; p++;
    histoNames[p]="psi012004030503"        ; histoTitles[p]="#psi(1,2.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.08 ; p++;
    histoNames[p]="psi021004030503"        ; histoTitles[p]="#psi(2,1.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi020504010502"        ; histoTitles[p]="#psi(2,0.5,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.08 ; p++;
    histoNames[p]="psi022004030503"        ; histoTitles[p]="#psi(2,2.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2    ; p++;
    histoNames[p]="psi011004010502"        ; histoTitles[p]="#psi(1,1.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.06 ; p++;
    histoNames[p]="psi030503010502"        ; histoTitles[p]="#psi(3,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2    ; p++;
    histoNames[p]="psi021003010502"        ; histoTitles[p]="#psi(2,1.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 5    ; p++;
    histoNames[p]="psi012004011003"        ; histoTitles[p]="#psi(1,2.0,4,1,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 1    ; p++;
    histoNames[p]="psi012004021004"        ; histoTitles[p]="#psi(1,2.0,4,2,1.0,4)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.7  ; p++;
    histoNames[p]="psi021004012004"        ; histoTitles[p]="#psi(2,1.0,4,1,2.0,4)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3    ; p++;
    histoNames[p]="psi011003020503"        ; histoTitles[p]="#psi(1,1.0,3,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.9  ; p++;
    histoNames[p]="psi020503011003"        ; histoTitles[p]="#psi(2,0.5,3,1,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3    ; p++;
    histoNames[p]="psi012003030503"        ; histoTitles[p]="#psi(1,2.0,3,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 1    ; p++;
    histoNames[p]="psi030503012003"        ; histoTitles[p]="#psi(3,0.5,3,1,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2    ; p++;
    histoNames[p]="psi010503010502"        ; histoTitles[p]="#psi(1,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.3  ; p++;
    histoNames[p]="psi020503011002"        ; histoTitles[p]="#psi(2,0.5,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.3  ; p++;
    histoNames[p]="psi011003011002"        ; histoTitles[p]="#psi(1,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.2  ; p++;
    histoNames[p]="psi020504011002"        ; histoTitles[p]="#psi(2,0.5,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.04 ; p++;
    histoNames[p]="psi012004021003"        ; histoTitles[p]="#psi(1,2.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.06 ; p++;
    histoNames[p]="psi012004012002"        ; histoTitles[p]="#psi(1,2.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.007; p++;
    histoNames[p]="psi020504020503"        ; histoTitles[p]="#psi(2,0.5,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.2  ; p++;
    histoNames[p]="psi011004011002"        ; histoTitles[p]="#psi(1,1.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.02 ; p++;
    histoNames[p]="psi021004012002"        ; histoTitles[p]="#psi(2,1.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.01 ; p++;
    histoNames[p]="psi021003012002"        ; histoTitles[p]="#psi(2,1.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.2  ; p++;
    histoNames[p]="psi011004020503"        ; histoTitles[p]="#psi(1,1.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi010504010502"        ; histoTitles[p]="#psi(1,0.5,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.05 ; p++;
    histoNames[p]="psi021004021003"        ; histoTitles[p]="#psi(2,1.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi012003012002"        ; histoTitles[p]="#psi(1,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi022004022003"        ; histoTitles[p]="#psi(2,2.0,4,2,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.04 ; p++;
    histoNames[p]="bdtValue"               ; histoTitles[p]="BDT output"                                 ; nbins[p]=  40; xmin[p]=   MVAVarType==3?-1.:0; xmax[p]=    1.; p++; 
    histoNames[p]="MVAVar"                 ; histoTitles[p]=MVAVarName; pMVAVar=p; p++;
  } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) { int p=0;
    histoNames[p]="nIsojet"                ; histoTitles[p]="N isolated central jets"                    ; nbins[p]=   7; xmin[p]=    0.; xmax[p]=    7.; p++;
    histoNames[p]="isojetNBtags"           ; histoTitles[p]="B-tagged isolated jets"                     ; nbins[p]=   5; xmin[p]=    0.; xmax[p]=    5.; p++;
    histoNames[p]="fj1Pt"                  ; histoTitles[p]="FJ p_{T} [GeV]"                             ; nbins[p]=  20; xmin[p]=  200.; xmax[p]=  600.; p++;
    histoNames[p]="fj1Eta"                 ; histoTitles[p]="FJ #eta"                                    ; nbins[p]=  20; xmin[p]=  -2.5; xmax[p]=   2.5; p++;
    histoNames[p]="maxSubjetCSV"           ; histoTitles[p]="FJ greater subjet B-tag (CSV)"              ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="minSubjetCSV"           ; histoTitles[p]="FJ lesser subjet B-tag (CSV)"               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1DoubleCSV"           ; histoTitles[p]="FJ double B-tag"                            ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++;
    histoNames[p]="fj1HTTMass"             ; histoTitles[p]="FJ HTT mass [GeV]"                          ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=  200.; p++;
    histoNames[p]="fj1HTTFRec"             ; histoTitles[p]="FJ HTT f_{rec}"                             ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=   0.4; p++;
    histoNames[p]="fj1MSD_corr"            ; histoTitles[p]="FJ soft drop mass"                          ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=  300.; p++;
    histoNames[p]="fj1Tau32"               ; histoTitles[p]="FJ #tau_{32}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau32SD"             ; histoTitles[p]="FJ #tau_{32}^{SD}"                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau21"               ; histoTitles[p]="FJ #tau_{21}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau21SD"             ; histoTitles[p]="FJ #tau_{21}^{SD}"                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="lepton1Pt"              ; histoTitles[p]="Lepton p_{T} [GeV]"                         ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  500.; p++;
    histoNames[p]="lepton1Flav"            ; histoTitles[p]="Lepton flavor"                              ; nbins[p]=   3; xmin[p]=    11; xmax[p]=    14; p++; 
    histoNames[p]="lepton1Charge"          ; histoTitles[p]="Lepton charge"                              ; nbins[p]=   3; xmin[p]=    -1; xmax[p]=     2; p++; 
    histoNames[p]="deltaPhiVH"             ; histoTitles[p]="#Delta#phi(W,H) [Rad]"                      ; nbins[p]=  20; xmin[p]= 1.571; xmax[p]= 3.142; p++;
    histoNames[p]="deltaPhiLep1Met"        ; histoTitles[p]="#Delta#phi(lepton,p_{T}^{miss}) [Rad]"      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3.142; p++;
    histoNames[p]="topWBosonPt"            ; histoTitles[p]="W boson p_{T} [GeV]"                        ; nbins[p]=  20; xmin[p]=  140.; xmax[p]=  540.; p++;
    histoNames[p]="mT"                     ; histoTitles[p]="W boson m_{T} [GeV]"                        ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  200.; p++;
    histoNames[p]="pfmet"                  ; histoTitles[p]="p_{T}^{miss} [GeV]"                         ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  500.; p++;
    histoNames[p]="pfmetphi"               ; histoTitles[p]="#phi(p_{T}^{miss}) [Rad]"                   ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="pfmetsig"               ; histoTitles[p]="Significance of p_{T}^{miss}"               ; nbins[p]=  30; xmin[p]=     0; xmax[p]=    30; p++; 
    histoNames[p]="dPhil1W"                ; histoTitles[p]="#Delta#phi(lep,W) [Rad]"                    ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 1.571; p++;
    histoNames[p]="dPhil1fj1"              ; histoTitles[p]="#Delta#phi(lep,FJ) [Rad]"                   ; nbins[p]=  20; xmin[p]= 1.571; xmax[p]= 3.142; p++;
    histoNames[p]="dPhiWfj1"               ; histoTitles[p]="#Delta#phi(W,FJ) [Rad]"                     ; nbins[p]=  20; xmin[p]= 1.571; xmax[p]= 3.142; p++;
    histoNames[p]="dEtal1fj1"              ; histoTitles[p]="|#Delta#eta(lep,FJ)|"                       ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    5.; p++;
    histoNames[p]="psi021004010502"        ; histoTitles[p]="#psi(2,1.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.4  ; p++;
    histoNames[p]="psi012004010502"        ; histoTitles[p]="#psi(1,2.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi021003011002"        ; histoTitles[p]="#psi(2,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.6  ; p++;
    histoNames[p]="psi022004011002"        ; histoTitles[p]="#psi(2,2.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi020503010502"        ; histoTitles[p]="#psi(2,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.6  ; p++;
    histoNames[p]="psi022003012002"        ; histoTitles[p]="#psi(2,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.4  ; p++;
    histoNames[p]="psi012004020503"        ; histoTitles[p]="#psi(1,2.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.6  ; p++;
    histoNames[p]="psi032003012002"        ; histoTitles[p]="#psi(3,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 4    ; p++;
    histoNames[p]="psi012003011002"        ; histoTitles[p]="#psi(1,2.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.4  ; p++;
    histoNames[p]="psi011003010502"        ; histoTitles[p]="#psi(1,1.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.5  ; p++;
    histoNames[p]="psi031003011002"        ; histoTitles[p]="#psi(3,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 5    ; p++;
    histoNames[p]="psi031003012002"        ; histoTitles[p]="#psi(3,1.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.4  ; p++;
    histoNames[p]="psi030503011002"        ; histoTitles[p]="#psi(3,0.5,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.5  ; p++;
    histoNames[p]="psi022004012002"        ; histoTitles[p]="#psi(2,2.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.01 ; p++;
    histoNames[p]="psi021004011002"        ; histoTitles[p]="#psi(2,1.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.04 ; p++;
    histoNames[p]="psi012004011002"        ; histoTitles[p]="#psi(1,2.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.02 ; p++;
    histoNames[p]="psi022004021003"        ; histoTitles[p]="#psi(2,2.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.8  ; p++;
    histoNames[p]="psi012003010502"        ; histoTitles[p]="#psi(1,2.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3    ; p++;
    histoNames[p]="psi022003011002"        ; histoTitles[p]="#psi(2,2.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 7    ; p++;
    histoNames[p]="psi022004031003"        ; histoTitles[p]="#psi(2,2.0,4,3,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.05 ; p++;
    histoNames[p]="psi012004030503"        ; histoTitles[p]="#psi(1,2.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.08 ; p++;
    histoNames[p]="psi021004030503"        ; histoTitles[p]="#psi(2,1.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi020504010502"        ; histoTitles[p]="#psi(2,0.5,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.08 ; p++;
    histoNames[p]="psi022004030503"        ; histoTitles[p]="#psi(2,2.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2    ; p++;
    histoNames[p]="psi011004010502"        ; histoTitles[p]="#psi(1,1.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.06 ; p++;
    histoNames[p]="psi030503010502"        ; histoTitles[p]="#psi(3,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3    ; p++;
    histoNames[p]="psi021003010502"        ; histoTitles[p]="#psi(2,1.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 5    ; p++;
    histoNames[p]="psi012004021004"        ; histoTitles[p]="#psi(1,2.0,4,2,1.0,4)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.7  ; p++;
    histoNames[p]="psi021004012004"        ; histoTitles[p]="#psi(2,1.0,4,1,2.0,4)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3    ; p++;
    histoNames[p]="psi011003020503"        ; histoTitles[p]="#psi(1,1.0,3,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.9  ; p++;
    histoNames[p]="psi020503011003"        ; histoTitles[p]="#psi(2,0.5,3,1,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3    ; p++;
    histoNames[p]="psi012003030503"        ; histoTitles[p]="#psi(1,2.0,3,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 1    ; p++;
    histoNames[p]="psi030503012003"        ; histoTitles[p]="#psi(3,0.5,3,1,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2    ; p++;
    histoNames[p]="psi010503010502"        ; histoTitles[p]="#psi(1,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.3  ; p++;
    histoNames[p]="psi020503011002"        ; histoTitles[p]="#psi(2,0.5,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.3  ; p++;
    histoNames[p]="psi011003011002"        ; histoTitles[p]="#psi(1,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.2  ; p++;
    histoNames[p]="psi020504011002"        ; histoTitles[p]="#psi(2,0.5,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.04 ; p++;
    histoNames[p]="psi012004021003"        ; histoTitles[p]="#psi(1,2.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.06 ; p++;
    histoNames[p]="psi012004012002"        ; histoTitles[p]="#psi(1,2.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.007; p++;
    histoNames[p]="psi020504020503"        ; histoTitles[p]="#psi(2,0.5,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.2  ; p++;
    histoNames[p]="psi011004011002"        ; histoTitles[p]="#psi(1,1.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.02 ; p++;
    histoNames[p]="psi021004012002"        ; histoTitles[p]="#psi(2,1.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.01 ; p++;
    histoNames[p]="psi021003012002"        ; histoTitles[p]="#psi(2,1.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.2  ; p++;
    histoNames[p]="psi011004020503"        ; histoTitles[p]="#psi(1,1.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi010504010502"        ; histoTitles[p]="#psi(1,0.5,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.05 ; p++;
    histoNames[p]="psi021004021003"        ; histoTitles[p]="#psi(2,1.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi012003012002"        ; histoTitles[p]="#psi(1,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.1  ; p++;
    histoNames[p]="psi022004022003"        ; histoTitles[p]="#psi(2,2.0,4,2,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 0.04 ; p++;

    histoNames[p]="psi012002011002"        ; histoTitles[p]="#psi(1,2.0,2,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=   2.5; xmax[p]=   7.5; p++;
    histoNames[p]="psi014003022003"        ; histoTitles[p]="#psi(1,4.0,3,2,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1 ; p++;
    histoNames[p]="psi031003014003"        ; histoTitles[p]="#psi(3,1.0,3,1,4.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    5 ; p++;
    histoNames[p]="psi031003022003"        ; histoTitles[p]="#psi(3,1.0,3,2,2.0,3)"                      ; nbins[p]=  20; xmin[p]=   0.3; xmax[p]=   1.5; p++;
    histoNames[p]="psi032003034003"        ; histoTitles[p]="#psi(3,2.0,3,3,4.0,3)"                      ; nbins[p]=  20; xmin[p]=    0 ; xmax[p]=  0.25; p++;
    histoNames[p]="psi012004011003"        ; histoTitles[p]="#psi(1,2.0,4,1,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 1    ; p++;
    histoNames[p]="psi014004012003"        ; histoTitles[p]="#psi(1,4.0,4,1,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    4 ; p++;
    histoNames[p]="psi020504010503"        ; histoTitles[p]="#psi(2,0.5,4,1,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    4 ; p++;
    histoNames[p]="psi021004011003"        ; histoTitles[p]="#psi(2,1.0,4,1,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    5 ; p++;
    histoNames[p]="psi021004020503"        ; histoTitles[p]="#psi(2,1.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    2 ; p++;
    histoNames[p]="psi022004012003"        ; histoTitles[p]="#psi(2,2.0,4,1,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    6 ; p++;
    histoNames[p]="bdtValue"               ; histoTitles[p]="BDT output"                                 ; nbins[p]=  40; xmin[p]=   MVAVarType==3?-1.:0; xmax[p]=    1.; p++; 
    histoNames[p]="MVAVar"                 ; histoTitles[p]=MVAVarName; pMVAVar=p; p++;
  }
  assert(pMVAVar!=-1);

  TH1D *histos[512][nPlotCategories];
  for(int p=0; p<nPlots; p++) {
    if(histoNames[p]=="") continue;
    for(theCategory=0; theCategory<nPlotCategories; theCategory++) {

      if(histoNames[p]=="MVAVar") {
        histos[p][theCategory] = new TH1D(Form("histo%d",theCategory), histoTitles[p], MVAbins.size()-1, MVAbins.data());
      } else
        histos[p][theCategory] = new TH1D(Form("histo%d",theCategory), histoTitles[p], nbins[p], xmin[p], xmax[p]);
      histos[p][theCategory]->Sumw2();
      histos[p][theCategory]->SetDirectory(0);
    }
  }
  TH1D *histo_VHCorrUp    [nPlotCategories];
  TH1D *histo_VHCorrDown  [nPlotCategories];
  //TH1D *histo_pdfUp    [nPlotCategories];
  //TH1D *histo_pdfDown  [nPlotCategories];
  TH1D *histo_QCDr1f2  [nPlotCategories];
  TH1D *histo_QCDr1f5  [nPlotCategories];
  TH1D *histo_QCDr2f1  [nPlotCategories];
  TH1D *histo_QCDr2f2  [nPlotCategories];
  TH1D *histo_QCDr5f1  [nPlotCategories];
  TH1D *histo_QCDr5f5  [nPlotCategories];
  TH1D *histo_NLOQCDrUp   [nPlotCategories];
  TH1D *histo_NLOQCDrDown [nPlotCategories];
  TH1D *histo_NLOQCDfUp   [nPlotCategories];
  TH1D *histo_NLOQCDfDown [nPlotCategories];
  TH1D *histo_QCDrScaleUp   [nPlotCategories];
  TH1D *histo_QCDrScaleDown [nPlotCategories];
  TH1D *histo_QCDfScaleUp   [nPlotCategories];
  TH1D *histo_QCDfScaleDown [nPlotCategories];
  TH1D *histo_eleSFUp  [nPlotCategories],*histo_eleSFDown[nPlotCategories];
  TH1D *histo_muSFUp   [nPlotCategories],*histo_muSFDown [nPlotCategories];
  TH1D *histo_cmvaJESUp     [5][3][nPlotCategories], *histo_cmvaJESDown     [5][3][nPlotCategories]; 
  TH1D *histo_cmvaLFUp      [5][3][nPlotCategories], *histo_cmvaLFDown      [5][3][nPlotCategories]; 
  TH1D *histo_cmvaHFUp      [5][3][nPlotCategories], *histo_cmvaHFDown      [5][3][nPlotCategories]; 
  TH1D *histo_cmvaHFStats1Up[5][3][nPlotCategories], *histo_cmvaHFStats1Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaHFStats2Up[5][3][nPlotCategories], *histo_cmvaHFStats2Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaLFStats1Up[5][3][nPlotCategories], *histo_cmvaLFStats1Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaLFStats2Up[5][3][nPlotCategories], *histo_cmvaLFStats2Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaCErr1Up   [5][3][nPlotCategories], *histo_cmvaCErr1Down   [5][3][nPlotCategories]; 
  TH1D *histo_cmvaCErr2Up   [5][3][nPlotCategories], *histo_cmvaCErr2Down   [5][3][nPlotCategories]; 
  TH1D *histo_VGluUp   [nPlotCategories];
  TH1D *histo_VGluDown [nPlotCategories];
  TH1D *histo_wptCorrUp   [nPlotCategories];
  TH1D *histo_wptCorrDown [nPlotCategories];
  TH1D *histo_wptCorrNone [nPlotCategories];
  TH1D *histo_jesUp   [nPlotCategories];
  TH1D *histo_jesDown [nPlotCategories];
  for(theCategory=0; theCategory<nPlotCategories; theCategory++) {
    histo_VHCorrUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_VHCorrUp"  , theCategory)); histo_VHCorrUp  [theCategory]->SetDirectory(0);
    histo_VHCorrDown[theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_VHCorrDown", theCategory)); histo_VHCorrDown[theCategory]->SetDirectory(0);
    //histo_pdfUp      [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_pdfUp"    , theCategory)); histo_pdfUp     [theCategory]->SetDirectory(0);
    //histo_pdfDown    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_pdfDown"  , theCategory)); histo_pdfDown   [theCategory]->SetDirectory(0);
    histo_QCDr1f2    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr1f2"  , theCategory)); histo_QCDr1f2   [theCategory]->SetDirectory(0);
    histo_QCDr1f5    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr1f5"  , theCategory)); histo_QCDr1f5   [theCategory]->SetDirectory(0);
    histo_QCDr2f1    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr2f1"  , theCategory)); histo_QCDr2f1   [theCategory]->SetDirectory(0);
    histo_QCDr2f2    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr2f2"  , theCategory)); histo_QCDr2f2   [theCategory]->SetDirectory(0);
    histo_QCDr5f1    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr5f1"  , theCategory)); histo_QCDr5f1   [theCategory]->SetDirectory(0);
    histo_QCDr5f5    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr5f5"  , theCategory)); histo_QCDr5f5   [theCategory]->SetDirectory(0);
    histo_QCDrScaleUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDrScaleUp"   , theCategory)); histo_QCDrScaleUp  [theCategory]->SetDirectory(0);
    histo_QCDrScaleDown[theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDrScaleDown" , theCategory)); histo_QCDrScaleDown[theCategory]->SetDirectory(0);
    histo_QCDfScaleUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDfScaleUp"   , theCategory)); histo_QCDfScaleUp  [theCategory]->SetDirectory(0);
    histo_QCDfScaleDown[theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDfScaleDown" , theCategory)); histo_QCDfScaleDown[theCategory]->SetDirectory(0);
    histo_NLOQCDrUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_NLOQCDrUp"   , theCategory)); histo_NLOQCDrUp  [theCategory]->SetDirectory(0);
    histo_NLOQCDrDown[theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_NLOQCDrDown" , theCategory)); histo_NLOQCDrDown[theCategory]->SetDirectory(0);
    histo_NLOQCDfUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_NLOQCDfUp"   , theCategory)); histo_NLOQCDfUp  [theCategory]->SetDirectory(0);
    histo_NLOQCDfDown[theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_NLOQCDfDown" , theCategory)); histo_NLOQCDfDown[theCategory]->SetDirectory(0);
    histo_eleSFUp    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_eleSFUp"  , theCategory)); histo_eleSFUp   [theCategory]->SetDirectory(0);
    histo_eleSFDown  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_eleSFDown", theCategory)); histo_eleSFDown [theCategory]->SetDirectory(0);
    histo_muSFUp     [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_muSFUp"   , theCategory)); histo_muSFUp    [theCategory]->SetDirectory(0);
    histo_muSFDown   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_muSFDown" , theCategory)); histo_muSFDown  [theCategory]->SetDirectory(0);
    for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
      histo_cmvaJESUp       [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaJESUp_pt%d_eta%d"       , theCategory, iPt, iEta)); histo_cmvaJESUp        [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFUp        [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFUp_pt%d_eta%d"        , theCategory, iPt, iEta)); histo_cmvaLFUp         [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFUp        [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFUp_pt%d_eta%d"        , theCategory, iPt, iEta)); histo_cmvaHFUp         [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats1Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats1Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaHFStats1Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats2Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats2Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaHFStats2Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats1Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats1Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaLFStats1Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats2Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats2Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaLFStats2Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr1Up     [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr1Up_pt%d_eta%d"     , theCategory, iPt, iEta)); histo_cmvaCErr1Up      [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr2Up     [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr2Up_pt%d_eta%d"     , theCategory, iPt, iEta)); histo_cmvaCErr2Up      [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaJESDown     [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaJESDown_pt%d_eta%d"     , theCategory, iPt, iEta)); histo_cmvaJESDown      [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFDown      [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFDown_pt%d_eta%d"      , theCategory, iPt, iEta)); histo_cmvaLFDown       [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFDown      [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFDown_pt%d_eta%d"      , theCategory, iPt, iEta)); histo_cmvaHFDown       [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats1Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats1Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaHFStats1Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats2Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats2Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaHFStats2Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats1Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats1Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaLFStats1Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats2Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats2Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaLFStats2Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr1Down   [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr1Down_pt%d_eta%d"   , theCategory, iPt, iEta)); histo_cmvaCErr1Down    [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr2Down   [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr2Down_pt%d_eta%d"   , theCategory, iPt, iEta)); histo_cmvaCErr2Down    [iPt][iEta][theCategory]->SetDirectory(0);
    }
    histo_VGluUp      [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_VGluUp"      , theCategory)); histo_VGluUp      [theCategory]->SetDirectory(0);
    histo_VGluDown    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_VGluDown"    , theCategory)); histo_VGluDown    [theCategory]->SetDirectory(0);
    histo_wptCorrUp   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_wptCorrUp"   , theCategory)); histo_wptCorrUp   [theCategory]->SetDirectory(0);
    histo_wptCorrDown [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_wptCorrDown" , theCategory)); histo_wptCorrDown [theCategory]->SetDirectory(0);
    histo_wptCorrNone [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_wptCorrNone" , theCategory)); histo_wptCorrNone [theCategory]->SetDirectory(0);
    histo_jesUp   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_jesUp"   , theCategory)); histo_jesUp   [theCategory]->SetDirectory(0);
    histo_jesDown [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_jesDown" , theCategory)); histo_jesDown [theCategory]->SetDirectory(0);
  }
  // done constructing histograms

  // Initialize mva reader
  std::map <string, float*> mvaInputRefs, mvaInputRefs_jesUp, mvaInputRefs_jesDown;
  string str_maxSubjetCSV = "fj1MaxCSV!=fj1MaxCSV?0.:TMath::Min(1.,fj1MaxCSV)";
  string str_minSubjetCSV = "fj1MinCSV!=fj1MinCSV?0.:TMath::Min(1.,fj1MinCSV)";
  mvaInputRefs["hbbDijetMass"                                          ]=&hbbDijetMass          ;    
  mvaInputRefs["hbbDijetPt"                                            ]=&hbbDijetPt            ;    
  mvaInputRefs["topWBosonPt"                                           ]=&topWBosonPt           ;    
  mvaInputRefs["bDiscrMin"                                             ]=&bDiscrMin             ;    
  mvaInputRefs["topMassLep1Met"                                        ]=&topMassLep1Met        ;    
  mvaInputRefs["nJet-2"                                                ]=&nAddJet               ;    
  mvaInputRefs["deltaPhiVH"                                            ]=&deltaPhiVH            ;    
  mvaInputRefs["deltaPhiLep1Met"                                       ]=&deltaPhiLep1Met       ;    
  mvaInputRefs["mT"                                                    ]=&mT                    ;    
  mvaInputRefs["pfmet"                                                 ]=&pfmet                 ;    
  mvaInputRefs["pfmetsig"                                              ]=&pfmetsig              ;    
  mvaInputRefs["nSoft5"                                                ]=&mvaNSoft5             ;    
  mvaInputRefs["lepton1Pt"                                             ]=&lepton1Pt             ;
  mvaInputRefs["lepton1Flav"                                           ]=&mva_lepton1Flav       ;
  mvaInputRefs["lepton1Charge"                                         ]=&mva_lepton1Charge     ;
  mvaInputRefs["bDiscrMax"                                             ]=&bDiscrMax             ;    
  mvaInputRefs["sumEtSoft1"                                            ]=&sumEtSoft1            ;    
  mvaInputRefs["nSoft2"                                                ]=&mvaNSoft2             ;    
  mvaInputRefs["nSoft10"                                               ]=&mvaNSoft10            ;    
  mvaInputRefs["topWBosonCosThetaCS"                                   ]=&topWBosonCosThetaCS   ;    
  mvaInputRefs["hbbCosThetaJJ"                                         ]=&hbbCosThetaJJ         ;    
  mvaInputRefs["hbbCosThetaCSJ1"                                       ]=&hbbCosThetaCSJ1       ;    
  mvaInputRefs["hbbJet1Pt"                                             ]=&hbbJet1Pt             ;    
  mvaInputRefs["hbbJet2Pt"                                             ]=&hbbJet2Pt             ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-topWBosonPhi))"   ]=&dPhil1W               ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-hbbJet1Phi))"     ]=&dPhil1b1              ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-hbbJet2Phi))"     ]=&dPhil1b2              ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet1Phi))"   ]=&dPhiWb1               ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet2Phi))"   ]=&dPhiWb2               ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(hbbJet1Phi-hbbJet2Phi))"     ]=&dPhib1b2              ;    
  mvaInputRefs["fabs(lepton1Eta-topWBosonEta)"                         ]=&dEtal1W               ;    
  mvaInputRefs["fabs(lepton1Eta-hbbJet1Eta)"                           ]=&dEtal1b1              ;    
  mvaInputRefs["fabs(lepton1Eta-hbbJet2Eta)"                           ]=&dEtal1b2              ;    
  mvaInputRefs["fabs(topWBosonEta-hbbJet1Eta)"                         ]=&dEtaWb1               ;    
  mvaInputRefs["fabs(topWBosonEta-hbbJet2Eta)"                         ]=&dEtaWb2               ;    
  mvaInputRefs["fabs(hbbJet1Eta-hbbJet2Eta)"                           ]=&dEtab1b2              ;    
  mvaInputRefs["nIsojet"                                               ]=&mva_nIsojet           ;
  mvaInputRefs["fj1Pt"                                                 ]=&fj1Pt                 ;
  mvaInputRefs["(nFatjet==0)*(-5)+(nFatjet!=0)*(fj1Eta)"               ]=&adjustedFatjetEta     ;
  mvaInputRefs[str_maxSubjetCSV                                        ]=&fj1MaxCSV             ;
  mvaInputRefs[str_minSubjetCSV                                        ]=&fj1MinCSV             ;
  mvaInputRefs["fj1DoubleCSV"                                          ]=&fj1DoubleCSV          ;
  mvaInputRefs["fj1HTTMass"                                            ]=&fj1HTTMass            ;
  mvaInputRefs["fj1HTTFRec"                                            ]=&fj1HTTFRec            ;
  mvaInputRefs["fj1MSD_corr"                                           ]=&fj1MSD_corr           ;
  mvaInputRefs["fj1Tau32"                                              ]=&fj1Tau32              ;
  mvaInputRefs["fj1Tau32SD"                                            ]=&fj1Tau32SD            ;
  mvaInputRefs["fj1Tau21"                                              ]=&fj1Tau21              ;
  mvaInputRefs["fj1Tau21SD"                                            ]=&fj1Tau21SD            ;
  mvaInputRefs["(nFatjet==0)*(-1)+(nFatjet!=0)*fabs(TVector2::Phi_mpi_pi(lepton1Phi-fj1Phi))"]=&dPhil1fj1;
  mvaInputRefs["(nFatjet==0)*(-1)+(nFatjet!=0)*fabs(TVector2::Phi_mpi_pi(topWBosonPhi-fj1Phi))"]=&dPhiWfj1;
  mvaInputRefs["(nFatjet==0)*(-1)+(nFatjet!=0)*fabs(lepton1Eta-fj1Eta)"]=&dEtal1fj1;
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)"]=&mvaPsi["021004010502"]; 
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)"]=&mvaPsi["012004010502"]; 
  mvaInputRefs["fj1ECFN_2_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)"]=&mvaPsi["021003011002"]; 
  mvaInputRefs["fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),4.00)"]=&mvaPsi["022004011002"]; 
  mvaInputRefs["fj1ECFN_2_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)"]=&mvaPsi["020503010502"]; 
  mvaInputRefs["fj1ECFN_2_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),2.00)"]=&mvaPsi["022003012002"]; 
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_05),2.00)"]=&mvaPsi["012004020503"]; 
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_2_3_05),2.00)"]=&mvaPsi["021004020503"]; 
  mvaInputRefs["fj1ECFN_3_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),3.00)"]=&mvaPsi["032003012002"]; 
  mvaInputRefs["fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)"]=&mvaPsi["012003011002"]; 
  mvaInputRefs["fj1ECFN_1_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)"]=&mvaPsi["011003010502"]; 
  mvaInputRefs["fj1ECFN_3_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),3.00)"]=&mvaPsi["031003011002"]; 
  mvaInputRefs["fj1ECFN_3_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.50)"]=&mvaPsi["031003012002"]; 
  mvaInputRefs["fj1ECFN_3_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.50)"]=&mvaPsi["030503011002"]; 
  mvaInputRefs["fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),2.00)"]=&mvaPsi["022004012002"]; 
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)"]=&mvaPsi["021004011002"]; 
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)"]=&mvaPsi["012004011002"]; 
  mvaInputRefs["fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_10),2.00)"]=&mvaPsi["022004021003"]; 
  mvaInputRefs["fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)"]=&mvaPsi["012003010502"]; 
  mvaInputRefs["fj1ECFN_2_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),4.00)"]=&mvaPsi["022003011002"]; 
  mvaInputRefs["fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_3_3_10),1.33)"]=&mvaPsi["022004031003"]; 
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_3_3_05),1.33)"]=&mvaPsi["012004030503"]; 
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_3_3_05),1.33)"]=&mvaPsi["021004030503"]; 
  mvaInputRefs["fj1ECFN_2_4_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)"]=&mvaPsi["020504010502"]; 
  mvaInputRefs["fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_3_3_05),2.67)"]=&mvaPsi["022004030503"]; 
  mvaInputRefs["fj1ECFN_1_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)"]=&mvaPsi["011004010502"]; 
  mvaInputRefs["fj1ECFN_3_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),3.00)"]=&mvaPsi["030503010502"]; 
  mvaInputRefs["fj1ECFN_2_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)"]=&mvaPsi["021003010502"]; 
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_2_4_10),1.00)"]=&mvaPsi["012004021004"]; 
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_4_20),1.00)"]=&mvaPsi["021004012004"]; 
  mvaInputRefs["fj1ECFN_1_3_10/pow(TMath::Max(0.,fj1ECFN_2_3_05),1.00)"]=&mvaPsi["011003020503"]; 
  mvaInputRefs["fj1ECFN_2_3_05/pow(TMath::Max(0.,fj1ECFN_1_3_10),1.00)"]=&mvaPsi["020503011003"]; 
  mvaInputRefs["fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_3_3_05),1.33)"]=&mvaPsi["012003030503"]; 
  mvaInputRefs["fj1ECFN_3_3_05/pow(TMath::Max(0.,fj1ECFN_1_3_20),0.75)"]=&mvaPsi["030503012003"]; 
  mvaInputRefs["fj1ECFN_1_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),1.00)"]=&mvaPsi["010503010502"]; 
  mvaInputRefs["fj1ECFN_2_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)"]=&mvaPsi["020503011002"]; 
  mvaInputRefs["fj1ECFN_1_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)"]=&mvaPsi["011003011002"]; 
  mvaInputRefs["fj1ECFN_2_4_05/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)"]=&mvaPsi["020504011002"]; 
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_10),1.00)"]=&mvaPsi["012004021003"]; 
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)"]=&mvaPsi["012004012002"]; 
  mvaInputRefs["fj1ECFN_2_4_05/pow(TMath::Max(0.,fj1ECFN_2_3_05),1.00)"]=&mvaPsi["020504020503"]; 
  mvaInputRefs["fj1ECFN_1_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)"]=&mvaPsi["011004011002"]; 
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)"]=&mvaPsi["021004012002"]; 
  mvaInputRefs["fj1ECFN_2_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)"]=&mvaPsi["021003012002"]; 
  mvaInputRefs["fj1ECFN_1_4_10/pow(TMath::Max(0.,fj1ECFN_2_3_05),1.00)"]=&mvaPsi["011004020503"]; 
  mvaInputRefs["fj1ECFN_1_4_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),1.00)"]=&mvaPsi["010504010502"]; 
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_2_3_10),1.00)"]=&mvaPsi["021004021003"]; 
  mvaInputRefs["fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)"]=&mvaPsi["012003012002"]; 
  mvaInputRefs["fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_20),1.00)"]=&mvaPsi["022004022003"]; 

  mvaInputRefs["fj1ECFN_1_2_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)"]=&mvaPsi["012002011002"];
  mvaInputRefs["fj1ECFN_1_3_40/pow(TMath::Max(0.,fj1ECFN_2_3_20),1.00)"]=&mvaPsi["014003022003"];
  mvaInputRefs["fj1ECFN_3_3_10/pow(TMath::Max(0.,fj1ECFN_1_3_40),0.75)"]=&mvaPsi["031003014003"];
  mvaInputRefs["fj1ECFN_3_3_10/pow(TMath::Max(0.,fj1ECFN_2_3_20),0.75)"]=&mvaPsi["031003022003"];
  mvaInputRefs["fj1ECFN_3_3_20/pow(TMath::Max(0.,fj1ECFN_3_3_40),0.50)"]=&mvaPsi["032003034003"];
  mvaInputRefs["fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_3_10),2.00)"]=&mvaPsi["012004011003"];
  mvaInputRefs["fj1ECFN_1_4_40/pow(TMath::Max(0.,fj1ECFN_1_3_20),2.00)"]=&mvaPsi["014004012003"];
  mvaInputRefs["fj1ECFN_2_4_05/pow(TMath::Max(0.,fj1ECFN_1_3_05),2.00)"]=&mvaPsi["020504010503"];
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_3_10),2.00)"]=&mvaPsi["021004011003"];
  mvaInputRefs["fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_2_3_05),2.00)"]=&mvaPsi["021004020503"];
  mvaInputRefs["fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_1_3_20),2.00)"]=&mvaPsi["022004012003"];
 
  mvaInputRefs_jesUp = mvaInputRefs;
  mvaInputRefs_jesUp["hbbDijetMass"                                          ]=&hbbDijetMassUp        ;    
  mvaInputRefs_jesUp["hbbDijetPt"                                            ]=&hbbDijetPtUp          ;    
  mvaInputRefs_jesUp["topWBosonPt"                                           ]=&topWBosonPt_jesUp     ;    
  mvaInputRefs_jesUp["topMassLep1Met"                                        ]=&topMassLep1Met_jesUp  ;    
  mvaInputRefs_jesUp["nJet-2"                                                ]=&nAddJet_jesUp         ;    
  mvaInputRefs_jesUp["deltaPhiVH"                                            ]=&deltaPhiVH_jesUp      ;    
  mvaInputRefs_jesUp["mT"                                                    ]=&mT_jesUp              ;    
  mvaInputRefs_jesUp["pfmet"                                                 ]=&pfmetUp               ;    
  mvaInputRefs_jesUp["hbbJet1Pt"                                             ]=&hbbJet1PtUp           ;    
  mvaInputRefs_jesUp["hbbJet2Pt"                                             ]=&hbbJet2PtUp           ;    
  mvaInputRefs_jesUp["fabs(TVector2::Phi_mpi_pi(lepton1Phi-topWBosonPhi))"   ]=&dPhil1W_jesUp         ;    
  mvaInputRefs_jesUp["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet1Phi))"   ]=&dPhiWb1_jesUp         ;    
  mvaInputRefs_jesUp["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet2Phi))"   ]=&dPhiWb2_jesUp         ;    
  mvaInputRefs_jesUp["fj1Pt"                                                 ]=&fj1PtScaleUp          ;
  mvaInputRefs_jesUp["fj1MSD_corr"                                           ]=&fj1MSD_corr_jesUp     ;

  mvaInputRefs_jesDown = mvaInputRefs;
  mvaInputRefs_jesDown["hbbDijetMass"                                          ]=&hbbDijetMassDown      ;    
  mvaInputRefs_jesDown["hbbDijetPt"                                            ]=&hbbDijetPtDown        ;    
  mvaInputRefs_jesDown["topWBosonPt"                                           ]=&topWBosonPt_jesDown   ;    
  mvaInputRefs_jesDown["topMassLep1Met"                                        ]=&topMassLep1Met_jesDown;    
  mvaInputRefs_jesDown["nJet-2"                                                ]=&nAddJet_jesDown       ;    
  mvaInputRefs_jesDown["deltaPhiVH"                                            ]=&deltaPhiVH_jesUp      ;    
  mvaInputRefs_jesDown["mT"                                                    ]=&mT_jesDown            ;    
  mvaInputRefs_jesDown["pfmet"                                                 ]=&pfmetDown             ;    
  mvaInputRefs_jesDown["hbbJet1Pt"                                             ]=&hbbJet1PtDown         ;    
  mvaInputRefs_jesDown["hbbJet2Pt"                                             ]=&hbbJet2PtDown         ;    
  mvaInputRefs_jesDown["fabs(TVector2::Phi_mpi_pi(lepton1Phi-topWBosonPhi))"   ]=&dPhil1W_jesDown       ;    
  mvaInputRefs_jesDown["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet1Phi))"   ]=&dPhiWb1_jesDown       ;    
  mvaInputRefs_jesDown["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet2Phi))"   ]=&dPhiWb2_jesDown       ;    
  mvaInputRefs_jesDown["fj1Pt"                                                 ]=&fj1PtScaleDown        ;
  mvaInputRefs_jesDown["fj1MSD_corr"                                           ]=&fj1MSD_corr_jesDown   ;
  TMVA::Reader *reader=0, *reader_jesUp=0, *reader_jesDown=0;
  // to do - not hardcoded
  if(MVAVarType==2 || MVAVarType==3) {
    reader=new TMVA::Reader();
    reader_jesUp=new TMVA::Reader("Silent");
    reader_jesDown=new TMVA::Reader("Silent");
    assert(bdtWeights!="");
    bdtWeights = "MitVHBBAnalysis/weights/"+bdtWeights;
  }
  if(bdtWeights!="") {
    vector<pair<string, string> > varLabelExprs=parseTmvaXmlFile(string(bdtWeights.Data()));
    for(unsigned iVar=0; iVar<(unsigned)varLabelExprs.size(); iVar++)
      if(mvaInputRefs.find( varLabelExprs[iVar].second ) != mvaInputRefs.end()) {
        printf("\tAdding variable #%d with label \"%s\"\n", iVar, varLabelExprs[iVar].first.c_str());
        reader        ->AddVariable( varLabelExprs[iVar].second, mvaInputRefs        [varLabelExprs[iVar].second]);
        reader_jesUp  ->AddVariable( varLabelExprs[iVar].second, mvaInputRefs_jesUp  [varLabelExprs[iVar].second]);
        reader_jesDown->AddVariable( varLabelExprs[iVar].second, mvaInputRefs_jesDown[varLabelExprs[iVar].second]);
      }
    reader        ->BookMVA("BDT", bdtWeights.Data());
    reader_jesUp  ->BookMVA("BDT", bdtWeights.Data());
    reader_jesDown->BookMVA("BDT", bdtWeights.Data());
  }

  float sf_training=1;
  // train on the 30% of events with event number ending in 0,1,2
  // reweight the rest to have 43% more weight
  if((selection==kWHSR || selection==kWHFJSR)&& (MVAVarType==2 || MVAVarType==3)) sf_training=1.4286;

  // begin plot tree loop
  Long64_t nentries = plotTree->GetEntries();
  Long64_t oneTenth = nentries/10;
  for(Long64_t ientry=0; ientry<nentries; ientry++) {
    if(ientry%oneTenth==0) printf("######## Reading entry %lld/%lld ########################################################\n",ientry,nentries); 
    plotTree->GetBranch("typeLepSel")->GetEntry(ientry); if(theLepSel!=-1 && typeLepSel!=theLepSel) continue;
    plotTree->GetBranch("nMinusOneBits")->GetEntry(ientry); if((nMinusOneBits & selection)==0) continue;
    plotTree->GetBranch("theCategory")->GetEntry(ientry); if(theCategory>=nPlotCategories) continue;
    if(isBlinded && theCategory==kPlotData && (selection==kWHSR)) continue;
    plotTree->GetEntry(ientry);
    if(weight>500.) continue;
    if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR && fj1MSD_corr<40) continue;
    if(theCategory!=kPlotData && theCategory!=kPlotQCD && (MVAVarType==2 || MVAVarType==3)) { // veto on training events
      if((selection==kWHSR || selection==kWHFJSR) &&(eventNumber % 10) < 3) continue; 
    }
    // physics calculations
    mT = TMath::Sqrt(2.*pfmet*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-pfmetphi))));
    if(selection>=kWHLightFlavorCR && selection<=kWHSR) {
      balanceVH = hbbDijetPt/topWBosonPt;  
      dPhiHbbJet1MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet1Phi - pfmetphi  ));
      dPhiHbbJet2MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet2Phi - pfmetphi  ));
      nAddJet = nJet-2;
      nAddJet_jesUp   = nJet-2; //FIX THIS
      nAddJet_jesDown = nJet-2; //FIX THIS
      mvaNSoft2=nSoft2;
      mvaNSoft5=nSoft5;
      mvaNSoft10=nSoft10;
      dPhil1W  = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - topWBosonPhi));
      dPhil1b1 = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - hbbJet1Phi  ));
      dPhil1b2 = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - hbbJet2Phi  ));
      dPhiWb1  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi - hbbJet1Phi  ));
      dPhiWb2  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi - hbbJet2Phi  ));
      dPhib1b2 = fabs(TVector2::Phi_mpi_pi(hbbJet1Phi   - hbbJet2Phi  ));
      dEtal1W  = fabs(lepton1Eta   - topWBosonEta);
      dEtal1b1 = fabs(lepton1Eta   - hbbJet1Eta  );
      dEtal1b2 = fabs(lepton1Eta   - hbbJet2Eta  );
      dEtaWb1  = fabs(topWBosonEta - hbbJet1Eta  );
      dEtaWb2  = fabs(topWBosonEta - hbbJet2Eta  );
      dEtab1b2 = fabs(hbbJet1Eta   - hbbJet2Eta  );
      dPhil1W_jesUp  = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - topWBosonPhi_jesUp));
      dPhiWb1_jesUp  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi_jesUp - hbbJet1Phi  ));
      dPhiWb2_jesUp  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi_jesUp - hbbJet2Phi  ));
      dPhil1W_jesDown  = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - topWBosonPhi_jesDown));
      dPhiWb1_jesDown  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi_jesDown - hbbJet1Phi  ));
      dPhiWb2_jesDown  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi_jesDown - hbbJet2Phi  ));
      mva_lepton1Flav = lepton1Flav;
      mva_lepton1Charge = lepton1Charge; 

      // Calculate Fatjet variables in the inclusive selection
      /*
      adjustedFatjetEta = (nFatjet==0)*(-5) + (nFatjet!=0)*(fj1Eta);
      dPhil1fj1 = (nFatjet==0)*(-1) + (nFatjet!=0)*fabs(TVector2::Phi_mpi_pi(lepton1Phi   - fj1Phi  ));
      dPhiWfj1  = (nFatjet==0)*(-1) + (nFatjet!=0)*fabs(TVector2::Phi_mpi_pi(topWBosonPhi - fj1Phi  ));
      dEtal1fj1 = (nFatjet==0)*(-1) + (nFatjet!=0)*fabs(lepton1Eta   - fj1Eta );
      mva_nIsojet = nIsojet;
      mvaPsi["021004010502"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["012004010502"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["021003011002"] = fj1ECFN["2_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["022004011002"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),4.00);
      mvaPsi["020503010502"] = fj1ECFN["2_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["022003012002"] = fj1ECFN["2_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),2.00);
      mvaPsi["012004020503"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),2.00);
      mvaPsi["021004020503"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),2.00);
      mvaPsi["032003012002"] = fj1ECFN["3_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),3.00);
      mvaPsi["012003011002"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["011003010502"] = fj1ECFN["1_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["031003011002"] = fj1ECFN["3_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),3.00);
      mvaPsi["031003012002"] = fj1ECFN["3_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.50);
      mvaPsi["030503011002"] = fj1ECFN["3_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.50);
      mvaPsi["022004012002"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),2.00);
      mvaPsi["021004011002"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["012004011002"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["022004021003"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_10"]),2.00);
      mvaPsi["012003010502"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["022003011002"] = fj1ECFN["2_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),4.00);
      mvaPsi["022004031003"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_10"]),1.33);
      mvaPsi["012004030503"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),1.33);
      mvaPsi["021004030503"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),1.33);
      mvaPsi["020504010502"] = fj1ECFN["2_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["022004030503"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),2.67);
      mvaPsi["011004010502"] = fj1ECFN["1_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["030503010502"] = fj1ECFN["3_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),3.00);
      mvaPsi["021003010502"] = fj1ECFN["2_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["012004011003"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_10"]),2.00);
      mvaPsi["012004021004"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_4_10"]),1.00);
      mvaPsi["021004012004"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_4_20"]),1.00);
      mvaPsi["011003020503"] = fj1ECFN["1_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),1.00);
      mvaPsi["020503011003"] = fj1ECFN["2_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_10"]),1.00);
      mvaPsi["012003030503"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),1.33);
      mvaPsi["030503012003"] = fj1ECFN["3_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_20"]),0.75);
      mvaPsi["010503010502"] = fj1ECFN["1_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),1.00);
      mvaPsi["020503011002"] = fj1ECFN["2_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["011003011002"] = fj1ECFN["1_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["020504011002"] = fj1ECFN["2_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["012004021003"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_10"]),1.00);
      mvaPsi["012004012002"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["020504020503"] = fj1ECFN["2_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),1.00);
      mvaPsi["011004011002"] = fj1ECFN["1_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["021004012002"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["021003012002"] = fj1ECFN["2_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["011004020503"] = fj1ECFN["1_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),1.00);
      mvaPsi["010504010502"] = fj1ECFN["1_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),1.00);
      mvaPsi["021004021003"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_10"]),1.00);
      mvaPsi["012003012002"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["022004022003"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_20"]),1.00);
      */
    } else if (selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
      dPhil1W   = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - topWBosonPhi));
      dPhil1fj1 = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - fj1Phi  ));
      dPhiWfj1  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi - fj1Phi  ));
      dEtal1fj1 = fabs(lepton1Eta   - fj1Eta );
      mva_lepton1Flav = lepton1Flav;
      mva_lepton1Charge = lepton1Charge; 
      mva_nIsojet = nIsojet;
      mvaPsi["021004010502"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["012004010502"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["021003011002"] = fj1ECFN["2_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["022004011002"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),4.00);
      mvaPsi["020503010502"] = fj1ECFN["2_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["022003012002"] = fj1ECFN["2_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),2.00);
      mvaPsi["012004020503"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),2.00);
      mvaPsi["021004020503"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),2.00);
      mvaPsi["032003012002"] = fj1ECFN["3_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),3.00);
      mvaPsi["012003011002"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["011003010502"] = fj1ECFN["1_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["031003011002"] = fj1ECFN["3_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),3.00);
      mvaPsi["031003012002"] = fj1ECFN["3_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.50);
      mvaPsi["030503011002"] = fj1ECFN["3_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.50);
      mvaPsi["022004012002"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),2.00);
      mvaPsi["021004011002"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["012004011002"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["022004021003"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_10"]),2.00);
      mvaPsi["012003010502"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["022003011002"] = fj1ECFN["2_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),4.00);
      mvaPsi["022004031003"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_10"]),1.33);
      mvaPsi["012004030503"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),1.33);
      mvaPsi["021004030503"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),1.33);
      mvaPsi["020504010502"] = fj1ECFN["2_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["022004030503"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),2.67);
      mvaPsi["011004010502"] = fj1ECFN["1_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),2.00);
      mvaPsi["030503010502"] = fj1ECFN["3_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),3.00);
      mvaPsi["021003010502"] = fj1ECFN["2_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),4.00);
      mvaPsi["012004021004"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_4_10"]),1.00);
      mvaPsi["021004012004"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_4_20"]),1.00);
      mvaPsi["011003020503"] = fj1ECFN["1_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),1.00);
      mvaPsi["020503011003"] = fj1ECFN["2_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_10"]),1.00);
      mvaPsi["012003030503"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_05"]),1.33);
      mvaPsi["030503012003"] = fj1ECFN["3_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_20"]),0.75);
      mvaPsi["010503010502"] = fj1ECFN["1_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),1.00);
      mvaPsi["020503011002"] = fj1ECFN["2_3_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["011003011002"] = fj1ECFN["1_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["020504011002"] = fj1ECFN["2_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["012004021003"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_10"]),1.00);
      mvaPsi["012004012002"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["020504020503"] = fj1ECFN["2_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),1.00);
      mvaPsi["011004011002"] = fj1ECFN["1_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),1.00);
      mvaPsi["021004012002"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["021003012002"] = fj1ECFN["2_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["011004020503"] = fj1ECFN["1_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),1.00);
      mvaPsi["010504010502"] = fj1ECFN["1_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_05"]),1.00);
      mvaPsi["021004021003"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_10"]),1.00);
      mvaPsi["012003012002"] = fj1ECFN["1_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_20"]),1.00);
      mvaPsi["022004022003"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_20"]),1.00);

      mvaPsi["012002011002"] = fj1ECFN["1_2_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_2_10"]),2.00);
      mvaPsi["014003022003"] = fj1ECFN["1_3_40"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_20"]),1.00);
      mvaPsi["031003014003"] = fj1ECFN["3_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_40"]),0.75);
      mvaPsi["031003022003"] = fj1ECFN["3_3_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_20"]),0.75);
      mvaPsi["032003034003"] = fj1ECFN["3_3_20"]/pow(TMath::Max(0.0f,fj1ECFN["3_3_40"]),0.50);
      mvaPsi["012004011003"] = fj1ECFN["1_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_10"]),2.00);
      mvaPsi["014004012003"] = fj1ECFN["1_4_40"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_20"]),2.00);
      mvaPsi["020504010503"] = fj1ECFN["2_4_05"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_05"]),2.00);
      mvaPsi["021004011003"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_10"]),2.00);
      mvaPsi["021004020503"] = fj1ECFN["2_4_10"]/pow(TMath::Max(0.0f,fj1ECFN["2_3_05"]),2.00);
      mvaPsi["022004012003"] = fj1ECFN["2_4_20"]/pow(TMath::Max(0.0f,fj1ECFN["1_3_20"]),2.00);
    }

    // Wpt correction
    float weight_wptCorrUp   = weight;
    float weight_wptCorrDown = weight;
    float weight_wptCorrNone = weight;
    if(useWptCorr && (theCategory==kPlotWbb || theCategory==kPlotWb || theCategory==kPlotWLF || theCategory==kPlotTop || theCategory==kPlotTT)) {
      float wptCorrFactor=1, wptCorrError=1;
      if(wptCorrType==0) {
        wptCorrFactor=theWptCorr->Eval(TMath::Min(topWBosonPt,(float)499.99));
        wptCorrError=1.+(theWptCorrHist->GetBinError( theWptCorrHist->FindBin(TMath::Min(topWBosonPt,(float)499.99)) ) / wptCorrFactor);
        weight *= wptCorrFactor;
        if(wptCorrError>0) weight_wptCorrUp   *= wptCorrError;
        if(wptCorrError>0) weight_wptCorrDown /= wptCorrError;
      } else if(wptCorrType==1) {
        if(theCategory==kPlotWLF) {
          weight_wptCorrUp   = (weight_wptCorrNone) * (1 - 5.29e-4 * topWBosonPt); 
          weight             = (weight_wptCorrNone) * (1 - 5.75e-4 * topWBosonPt); 
          weight_wptCorrDown = (weight_wptCorrNone) * (1 - 6.21e-4 * topWBosonPt); 
        } else if(theCategory==kPlotTT) {
          weight_wptCorrUp   = (weight_wptCorrNone) * (1 - 2.91e-4 * topWBosonPt); 
          weight             = (weight_wptCorrNone) * (1 - 3.80e-4 * topWBosonPt); 
          weight_wptCorrDown = (weight_wptCorrNone) * (1 - 4.69e-4 * topWBosonPt); 
        } else {
          weight_wptCorrUp   = (weight_wptCorrNone) * (1 - 1.54e-3 * topWBosonPt); 
          weight             = (weight_wptCorrNone) * (1 - 1.67e-3 * topWBosonPt); 
          weight_wptCorrDown = (weight_wptCorrNone) * (1 - 1.80e-3 * topWBosonPt); 
        }
        if(weight_wptCorrNone!=0) wptCorrFactor = weight/weight_wptCorrNone;
      }
      if(wptCorrFactor!=1) {
        weight_VHCorrUp    *= wptCorrFactor;
        weight_VHCorrDown  *= wptCorrFactor;
        //weight_pdfUp       *= wptCorrFactor;
        //weight_pdfDown     *= wptCorrFactor;
        weight_QCDr1f2     *= wptCorrFactor;
        weight_QCDr1f5     *= wptCorrFactor;
        weight_QCDr2f1     *= wptCorrFactor;
        weight_QCDr2f2     *= wptCorrFactor;
        weight_QCDr5f1     *= wptCorrFactor;
        weight_QCDr5f5     *= wptCorrFactor;
        weight_NLOQCDfUp   *= wptCorrFactor;
        weight_NLOQCDfDown *= wptCorrFactor;
        weight_NLOQCDrUp   *= wptCorrFactor;
        weight_NLOQCDrDown *= wptCorrFactor;
        weight_lepSFUp     *= wptCorrFactor;
        // VGluUp and VGluDown not corrected by the wptCorrFactor here because
        // the resulting shape histo gets renormalized to the nominal sum of weights later
        for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
          weight_cmvaLFUp        [iPt][iEta] *= wptCorrFactor;
          weight_cmvaHFUp        [iPt][iEta] *= wptCorrFactor;
          weight_cmvaHFStats1Up  [iPt][iEta] *= wptCorrFactor;
          weight_cmvaHFStats2Up  [iPt][iEta] *= wptCorrFactor;
          weight_cmvaLFStats1Up  [iPt][iEta] *= wptCorrFactor;
          weight_cmvaLFStats2Up  [iPt][iEta] *= wptCorrFactor;
          weight_cmvaCErr1Up     [iPt][iEta] *= wptCorrFactor;
          weight_cmvaCErr2Up     [iPt][iEta] *= wptCorrFactor;
          weight_cmvaLFDown      [iPt][iEta] *= wptCorrFactor;
          weight_cmvaHFDown      [iPt][iEta] *= wptCorrFactor;
          weight_cmvaHFStats1Down[iPt][iEta] *= wptCorrFactor;
          weight_cmvaHFStats2Down[iPt][iEta] *= wptCorrFactor;
          weight_cmvaLFStats1Down[iPt][iEta] *= wptCorrFactor;
          weight_cmvaLFStats2Down[iPt][iEta] *= wptCorrFactor;
          weight_cmvaCErr1Down   [iPt][iEta] *= wptCorrFactor;
          weight_cmvaCErr2Down   [iPt][iEta] *= wptCorrFactor;
        }
      }
    }

    float weight_VGluUp   = weight;
    float weight_VGluDown = weight;
    if(fj1HighestPtGen==21 && (
      theCategory==kPlotWbb || theCategory==kPlotWb || theCategory==kPlotWLF ||
      theCategory==kPlotZbb || theCategory==kPlotZb || theCategory==kPlotZLF)
    ) { weight_VGluUp *= 1.2; weight_VGluDown *= 0.8; }

    float theVar; 
    bool passFullSel         = (selectionBits & selection)!=0; 
    bool passFullSel_jesUp   = (selectionBits_jesUp & selection)!=0; 
    bool passFullSel_jesDown = (selectionBits_jesDown & selection)!=0; 
    bool passMassSplit=true;
    if(selection==kWHHeavyFlavorCR) { // Mass splitting for the WH HF CR
      if(lowMass  && hbbDijetMass>=90.) passMassSplit=false;
      if(!lowMass && hbbDijetMass<150.) passMassSplit=false;
    }
    float bdtValue=-99,bdtValue_jesUp=-99,bdtValue_jesDown=-99;
    switch(MVAVarType) {
      case 1:
      default:
        if(selection>=kWHLightFlavorCR && selection>=kWHPresel) {
          MVAVar=hbbDijetPt;
          MVAVar_jesUp=hbbDijetPtUp;
          MVAVar_jesDown=hbbDijetPtDown;
        } else if(selection<=kWHLightFlavorFJCR && selection<=kWHFJPresel) {
          MVAVar=fj1Pt;
          MVAVar_jesUp=fj1PtScaleUp;
          MVAVar_jesDown=fj1PtScaleDown;
        } break;
      case 2:
        if(selection==kWHLightFlavorCR || selection==kWHLightFlavorFJCR) {
          if(passFullSel        ) bdtValue        =(reader        ->EvaluateMulticlass("BDT")[1]);
          if(passFullSel_jesUp  ) bdtValue_jesUp  =(reader_jesUp  ->EvaluateMulticlass("BDT")[1]);
          if(passFullSel_jesDown) bdtValue_jesDown=(reader_jesDown->EvaluateMulticlass("BDT")[1]);
        } else if(selection==kWHHeavyFlavorCR || selection==kWHHeavyFlavorFJCR) {
          if(passFullSel        ) bdtValue        =(reader        ->EvaluateMulticlass("BDT")[2]);
          if(passFullSel_jesUp  ) bdtValue_jesUp  =(reader_jesUp  ->EvaluateMulticlass("BDT")[2]);
          if(passFullSel_jesDown) bdtValue_jesDown=(reader_jesDown->EvaluateMulticlass("BDT")[2]);
        } else if(selection==kWH2TopCR || selection==kWHTT2bFJCR || selection==kWHTT1bFJCR) {
          if(passFullSel        ) bdtValue        =(reader->EvaluateMulticlass("BDT")[3]);
          if(passFullSel_jesUp  ) bdtValue_jesUp  =(reader->EvaluateMulticlass("BDT")[3]);
          if(passFullSel_jesDown) bdtValue_jesDown=(reader->EvaluateMulticlass("BDT")[3]);
        } else {
          if(passFullSel        ) bdtValue        =(reader->EvaluateMulticlass("BDT")[0]);
          if(passFullSel_jesUp  ) bdtValue_jesUp  =(reader->EvaluateMulticlass("BDT")[0]);
          if(passFullSel_jesDown) bdtValue_jesDown=(reader->EvaluateMulticlass("BDT")[0]);
        }
        if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) {
          MVAVar=bDiscrMin; MVAVar_jesUp=bDiscrMin; MVAVar_jesDown=bDiscrMin;
        } else if(selection==kWHSR) {
          MVAVar=bdtValue; MVAVar_jesUp=bdtValue_jesUp; MVAVar_jesDown=bdtValue_jesDown;
        } else if(selection>=kWHLightFlavorFJCR && selection<=kWHTT1bFJCR) {
          MVAVar=fj1MSD_corr;
          MVAVar_jesUp=fj1MSD_corr_jesUp;
          MVAVar_jesDown=fj1MSD_corr_jesDown;
        } else if(selection==kWHFJSR) {
          MVAVar=bdtValue; MVAVar_jesUp=bdtValue_jesUp; MVAVar_jesDown=bdtValue_jesDown;
        } break;
      case 3:
        if(passFullSel        ) bdtValue        =reader        ->EvaluateMVA("BDT");
        if(passFullSel_jesUp  ) bdtValue_jesUp  =reader_jesUp  ->EvaluateMVA("BDT");
        if(passFullSel_jesDown) bdtValue_jesDown=reader_jesDown->EvaluateMVA("BDT");
        if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) {
          MVAVar=bDiscrMin; MVAVar_jesUp=bDiscrMin; MVAVar_jesDown=bDiscrMin;
        } else if(selection==kWHSR) {
          MVAVar=bdtValue; MVAVar_jesUp=bdtValue_jesUp; MVAVar_jesDown=bdtValue_jesDown;
        } else if(selection>=kWHLightFlavorFJCR && selection<=kWHTT1bFJCR) {
          MVAVar=fj1MSD_corr;
          MVAVar_jesUp=fj1MSD_corr_jesUp;
          MVAVar_jesDown=fj1MSD_corr_jesDown;
        } else if(selection==kWHFJSR) {
          MVAVar=bdtValue; MVAVar_jesUp=bdtValue_jesUp; MVAVar_jesDown=bdtValue_jesDown;
        }
        break;
    }
    MVAVar=TMath::Min(MVAVar, (float)(MVAbins[MVAbins.size()-1]-0.001));
    bool makePlot;
    if(selection>=kWHLightFlavorCR && selection<=kWHSR) { for(int p=0; p<nPlots; p++) {
      makePlot=false;
      // Variables -- change the makePlot for n-1 later
      if      (histoNames[p]=="pTH"                 ) { theVar = TMath::Min(hbbDijetPt    , float(xmax[p]- 1e-6)); makePlot = passFullSel != ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="mH"                  ) { theVar = TMath::Min(hbbDijetMass  , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="WpT"                 ) { theVar = TMath::Min(topWBosonPt   , float(xmax[p]- 1e-6)); makePlot = passFullSel != ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="deltaPhiVH"          ) { theVar = TMath::Min(deltaPhiVH    , float(xmax[p]- 1e-6)); makePlot = passFullSel != (selection==kWHSR && (nMinusOneBits & selection)!=0 && theVar<2.5); }
      else if (histoNames[p]=="pTBalanceDijetW"     ) { theVar = TMath::Min(balanceVH     , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Pt"           ) { theVar = TMath::Min(lepton1Pt     , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Flav"         ) { theVar = lepton1Flav                                     ; makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Charge"       ) { theVar = lepton1Charge                                   ; makePlot = passFullSel; }
      else if (histoNames[p]=="pfmet"               ) { theVar = TMath::Min(pfmet         , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmetphi"            ) { theVar = TMath::Min(pfmetphi      , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmetsig"            ) { theVar = TMath::Min(pfmetsig      , float(xmax[p]- 1e-6)); makePlot = passFullSel != ((selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR) && (nMinusOneBits & selection)!=0 && theVar<=2.); }
      else if (histoNames[p]=="mTW"                 ) { theVar = TMath::Min(mT            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet1Pt"            ) { theVar = TMath::Min(hbbJet1Pt     , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet2Pt"            ) { theVar = TMath::Min(hbbJet2Pt     , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet1MET"      ) { theVar = TMath::Min(dPhiHbbJet1MET, float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet2MET"      ) { theVar = TMath::Min(dPhiHbbJet2MET, float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="bDiscrMin"           ) { theVar = TMath::Min(bDiscrMin     , float(xmax[p]- 1e-6)); makePlot = passFullSel != ((nMinusOneBits & selection)!=0 && (selection==kWHSR && theVar<bDiscrLoose )); }
      else if (histoNames[p]=="bDiscrMax"           ) { theVar = TMath::Min(bDiscrMax     , float(xmax[p]- 1e-6)); makePlot = passFullSel != ((nMinusOneBits & selection)!=0 && ((selection==kWHLightFlavorCR && (theVar<bDiscrLoose || theVar>bDiscrMedium) ) || (selection!=kWHLightFlavorCR && theVar<bDiscrTight)) ); }
      else if (histoNames[p]=="topMassLep1Met"      ) { theVar = TMath::Min(topMassLep1Met, float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="nJet"                ) { theVar = TMath::Min((float)nJet   , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="sumEtSoft1"          ) { theVar = TMath::Min(sumEtSoft1    , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="nSoft2"              ) { theVar = nSoft2                                          ; makePlot = passFullSel; }
      else if (histoNames[p]=="nSoft5"              ) { theVar = nSoft5                                          ; makePlot = passFullSel; }
      else if (histoNames[p]=="nSoft10"             ) { theVar = nSoft10                                         ; makePlot = passFullSel; }
      else if (histoNames[p]=="topWBosonCosThetaCS" ) { theVar = topWBosonCosThetaCS                             ; makePlot = passFullSel; }
      else if (histoNames[p]=="hbbCosThetaJJ"       ) { theVar = hbbCosThetaJJ                                   ; makePlot = passFullSel; }
      else if (histoNames[p]=="hbbCosThetaCSJ1"     ) { theVar = hbbCosThetaCSJ1                                 ; makePlot = passFullSel; }
      else if (histoNames[p]=="dPhil1W"             ) { theVar = dPhil1W                                         ; makePlot = passFullSel; }
      else if (histoNames[p]=="dPhil1b1"            ) { theVar = dPhil1b1                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="dPhil1b2"            ) { theVar = dPhil1b2                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiWb1"             ) { theVar = dPhiWb1                                         ; makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiWb2"             ) { theVar = dPhiWb2                                         ; makePlot = passFullSel; }
      else if (histoNames[p]=="dPhib1b2"            ) { theVar = dPhib1b2                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="dEtal1W"             ) { theVar = dEtal1W                                         ; makePlot = passFullSel; }
      else if (histoNames[p]=="dEtal1b1"            ) { theVar = dEtal1b1                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="dEtal1b2"            ) { theVar = dEtal1b2                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="dEtaWb1"             ) { theVar = dEtaWb1                                         ; makePlot = passFullSel; }
      else if (histoNames[p]=="dEtaWb2"             ) { theVar = dEtaWb2                                         ; makePlot = passFullSel; }
      else if (histoNames[p]=="dEtab1b2"            ) { theVar = dEtab1b2                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="nIsojet"             ) { theVar = mva_nIsojet                                             ; makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Pt"               ) { theVar = TMath::Min(fj1Pt                 , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Eta"              ) { theVar = TMath::Min(adjustedFatjetEta     , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="maxSubjetCSV"        ) { theVar = TMath::Min(fj1MaxCSV             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="minSubjetCSV"        ) { theVar = TMath::Min(fj1MinCSV             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1DoubleCSV"        ) { theVar = TMath::Min(fj1DoubleCSV          , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1HTTMass"          ) { theVar = TMath::Min(fj1HTTMass            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1HTTFRec"          ) { theVar = TMath::Min(fj1HTTFRec            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1MSD_corr"         ) { theVar = TMath::Min(fj1MSD_corr           , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Tau32"            ) { theVar = TMath::Min(fj1Tau32              , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Tau32SD"          ) { theVar = TMath::Min(fj1Tau32SD            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Tau21"            ) { theVar = TMath::Min(fj1Tau21              , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Tau21SD"          ) { theVar = TMath::Min(fj1Tau21SD            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhil1fj1"           ) { theVar = TMath::Min(dPhil1fj1             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiWfj1"            ) { theVar = TMath::Min(dPhiWfj1              , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dEtal1fj1"           ) { theVar = TMath::Min(dEtal1fj1             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004010502"     ) { theVar = TMath::Min(mvaPsi["021004010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004010502"     ) { theVar = TMath::Min(mvaPsi["012004010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021003011002"     ) { theVar = TMath::Min(mvaPsi["021003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004011002"     ) { theVar = TMath::Min(mvaPsi["022004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020503010502"     ) { theVar = TMath::Min(mvaPsi["020503010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022003012002"     ) { theVar = TMath::Min(mvaPsi["022003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004020503"     ) { theVar = TMath::Min(mvaPsi["012004020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004020503"     ) { theVar = TMath::Min(mvaPsi["021004020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi032003012002"     ) { theVar = TMath::Min(mvaPsi["032003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003011002"     ) { theVar = TMath::Min(mvaPsi["012003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011003010502"     ) { theVar = TMath::Min(mvaPsi["011003010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi031003011002"     ) { theVar = TMath::Min(mvaPsi["031003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi031003012002"     ) { theVar = TMath::Min(mvaPsi["031003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi030503011002"     ) { theVar = TMath::Min(mvaPsi["030503011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004012002"     ) { theVar = TMath::Min(mvaPsi["022004012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004011002"     ) { theVar = TMath::Min(mvaPsi["021004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004011002"     ) { theVar = TMath::Min(mvaPsi["012004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004021003"     ) { theVar = TMath::Min(mvaPsi["022004021003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003010502"     ) { theVar = TMath::Min(mvaPsi["012003010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022003011002"     ) { theVar = TMath::Min(mvaPsi["022003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004031003"     ) { theVar = TMath::Min(mvaPsi["022004031003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004030503"     ) { theVar = TMath::Min(mvaPsi["012004030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004030503"     ) { theVar = TMath::Min(mvaPsi["021004030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020504010502"     ) { theVar = TMath::Min(mvaPsi["020504010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004030503"     ) { theVar = TMath::Min(mvaPsi["022004030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011004010502"     ) { theVar = TMath::Min(mvaPsi["011004010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi030503010502"     ) { theVar = TMath::Min(mvaPsi["030503010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021003010502"     ) { theVar = TMath::Min(mvaPsi["021003010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004011003"     ) { theVar = TMath::Min(mvaPsi["012004011003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004021004"     ) { theVar = TMath::Min(mvaPsi["012004021004"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004012004"     ) { theVar = TMath::Min(mvaPsi["021004012004"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011003020503"     ) { theVar = TMath::Min(mvaPsi["011003020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020503011003"     ) { theVar = TMath::Min(mvaPsi["020503011003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003030503"     ) { theVar = TMath::Min(mvaPsi["012003030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi030503012003"     ) { theVar = TMath::Min(mvaPsi["030503012003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi010503010502"     ) { theVar = TMath::Min(mvaPsi["010503010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020503011002"     ) { theVar = TMath::Min(mvaPsi["020503011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011003011002"     ) { theVar = TMath::Min(mvaPsi["011003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020504011002"     ) { theVar = TMath::Min(mvaPsi["020504011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004021003"     ) { theVar = TMath::Min(mvaPsi["012004021003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004012002"     ) { theVar = TMath::Min(mvaPsi["012004012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020504020503"     ) { theVar = TMath::Min(mvaPsi["020504020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011004011002"     ) { theVar = TMath::Min(mvaPsi["011004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004012002"     ) { theVar = TMath::Min(mvaPsi["021004012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021003012002"     ) { theVar = TMath::Min(mvaPsi["021003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011004020503"     ) { theVar = TMath::Min(mvaPsi["011004020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi010504010502"     ) { theVar = TMath::Min(mvaPsi["010504010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004021003"     ) { theVar = TMath::Min(mvaPsi["021004021003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003012002"     ) { theVar = TMath::Min(mvaPsi["012003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004022003"     ) { theVar = TMath::Min(mvaPsi["022004022003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="bdtValue"            ) { theVar = bdtValue                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="MVAVar"              ) { theVar = MVAVar                                          ; makePlot = passFullSel && passMassSplit; }
      else continue;
      if(!makePlot) continue;
      histos[p][theCategory]->Fill(theVar, weight);
    }} else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJPresel) { for(int p=0; p<nPlots; p++) {
      if      (histoNames[p]=="lepton1Pt"           ) { theVar = TMath::Min(lepton1Pt             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Flav"         ) { theVar = lepton1Flav                                             ; makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Charge"       ) { theVar = lepton1Charge                                           ; makePlot = passFullSel; }
      else if (histoNames[p]=="deltaPhiVH"          ) { theVar = TMath::Min(deltaPhiVH            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="deltaPhiLep1Met"     ) { theVar = TMath::Min(deltaPhiLep1Met       , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="topWBosonPt"         ) { theVar = TMath::Min(topWBosonPt           , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="mT"                  ) { theVar = TMath::Min(mT                    , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmet"               ) { theVar = TMath::Min(pfmet         , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmetphi"            ) { theVar = TMath::Min(pfmetphi      , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmetsig"            ) { theVar = TMath::Min(pfmetsig      , float(xmax[p]- 1e-6)); makePlot = passFullSel != ((selection==kWHLightFlavorFJCR || selection==kWHHeavyFlavorFJCR) && (nMinusOneBits & selection)!=0 && theVar<=2.); }
      else if (histoNames[p]=="dPhil1W"             ) { theVar = TMath::Min(dPhil1W               , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="nIsojet"             ) { theVar = mva_nIsojet                                             ; makePlot = passFullSel; }
      else if (histoNames[p]=="isojetNBtags"        ) { theVar = isojetNBtags                                            ; makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Pt"               ) { theVar = TMath::Min(fj1Pt                 , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Eta"              ) { theVar = TMath::Min(fj1Eta                , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="maxSubjetCSV"        ) { theVar = TMath::Min(fj1MaxCSV             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="minSubjetCSV"        ) { theVar = TMath::Min(fj1MinCSV             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1DoubleCSV"        ) { theVar = TMath::Min(fj1DoubleCSV          , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1HTTMass"          ) { theVar = TMath::Min(fj1HTTMass            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1HTTFRec"          ) { theVar = TMath::Min(fj1HTTFRec            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1MSD_corr"         ) { theVar = TMath::Min(fj1MSD_corr           , float(xmax[p]- 1e-6)); makePlot = passFullSel != ((nMinusOneBits&selection)!=0 && theVar<40); }
      else if (histoNames[p]=="fj1Tau32"            ) { theVar = TMath::Min(fj1Tau32              , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Tau32SD"          ) { theVar = TMath::Min(fj1Tau32SD            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Tau21"            ) { theVar = TMath::Min(fj1Tau21              , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="fj1Tau21SD"          ) { theVar = TMath::Min(fj1Tau21SD            , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhil1fj1"           ) { theVar = TMath::Min(dPhil1fj1             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiWfj1"            ) { theVar = TMath::Min(dPhiWfj1              , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="dEtal1fj1"           ) { theVar = TMath::Min(dEtal1fj1             , float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004010502"     ) { theVar = TMath::Min(mvaPsi["021004010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004010502"     ) { theVar = TMath::Min(mvaPsi["012004010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021003011002"     ) { theVar = TMath::Min(mvaPsi["021003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004011002"     ) { theVar = TMath::Min(mvaPsi["022004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020503010502"     ) { theVar = TMath::Min(mvaPsi["020503010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022003012002"     ) { theVar = TMath::Min(mvaPsi["022003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004020503"     ) { theVar = TMath::Min(mvaPsi["012004020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi032003012002"     ) { theVar = TMath::Min(mvaPsi["032003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003011002"     ) { theVar = TMath::Min(mvaPsi["012003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011003010502"     ) { theVar = TMath::Min(mvaPsi["011003010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi031003011002"     ) { theVar = TMath::Min(mvaPsi["031003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi031003012002"     ) { theVar = TMath::Min(mvaPsi["031003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi030503011002"     ) { theVar = TMath::Min(mvaPsi["030503011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004012002"     ) { theVar = TMath::Min(mvaPsi["022004012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004011002"     ) { theVar = TMath::Min(mvaPsi["021004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004011002"     ) { theVar = TMath::Min(mvaPsi["012004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004021003"     ) { theVar = TMath::Min(mvaPsi["022004021003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003010502"     ) { theVar = TMath::Min(mvaPsi["012003010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022003011002"     ) { theVar = TMath::Min(mvaPsi["022003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004031003"     ) { theVar = TMath::Min(mvaPsi["022004031003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004030503"     ) { theVar = TMath::Min(mvaPsi["012004030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004030503"     ) { theVar = TMath::Min(mvaPsi["021004030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020504010502"     ) { theVar = TMath::Min(mvaPsi["020504010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004030503"     ) { theVar = TMath::Min(mvaPsi["022004030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011004010502"     ) { theVar = TMath::Min(mvaPsi["011004010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi030503010502"     ) { theVar = TMath::Min(mvaPsi["030503010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021003010502"     ) { theVar = TMath::Min(mvaPsi["021003010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004021004"     ) { theVar = TMath::Min(mvaPsi["012004021004"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004012004"     ) { theVar = TMath::Min(mvaPsi["021004012004"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011003020503"     ) { theVar = TMath::Min(mvaPsi["011003020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020503011003"     ) { theVar = TMath::Min(mvaPsi["020503011003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003030503"     ) { theVar = TMath::Min(mvaPsi["012003030503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi030503012003"     ) { theVar = TMath::Min(mvaPsi["030503012003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi010503010502"     ) { theVar = TMath::Min(mvaPsi["010503010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020503011002"     ) { theVar = TMath::Min(mvaPsi["020503011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011003011002"     ) { theVar = TMath::Min(mvaPsi["011003011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020504011002"     ) { theVar = TMath::Min(mvaPsi["020504011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004021003"     ) { theVar = TMath::Min(mvaPsi["012004021003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004012002"     ) { theVar = TMath::Min(mvaPsi["012004012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020504020503"     ) { theVar = TMath::Min(mvaPsi["020504020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011004011002"     ) { theVar = TMath::Min(mvaPsi["011004011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004012002"     ) { theVar = TMath::Min(mvaPsi["021004012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021003012002"     ) { theVar = TMath::Min(mvaPsi["021003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi011004020503"     ) { theVar = TMath::Min(mvaPsi["011004020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi010504010502"     ) { theVar = TMath::Min(mvaPsi["010504010502"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004021003"     ) { theVar = TMath::Min(mvaPsi["021004021003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012003012002"     ) { theVar = TMath::Min(mvaPsi["012003012002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004022003"     ) { theVar = TMath::Min(mvaPsi["022004022003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012002011002"     ) { theVar = TMath::Min(mvaPsi["012002011002"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi014003022003"     ) { theVar = TMath::Min(mvaPsi["014003022003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi031003014003"     ) { theVar = TMath::Min(mvaPsi["031003014003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi031003022003"     ) { theVar = TMath::Min(mvaPsi["031003022003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi032003034003"     ) { theVar = TMath::Min(mvaPsi["032003034003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi012004011003"     ) { theVar = TMath::Min(mvaPsi["012004011003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi014004012003"     ) { theVar = TMath::Min(mvaPsi["014004012003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi020504010503"     ) { theVar = TMath::Min(mvaPsi["020504010503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004011003"     ) { theVar = TMath::Min(mvaPsi["021004011003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi021004020503"     ) { theVar = TMath::Min(mvaPsi["021004020503"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="psi022004012003"     ) { theVar = TMath::Min(mvaPsi["022004012003"], float(xmax[p]- 1e-6)); makePlot = passFullSel; }
      else if (histoNames[p]=="bdtValue"            ) { theVar = bdtValue                                        ; makePlot = passFullSel; }
      else if (histoNames[p]=="MVAVar"              ) { theVar = MVAVar                                          ; makePlot = passFullSel; }
      else continue;
      if(!makePlot) continue;
      histos[p][theCategory]->Fill(theVar, weight);
    }}
    // Fill systematic shapes
    if(passFullSel && passMassSplit && theCategory!=kPlotData) {
      histo_VHCorrUp    [theCategory]->Fill(MVAVar, weight_VHCorrUp    ); 
      histo_VHCorrDown  [theCategory]->Fill(MVAVar, weight_VHCorrDown  ); 
      //histo_pdfUp       [theCategory]->Fill(MVAVar, weight_pdfUp    ); 
      //histo_pdfDown     [theCategory]->Fill(MVAVar, weight_pdfDown  ); 
      histo_QCDr1f2     [theCategory]->Fill(MVAVar, weight_QCDr1f2  ); 
      histo_QCDr1f5     [theCategory]->Fill(MVAVar, weight_QCDr1f5  ); 
      histo_QCDr2f1     [theCategory]->Fill(MVAVar, weight_QCDr2f1  ); 
      histo_QCDr2f2     [theCategory]->Fill(MVAVar, weight_QCDr2f2  ); 
      histo_QCDr5f1     [theCategory]->Fill(MVAVar, weight_QCDr5f1  ); 
      histo_QCDr5f5     [theCategory]->Fill(MVAVar, weight_QCDr5f5  ); 
      histo_NLOQCDfUp   [theCategory]->Fill(MVAVar, weight_NLOQCDfUp  ); 
      histo_NLOQCDfDown [theCategory]->Fill(MVAVar, weight_NLOQCDfDown); 
      histo_NLOQCDrUp   [theCategory]->Fill(MVAVar, weight_NLOQCDrUp  ); 
      histo_NLOQCDrDown [theCategory]->Fill(MVAVar, weight_NLOQCDrDown); 
      histo_muSFUp      [theCategory]->Fill(MVAVar, (typeLepSel==1&&weight_lepSFUp>0?weight_lepSFUp : weight)); 
      histo_muSFDown    [theCategory]->Fill(MVAVar, (typeLepSel==1&&weight_lepSFUp>0?weight*weight/weight_lepSFUp : weight)); 
      histo_eleSFUp     [theCategory]->Fill(MVAVar, (typeLepSel==2&&weight_lepSFUp>0?weight_lepSFUp : weight)); 
      histo_eleSFDown   [theCategory]->Fill(MVAVar, (typeLepSel==2&&weight_lepSFUp>0?weight*weight/weight_lepSFUp : weight)); 
      histo_VGluUp      [theCategory]->Fill(MVAVar, weight_VGluUp     );
      histo_VGluDown    [theCategory]->Fill(MVAVar, weight_VGluDown   );
      histo_wptCorrUp   [theCategory]->Fill(MVAVar, weight_wptCorrUp  );
      histo_wptCorrDown [theCategory]->Fill(MVAVar, weight_wptCorrDown);
      histo_wptCorrNone [theCategory]->Fill(MVAVar, weight_wptCorrNone);
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        histo_cmvaJESUp       [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFUp        [iPt][iEta]);
        histo_cmvaLFUp        [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFUp        [iPt][iEta]);
        histo_cmvaHFUp        [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaHFUp        [iPt][iEta]);
        histo_cmvaHFStats1Up  [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaHFStats1Up  [iPt][iEta]);
        histo_cmvaHFStats2Up  [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaHFStats2Up  [iPt][iEta]);
        histo_cmvaLFStats1Up  [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFStats1Up  [iPt][iEta]);
        histo_cmvaLFStats2Up  [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFStats2Up  [iPt][iEta]);
        histo_cmvaCErr1Up     [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaCErr1Up     [iPt][iEta]);
        histo_cmvaCErr2Up     [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaCErr2Up     [iPt][iEta]);
        histo_cmvaJESDown     [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFDown      [iPt][iEta]);
        histo_cmvaLFDown      [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFDown      [iPt][iEta]);
        histo_cmvaHFDown      [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaHFDown      [iPt][iEta]);
        histo_cmvaHFStats1Down[iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaHFStats1Down[iPt][iEta]);
        histo_cmvaHFStats2Down[iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaHFStats2Down[iPt][iEta]);
        histo_cmvaLFStats1Down[iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFStats1Down[iPt][iEta]);
        histo_cmvaLFStats2Down[iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaLFStats2Down[iPt][iEta]);
        histo_cmvaCErr1Down   [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaCErr1Down   [iPt][iEta]);
        histo_cmvaCErr2Down   [iPt][iEta][theCategory]->Fill(MVAVar, weight_cmvaCErr2Down   [iPt][iEta]);
      }
    }
    if(passFullSel_jesUp && theCategory!=kPlotData)
      histo_jesUp   [theCategory]->Fill(MVAVar_jesUp, weight); 
    if(passFullSel_jesDown && theCategory!=kPlotData)
      histo_jesDown [theCategory]->Fill(MVAVar_jesDown, weight);
  }
  if(sf_training!=1) for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
    if(ic==kPlotQCD) continue;
    for(int p=0; p<=pMVAVar; p++) histos[p][ic]->Scale(sf_training);
    histo_VHCorrUp    [ic]->Scale(sf_training);
    histo_VHCorrDown  [ic]->Scale(sf_training);
    //histo_pdfUp       [ic]->Scale(sf_training);
    //histo_pdfDown     [ic]->Scale(sf_training);
    histo_QCDr1f2     [ic]->Scale(sf_training);
    histo_QCDr1f5     [ic]->Scale(sf_training);
    histo_QCDr2f1     [ic]->Scale(sf_training);
    histo_QCDr2f2     [ic]->Scale(sf_training);
    histo_QCDr5f1     [ic]->Scale(sf_training);
    histo_QCDr5f5     [ic]->Scale(sf_training);
    histo_NLOQCDfUp   [ic]->Scale(sf_training);
    histo_NLOQCDfDown [ic]->Scale(sf_training);
    histo_NLOQCDrUp   [ic]->Scale(sf_training);
    histo_NLOQCDrDown [ic]->Scale(sf_training);
    histo_muSFUp      [ic]->Scale(sf_training);
    histo_muSFDown    [ic]->Scale(sf_training);
    histo_eleSFUp     [ic]->Scale(sf_training);
    histo_eleSFDown   [ic]->Scale(sf_training);
    histo_VGluUp      [ic]->Scale(sf_training);
    histo_VGluDown    [ic]->Scale(sf_training);
    histo_wptCorrUp   [ic]->Scale(sf_training);
    histo_wptCorrDown [ic]->Scale(sf_training);
    histo_wptCorrNone [ic]->Scale(sf_training);
    histo_jesUp       [ic]->Scale(sf_training);
    histo_jesDown     [ic]->Scale(sf_training);
    for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
      histo_cmvaJESUp       [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaLFUp        [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaHFUp        [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaHFStats1Up  [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaHFStats2Up  [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaLFStats1Up  [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaLFStats2Up  [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaCErr1Up     [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaCErr2Up     [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaJESDown     [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaLFDown      [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaHFDown      [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaHFStats1Down[iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaHFStats2Down[iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaLFStats1Down[iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaLFStats2Down[iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaCErr1Down   [iPt][iEta][ic]->Scale(sf_training);
      histo_cmvaCErr2Down   [iPt][iEta][ic]->Scale(sf_training);
    }

  }
  // Renormalize the W+glu and WPT shape uncertainties to the nominal shape
  printf("before Wpt renorm:\n");
  printf("  histos   [pMVAVar][kPlotWLF]->Integral() = %.2e\n", histos   [pMVAVar][kPlotWLF]->Integral());
  printf("  histo_VHCorrUp    [kPlotWLF]->Integral() = %.2e\n", histo_VHCorrUp    [kPlotWLF]->Integral());
  printf("  histo_VHCorrDown  [kPlotWLF]->Integral() = %.2e\n", histo_VHCorrDown  [kPlotWLF]->Integral());
  //printf("  histo_pdfUp       [kPlotWLF]->Integral() = %.2e\n", histo_pdfUp       [kPlotWLF]->Integral());
  //printf("  histo_pdfDown     [kPlotWLF]->Integral() = %.2e\n", histo_pdfDown     [kPlotWLF]->Integral());
  printf("  histo_QCDr1f2     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr1f2     [kPlotWLF]->Integral());
  printf("  histo_QCDr1f5     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr1f5     [kPlotWLF]->Integral());
  printf("  histo_QCDr2f1     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr2f1     [kPlotWLF]->Integral());
  printf("  histo_QCDr2f2     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr2f2     [kPlotWLF]->Integral());
  printf("  histo_QCDr5f1     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr5f1     [kPlotWLF]->Integral());
  printf("  histo_QCDr5f5     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr5f5     [kPlotWLF]->Integral());
  if(useWptCorr) for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
    if(ic==kPlotWbb || ic==kPlotWb || ic==kPlotWLF || ic==kPlotTop || ic==kPlotTT) {
      histo_wptCorrUp   [ic]->Scale(histo_wptCorrNone[ic]->Integral()/histo_wptCorrUp   [ic]->Integral());
      histo_wptCorrDown [ic]->Scale(histo_wptCorrNone[ic]->Integral()/histo_wptCorrDown [ic]->Integral());
      float maintainNormAfterWptCorr = histo_wptCorrNone[ic]->Integral()/histos[pMVAVar][ic]->Integral();
      histos   [pMVAVar][ic]->Scale(maintainNormAfterWptCorr);
      histo_VHCorrUp    [ic]->Scale(maintainNormAfterWptCorr);
      histo_VHCorrDown  [ic]->Scale(maintainNormAfterWptCorr);
      //histo_pdfUp       [ic]->Scale(maintainNormAfterWptCorr);
      //histo_pdfDown     [ic]->Scale(maintainNormAfterWptCorr);
      histo_QCDr1f2     [ic]->Scale(maintainNormAfterWptCorr);
      histo_QCDr1f5     [ic]->Scale(maintainNormAfterWptCorr);
      histo_QCDr2f1     [ic]->Scale(maintainNormAfterWptCorr);
      histo_QCDr2f2     [ic]->Scale(maintainNormAfterWptCorr);
      histo_QCDr5f1     [ic]->Scale(maintainNormAfterWptCorr);
      histo_QCDr5f5     [ic]->Scale(maintainNormAfterWptCorr);
      histo_NLOQCDfUp   [ic]->Scale(maintainNormAfterWptCorr);
      histo_NLOQCDfDown [ic]->Scale(maintainNormAfterWptCorr);
      histo_NLOQCDrUp   [ic]->Scale(maintainNormAfterWptCorr);
      histo_NLOQCDrDown [ic]->Scale(maintainNormAfterWptCorr);
      histo_muSFUp      [ic]->Scale(maintainNormAfterWptCorr);
      histo_muSFDown    [ic]->Scale(maintainNormAfterWptCorr);
      histo_eleSFUp     [ic]->Scale(maintainNormAfterWptCorr);
      histo_eleSFDown   [ic]->Scale(maintainNormAfterWptCorr);
      histo_VGluUp      [ic]->Scale(maintainNormAfterWptCorr);
      histo_VGluDown    [ic]->Scale(maintainNormAfterWptCorr);
      histo_jesUp       [ic]->Scale(maintainNormAfterWptCorr);
      histo_jesDown     [ic]->Scale(maintainNormAfterWptCorr);
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        histo_cmvaJESUp       [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaLFUp        [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaHFUp        [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaHFStats1Up  [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaHFStats2Up  [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaLFStats1Up  [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaLFStats2Up  [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaCErr1Up     [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaCErr2Up     [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaJESDown     [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaLFDown      [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaHFDown      [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaHFStats1Down[iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaHFStats2Down[iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaLFStats1Down[iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaLFStats2Down[iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaCErr1Down   [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
        histo_cmvaCErr2Down   [iPt][iEta][ic]->Scale(maintainNormAfterWptCorr);
      }
    }
  }
  printf("after Wpt renorm:\n");
  printf("  histos   [pMVAVar][kPlotWLF]->Integral() = %.2e\n", histos   [pMVAVar][kPlotWLF]->Integral());
  printf("  histo_VHCorrUp    [kPlotWLF]->Integral() = %.2e\n", histo_VHCorrUp    [kPlotWLF]->Integral());
  printf("  histo_VHCorrDown  [kPlotWLF]->Integral() = %.2e\n", histo_VHCorrDown  [kPlotWLF]->Integral());
  //printf("  histo_pdfUp       [kPlotWLF]->Integral() = %.2e\n", histo_pdfUp       [kPlotWLF]->Integral());
  //printf("  histo_pdfDown     [kPlotWLF]->Integral() = %.2e\n", histo_pdfDown     [kPlotWLF]->Integral());
  printf("  histo_QCDr1f2     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr1f2     [kPlotWLF]->Integral());
  printf("  histo_QCDr1f5     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr1f5     [kPlotWLF]->Integral());
  printf("  histo_QCDr2f1     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr2f1     [kPlotWLF]->Integral());
  printf("  histo_QCDr2f2     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr2f2     [kPlotWLF]->Integral());
  printf("  histo_QCDr5f1     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr5f1     [kPlotWLF]->Integral());
  printf("  histo_QCDr5f5     [kPlotWLF]->Integral() = %.2e\n", histo_QCDr5f5     [kPlotWLF]->Integral());
  
  // All normalizations done, start writing output


  char regionName[128];
  TString massSuffix="";
  if(selection==kWHHeavyFlavorCR) { // Mass splitting for the WH HF CR
    if(lowMass) massSuffix="LowMass";
    else        massSuffix="HighMass";
  }
  char leptonChar='l';
  if(theLepSel==1) leptonChar='m';
  if(theLepSel==2) leptonChar='e';
  sprintf(regionName, "W%cn%s%s", leptonChar,selectionNames[selection].Data(),massSuffix.Data());
  char outputFileName[512];
  sprintf(outputFileName,"MitVHBBAnalysis/datacards/%s/plots_%s.root",dataCardDir.Data(),regionName);
  TFile *output_plots = new TFile(outputFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  for(int p=0; p<nPlots; p++) {
    if(histoNames[p]=="") continue;
    TDirectory *plotDir = output_plots->mkdir(histoNames[p]);
    plotDir->cd();
    for(theCategory=kPlotData; theCategory!=nPlotCategories; theCategory++)
      histos[p][theCategory]->Write();
  }
  output_plots->Close();
  // Hack to symmetrize the JES uncertainties
  for(int ic=kPlotData+1; ic<nPlotCategories; ic++) for(int nb=1; nb<=histos[pMVAVar][kPlotData]->GetNbinsX(); nb++) {
    if(histo_jesUp[ic]->GetBinContent(nb)>0) 
      histo_jesDown[ic]->SetBinContent(nb,
        pow(histos[pMVAVar][ic]->GetBinContent(nb),2) / histo_jesUp[ic]->GetBinContent(nb)
      );
  }

  // Compute QCD Scale envelope and bound the systematic shapes to [0.5,2] the nominal
  for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
    for(int nb=1; nb<=histos[pMVAVar][kPlotData]->GetNbinsX(); nb++) {
      double nomYield=histos[pMVAVar][ic]->GetBinContent(nb);
      double halfNom=nomYield/2.,twiceNom=nomYield*2.;
      histo_QCDrScaleUp  [ic]->SetBinContent(nb, nomYield);
      histo_QCDrScaleDown[ic]->SetBinContent(nb, nomYield);
      histo_QCDfScaleUp  [ic]->SetBinContent(nb, nomYield);
      histo_QCDfScaleDown[ic]->SetBinContent(nb, nomYield);
      // Symmetrize this shape uncertainty
      if(fabs(histo_QCDr1f2[ic]->GetBinContent(nb)) > histo_QCDfScaleUp  [ic]->GetBinContent(nb)) histo_QCDfScaleUp  [ic]->SetBinContent(nb, fabs(histo_QCDr1f2[ic]->GetBinContent(nb)));
      if(fabs(histo_QCDr1f5[ic]->GetBinContent(nb)) > histo_QCDfScaleUp  [ic]->GetBinContent(nb)) histo_QCDfScaleUp  [ic]->SetBinContent(nb, fabs(histo_QCDr1f5[ic]->GetBinContent(nb)));
      if(fabs(histo_QCDr2f1[ic]->GetBinContent(nb)) > histo_QCDrScaleUp  [ic]->GetBinContent(nb)) histo_QCDrScaleUp  [ic]->SetBinContent(nb, fabs(histo_QCDr2f1[ic]->GetBinContent(nb)));
      if(fabs(histo_QCDr5f1[ic]->GetBinContent(nb)) > histo_QCDrScaleUp  [ic]->GetBinContent(nb)) histo_QCDrScaleUp  [ic]->SetBinContent(nb, fabs(histo_QCDr5f1[ic]->GetBinContent(nb)));
      if(histo_QCDfScaleUp[ic]->GetBinContent(nb)>0) histo_QCDfScaleDown[ic]->SetBinContent(nb, nomYield*nomYield/histo_QCDfScaleUp[ic]->GetBinContent(nb));
      if(histo_QCDrScaleUp[ic]->GetBinContent(nb)>0) histo_QCDrScaleDown[ic]->SetBinContent(nb, nomYield*nomYield/histo_QCDrScaleUp[ic]->GetBinContent(nb));
      if(ic==kPlotWLF) 
        printf("bin %d: nomYield %.2e, histo_QCDr2f1 %.2e, histo_QCDr5f1 %.2e, histo_QCDrScaleUp %.2e, histo_QCDrScaleDown %.2e\n",
          nb,
          nomYield,
          histo_QCDr2f1[ic]->GetBinContent(nb),
          histo_QCDr5f1[ic]->GetBinContent(nb),
          histo_QCDrScaleUp  [ic]->GetBinContent(nb),
          histo_QCDrScaleUp  [ic]->GetBinContent(nb)
        );
      // Bound between [0.5,2.0] (essentially does nothing)
      histo_QCDrScaleUp  [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_QCDrScaleUp  [ic]->GetBinContent(nb))));
      histo_QCDrScaleDown[ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_QCDrScaleDown[ic]->GetBinContent(nb))));
      histo_QCDfScaleUp  [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_QCDfScaleUp  [ic]->GetBinContent(nb))));
      histo_QCDfScaleDown[ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_QCDfScaleDown[ic]->GetBinContent(nb))));
      //histo_pdfUp       [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_pdfUp       [ic]->GetBinContent(nb))));
      //histo_pdfDown     [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_pdfDown     [ic]->GetBinContent(nb))));
      histo_VHCorrUp    [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_VHCorrUp    [ic]->GetBinContent(nb))));
      histo_VHCorrDown  [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_VHCorrDown  [ic]->GetBinContent(nb))));
      histo_muSFUp      [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_muSFUp      [ic]->GetBinContent(nb))));
      histo_muSFDown    [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_muSFDown    [ic]->GetBinContent(nb))));
      histo_eleSFUp     [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_eleSFUp     [ic]->GetBinContent(nb))));
      histo_eleSFDown   [ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_eleSFDown   [ic]->GetBinContent(nb))));
      for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
        histo_cmvaJESUp       [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaJESUp       [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaLFUp        [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaLFUp        [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaHFUp        [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaHFUp        [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaHFStats1Up  [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaHFStats1Up  [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaHFStats2Up  [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaHFStats2Up  [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaLFStats1Up  [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaLFStats1Up  [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaLFStats2Up  [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaLFStats2Up  [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaCErr1Up     [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaCErr1Up     [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaCErr2Up     [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaCErr2Up     [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaJESDown     [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaJESDown     [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaLFDown      [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaLFDown      [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaHFDown      [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaHFDown      [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaHFStats1Down[iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaHFStats1Down[iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaHFStats2Down[iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaHFStats2Down[iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaLFStats1Down[iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaLFStats1Down[iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaLFStats2Down[iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaLFStats2Down[iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaCErr1Down   [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaCErr1Down   [iPt][iEta][ic]->GetBinContent(nb))));
        histo_cmvaCErr2Down   [iPt][iEta][ic]->SetBinContent(nb, TMath::Max(halfNom,TMath::Min(twiceNom,histo_cmvaCErr2Down   [iPt][iEta][ic]->GetBinContent(nb))));
      }
    }
  }
  // Renormalize the W+glu shape uncertainties to the nominal shape
  for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
    if(ic==kPlotWbb || ic==kPlotWb || ic==kPlotWLF ||
       ic==kPlotZbb || ic==kPlotZb || ic==kPlotZLF) {
      histo_VGluUp      [ic]->Scale(histos[pMVAVar][ic]->Integral()/histo_VGluUp  [ic]->Integral());
      histo_VGluDown    [ic]->Scale(histos[pMVAVar][ic]->Integral()/histo_VGluDown[ic]->Integral());
    }
  }

  // Prepare the datacard C-('.' Q)
  int ic0 = kPlotQCD+1;
  const char *ttbar = plotBaseNames[kPlotTT ].Data(), *Wbb=plotBaseNames[kPlotWbb].Data(), *Wb=plotBaseNames[kPlotWb ].Data(), *WLF=plotBaseNames[kPlotWLF].Data();
  
  char outputLimitsShape[512];
  sprintf(outputLimitsShape,"MitVHBBAnalysis/datacards/%s/datacard_%s.txt",dataCardDir.Data(),regionName);
  ofstream card; card.open(outputLimitsShape);
  TFile *shapesFile = TFile::Open(Form("MitVHBBAnalysis/datacards/%s/hists_%s.root",dataCardDir.Data(),regionName),"recreate");
  card << Form("imax 1 number of channels\n");
  card << Form("jmax * number of background\n");
  card << Form("kmax * number of nuisance parameters\n");
  card << Form("shapes * %s hists_%s.root %s_%s_$PROCESS %s_%s_$PROCESS_$SYSTEMATIC",regionName,regionName,shapeType,regionName,shapeType,regionName)<<std::endl;
  card << Form("bin     "); card<<regionName<<std::endl;
  card << Form("Observation %d\n", (int)(isBlinded? 0 : histos[pMVAVar][kPlotData]->Integral())); 
  card << Form("bin     "); for(int ic=ic0; ic<nPlotCategories; ic++) card<<regionName<<" "; card<<std::endl;
  card << Form("process "); for(int ic=ic0; ic<nPlotCategories; ic++) card<<Form("%9s ",plotBaseNames[ic].Data()); card<<std::endl;
  card << Form("process "); for(int ic=ic0; ic<nPlotCategories; ic++) card<<Form("%9d ", ic==kPlotVH?0:ic); card<<std::endl;
  card << Form("rate    "); for(int ic=ic0; ic<nPlotCategories; ic++) card<<Form("%9.3f ", histos[pMVAVar][ic]->Integral()); card<<std::endl;
  // CHECKPOINT

  // Nominal Shape
  histos[pMVAVar][kPlotData]->Write(Form("%s_%s_data_obs",shapeType,regionName));
  vector<TString> shapeName(nPlotCategories);
  for(int ic=ic0; ic<nPlotCategories; ic++) {
    shapeName[ic]=Form("%s_%s_%s",shapeType,regionName,plotBaseNames[ic].Data());
    histos[pMVAVar][ic]->Write(shapeName[ic]);
  }
  // Systematics - Shape Uncertainties
  TString systName;
  
  // Stat Bounding
  TH1F *statBoundingUp=0, *statBoundingDown=0;
  if(selection==kWHSR || selection==kWHFJSR) for(int ic=ic0; ic<nPlotCategories; ic++) {
    for(int nb=1; nb<=histos[pMVAVar][ic]->GetNbinsX(); nb++) {
      //singleClassBDTShape_WenWHSR_$PROCESS_$SYSTEMATIC
      systName=Form("%sStatBounding_%s_bin%d_13TeV",plotBaseNames[ic].Data(),regionName,nb);
      statBoundingUp   = (TH1F*)histos[pMVAVar][ic]->Clone(Form("%s_%s_%s_%sUp"  ,shapeType,regionName,plotBaseNames[ic].Data(),systName.Data()));
      statBoundingDown = (TH1F*)histos[pMVAVar][ic]->Clone(Form("%s_%s_%s_%sDown",shapeType,regionName,plotBaseNames[ic].Data(),systName.Data()));
      if(histos[pMVAVar][ic]->GetBinContent(nb)>0) {
        statBoundingUp  ->SetBinContent(nb, TMath::Max(1e-6,histos[pMVAVar][ic]->GetBinContent(nb) + histos[pMVAVar][ic]->GetBinError(nb)));
        statBoundingDown->SetBinContent(nb, TMath::Max(1e-6,histos[pMVAVar][ic]->GetBinContent(nb) - histos[pMVAVar][ic]->GetBinError(nb)));
      } else {
        statBoundingUp  ->SetBinContent(nb, 1e-6f); 
        statBoundingDown->SetBinContent(nb, 1e-6f); 
      }
      card<<(systName+" shape ").Data(); for(int ic2=ic0; ic2<nPlotCategories; ic2++) { 
        if(ic2==ic) { 
          card<<"1.0 ";
          statBoundingUp  ->Write(statBoundingUp  ->GetName());
          statBoundingDown->Write(statBoundingDown->GetName());
        } else card<<"- ";
      }
      card<<std::endl;
    }
  }

  // VH EWK Corrections
  card<<Form("VH_EWKCorr shape "); for(int ic=ic0; ic<nPlotCategories; ic++) {
    if(ic==kPlotVH) { card<<"1.0 "; histo_VHCorrUp[ic]->Write(shapeName[ic]+"_VH_EWKCorrUp"); histo_VHCorrDown[ic]->Write(shapeName[ic]+"_VH_EWKCorrDown"); } 
    else card<<"- ";
  } card<<std::endl;
  // QCD Scale
  for(int ic=ic0; ic<nPlotCategories; ic++) {
    if(ic==kPlotTop) continue; // Bug in MiniAOD
    card<<Form("QCDrScale%s shape ",plotBaseNames[ic].Data());
    for(int jc=ic0; jc<nPlotCategories; jc++) {
      if(jc==ic) {
        card<<"1.0 ";
        histo_QCDrScaleUp  [ic]->Write(shapeName[ic]+Form("_QCDrScale%sUp"  ,plotBaseNames[ic].Data()));
        histo_QCDrScaleDown[ic]->Write(shapeName[ic]+Form("_QCDrScale%sDown",plotBaseNames[ic].Data()));
      } else card<<"- ";
    }
    card<<std::endl;
    card<<Form("QCDfScale%s shape ",plotBaseNames[ic].Data());
    for(int jc=ic0; jc<nPlotCategories; jc++) {
      if(jc==ic) {
        card<<"1.0 ";
        histo_QCDfScaleUp  [ic]->Write(shapeName[ic]+Form("_QCDfScale%sUp"  ,plotBaseNames[ic].Data()));
        histo_QCDfScaleDown[ic]->Write(shapeName[ic]+Form("_QCDfScale%sDown",plotBaseNames[ic].Data()));
      } else card<<"- ";
    }
    card<<std::endl;
  }
  // BTag Shape Uncertainty
  for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
    systName=Form("CMS_bTagWeightJES_Pt%d_Eta%d"     , iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaJESUp     [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaJESDown     [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightLF_Pt%d_Eta%d"      , iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaLFUp      [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaLFDown      [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightHF_Pt%d_Eta%d"      , iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaHFUp      [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaHFDown      [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightHFStats1_Pt%d_Eta%d", iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaHFStats1Up[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaHFStats1Down[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightHFStats2_Pt%d_Eta%d", iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaHFStats2Up[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaHFStats2Down[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightLFStats1_Pt%d_Eta%d", iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaLFStats1Up[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaLFStats1Down[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightLFStats2_Pt%d_Eta%d", iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaLFStats2Up[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaLFStats2Down[iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightCErr1_Pt%d_Eta%d"   , iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaCErr1Up   [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaCErr1Down   [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
    systName=Form("CMS_bTagWeightCErr2_Pt%d_Eta%d"   , iPt, iEta); card<<(systName+" shape ").Data(); for(int ic=ic0; ic<nPlotCategories; ic++) { card<<"1.0 "; histo_cmvaCErr2Up   [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Up"); histo_cmvaCErr2Down   [iPt][iEta][ic]->Write(shapeName[ic]+"_"+systName+"Down"); } card<<std::endl;
  }
  //// Total Jet Energy Scale Shape Uncertainty
  card<<Form("CMS_scale_j shape "); for(int ic=ic0; ic<nPlotCategories; ic++) {
    card<<"1.0 "; histo_jesUp[ic]->Write(shapeName[ic]+"_CMS_scale_jUp"); histo_jesDown[ic]->Write(shapeName[ic]+"_CMS_scale_jDown");
  } card<<std::endl;
  // WPT Corr Errors
  if(useWptCorr) { 
    card<<Form("empiricalWptCorrTT shape "); for(int ic=ic0; ic<nPlotCategories; ic++) {
      if(ic==kPlotTT) {
        card<<"1.0 "; histo_wptCorrUp[ic]->Write(shapeName[ic]+"_empiricalWptCorrTTUp"); histo_wptCorrDown[ic]->Write(shapeName[ic]+"_empiricalWptCorrTTDown");
      } else card<<"- ";
    } card<<std::endl;
    card<<Form("empiricalWptCorrWHF shape "); for(int ic=ic0; ic<nPlotCategories; ic++) {
      if(ic==kPlotWbb || ic==kPlotWb || ic==kPlotTop) {
        card<<"1.0 "; histo_wptCorrUp[ic]->Write(shapeName[ic]+"_empiricalWptCorrWHFUp"); histo_wptCorrDown[ic]->Write(shapeName[ic]+"_empiricalWptCorrWHFDown");
      } else card<<"- ";
    } card<<std::endl;
    card<<Form("empiricalWptCorrWLF shape "); for(int ic=ic0; ic<nPlotCategories; ic++) {
      if(ic==kPlotWLF) {
        card<<"1.0 "; histo_wptCorrUp[ic]->Write(shapeName[ic]+"_empiricalWptCorrWLFUp"); histo_wptCorrDown[ic]->Write(shapeName[ic]+"_empiricalWptCorrWLFDown");
      } else card<<"- ";
    } card<<std::endl;
  }
  // Muon SF Shape
  if(theLepSel!=2) { card<<Form("CMS_eff2016_m shape "); for(int ic=ic0; ic<nPlotCategories; ic++) {
    card<<"1.0 "; histo_muSFUp[ic]->Write(shapeName[ic]+"_CMS_eff2016_mUp"); histo_muSFDown[ic]->Write(shapeName[ic]+"_CMS_eff2016_mDown");
  } card<<std::endl; }
  // Electron SF Shape
  if(theLepSel!=1) { card<<Form("CMS_eff2016_e shape "); for(int ic=ic0; ic<nPlotCategories; ic++) {
    card<<"1.0 "; histo_eleSFUp[ic]->Write(shapeName[ic]+"_CMS_eff2016_eUp"); histo_eleSFDown[ic]->Write(shapeName[ic]+"_CMS_eff2016_eDown");
  } card<<std::endl; }

  // Systematics -- Normalization Uncertainties
  // PDF Acceptance
  for(int ic=ic0; ic<nPlotCategories; ic++) {
    if(!(ic==kPlotVZbb||ic==kPlotVVLF||ic==kPlotWbb||ic==kPlotWb||ic==kPlotWLF||ic==kPlotZbb||ic==kPlotZb||ic==kPlotZLF||ic==kPlotVH)) continue;
    card<<Form("pdf_ACCEPT_%s lnN ",plotBaseNames[ic].Data()); 
    for(int ic2=ic0; ic2<nPlotCategories; ic2++) {
      if(ic2==ic) card<<pdfAcceptUncs[ic]<<" "; 
      else card<<"- ";
    }
    card<<std::endl;
  } 
  // Muon Energy Scale Norm
  if(theLepSel!=2) { card<<Form("CMS_scale2016_m lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< "1.01 "; card <<std::endl; }
  // Electron Energy Scale Norm
  if(theLepSel!=1) { card<<Form("CMS_scale2016_e lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< "1.01 "; card<<std::endl; }
  // Trigger
  card<<Form("CMS_trigger2016 lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< "1.01 "; card << std::endl;
  // Lumi
  card<<Form("lumi_13TeV2016 lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< "1.025 "; card << std::endl;
  // VH Cross Section
  card<<Form("pdf_qqbar lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< ((ic==kPlotVZbb||ic==kPlotVVLF||ic==kPlotVH)?"1.01 ":"- "); card << std::endl;
  card<<Form("pdf_gg lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< ((ic==kPlotTop)?"1.01 ":"- "); card << std::endl;
  // Single Top Cross Section
  card<<Form("CMS_VH_TopNorm lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotTop? "1.15 ":"- "); card<<std::endl;
  // Diboson Cross Section
  card<<Form("CMS_VH_VVNorm lnN "); for(int ic=ic0; ic<nPlotCategories; ic++) card<< ((ic==kPlotVZbb||ic==kPlotVVLF)? "1.15 ":"- "); card<<std::endl;
  if(selection>=kWHLightFlavorCR && selection<=kWHSR) {
    card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",ttbar,ttbar) << std::endl;
    card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",Wbb  ,Wbb  ) << std::endl;
    card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",Wb   ,Wb   ) << std::endl;
    card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",WLF  ,WLF  ) << std::endl;
  } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
    // Additional uncertainty on W+LF to pass the MSD cut
    /*
      passMSD_WLF_Wln rateParam * WLF (@0*1.0) SF_passMSD_WLF
      SF_passMSD_WLF param 1.15 0.25
    */
    //card << Form("passMSD_%s_Wln rateParam * %s (@0*1.0) SF_passMSD_WLF", WLF,WLF) << std::endl;
    //card << Form("SF_passMSD_WLF param 0.87 0.25") << std::endl;

    // Double B-tag scale factor extraction in situ
    // A priori MC efficiencies for the double B cut
    
    // Mar 22 2018 (before switching to corrected MSD)
    //   Process      Eff.(0ijB)     Eff.(1+ijB)
    // -----------------------------------------
    //     ttbar           0.128           0.114
    //      W+bb           0.394           0.178
    //       W+b           0.229           0.104
    //   W+udcsg           0.021           0.026

    if(!useIgnorantVHFSFs) {
      card << "effDoubleB_TT   param 0.13  0.013" << std::endl; 
      card << "effDoubleB_Wbb  param 0.39  0.039"  << std::endl; 
      card << "effDoubleB_Wb   param 0.23  0.023"  << std::endl; 
      card << "effDoubleB_WLF  param 0.024 0.0024" << std::endl; 
      // Efficiency scale factors for the double B cut
      card << "effSFDoubleB_TT  extArg 1.0 [0.1,10]" << std::endl;
      card << "effSFDoubleB_Wbb extArg 1.0 [0.1,10]" << std::endl;
      card << "effSFDoubleB_Wb  extArg 1.0 [0.1,10]" << std::endl;
      card << "effSFDoubleB_WLF extArg 1.0 [0.1,10]" << std::endl;
      // Passing B-tag scale factor rateparam
      if(selection==kWHHeavyFlavorFJCR || selection==kWHTT2bFJCR || selection==kWHFJSR) {
        card << Form("passBB_TT  rateParam * %s (@0*1.0) effSFDoubleB_TT"  ,ttbar) << std::endl;
        card << Form("passBB_Wbb rateParam * %s (@0*1.0) effSFDoubleB_Wbb" ,Wbb  ) << std::endl;
        card << Form("passBB_Wb  rateParam * %s (@0*1.0) effSFDoubleB_Wb " ,Wb   ) << std::endl;
        card << Form("passBB_WLF rateParam * %s (@0*1.0) effSFDoubleB_WLF" ,WLF  ) << std::endl;
        
      } else if(selection==kWHLightFlavorFJCR || selection==kWHTT1bFJCR) {
        card << Form("failBB_TT  rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_TT,effDoubleB_TT",ttbar) << std::endl;
        card << Form("failBB_Wbb rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_Wbb,effDoubleB_Wbb",Wbb  ) << std::endl;
        card << Form("failBB_Wb  rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_Wb,effDoubleB_Wb"  ,Wb   ) << std::endl;
        card << Form("failBB_WLF rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_WLF,effDoubleB_WLF",WLF  ) << std::endl;
      }
    } else { // Ignorant scale factor method for W+HF.
      // 25% uncertainty on the passing, scale by efficiency for the failing
      if(selection==kWHHeavyFlavorFJCR || selection==kWHFJSR) {
        card << "ignorant_passBB_0ijb_Wbb lnN "; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWbb? "1.25/0.75 ":"- "); card<<std::endl;
        card << "ignorant_passBB_0ijb_Wb lnN " ; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWb ? "1.25/0.75 ":"- "); card<<std::endl;
      } else if(selection==kWHTT2bFJCR) {
        card << "ignorant_passBB_1ijb_Wbb lnN "; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWbb? "1.25/0.75 ":"- "); card<<std::endl;
        card << "ignorant_passBB_1ijb_Wb lnN " ; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWb ? "1.25/0.75 ":"- "); card<<std::endl;
      } else if(selection==kWHLightFlavorFJCR) {
        card << "ignorant_passBB_0ijb_Wbb lnN "; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWbb? "0.84/1.16 ":"- "); card<<std::endl;
        card << "ignorant_passBB_0ijb_Wb lnN " ; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWb ? "0.93/1.07 ":"- "); card<<std::endl;
      } else if(selection==kWHTT1bFJCR) {
        card << "ignorant_passBB_1ijb_Wbb lnN "; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWbb? "0.95/1.05 ":"- "); card<<std::endl;
        card << "ignorant_passBB_1ijb_Wb lnN " ; for(int ic=ic0; ic<nPlotCategories; ic++) card<< (ic==kPlotWb ? "0.97/1.03 ":"- "); card<<std::endl;
      }
      // Intelligent scale factors for tt, W+LF
      card << "effDoubleB_TT   param 0.13  0.013" << std::endl; 
      card << "effDoubleB_WLF  param 0.024 0.0024" << std::endl; 
      card << "effSFDoubleB_TT  extArg 1.0 [0.1,10]" << std::endl;
      card << "effSFDoubleB_WLF extArg 1.0 [0.1,10]" << std::endl;
      if(selection==kWHHeavyFlavorFJCR || selection==kWHTT2bFJCR || selection==kWHFJSR) {
        card << Form("passBB_TT  rateParam * %s (@0*1.0) effSFDoubleB_TT"  ,ttbar) << std::endl;
        card << Form("passBB_WLF rateParam * %s (@0*1.0) effSFDoubleB_WLF" ,WLF  ) << std::endl;
      } else if(selection==kWHLightFlavorFJCR || selection==kWHTT1bFJCR) {
        card << Form("failBB_TT  rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_TT,effDoubleB_TT",ttbar) << std::endl;
        card << Form("failBB_WLF rateParam * %s ((1.0-@0*@1)/(1.0-@1)) effSFDoubleB_WLF,effDoubleB_WLF",WLF  ) << std::endl;
      }
    }
  }
 

}

