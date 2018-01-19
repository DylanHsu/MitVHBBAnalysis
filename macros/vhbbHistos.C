#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PandaAnalysis/Utilities/src/CSVHelper.cc"
#include <Compression.h>
#include <TCanvas.h>
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

using namespace vhbbPlot;
void vhbbHistos(
  TString plotTreeFileName,
  TString outputFileName,
  TString dataCardDir,
  vhbbPlot::selectionType selection,
  int theLepSel=-1, // -1: any; 0, e-mu; 1, mu only; 2, ele only
  bool isBlinded=true,
  int MVAVarType=3, //1=Higgs pT; 2=Multiclass BDT; 3=Single-Class BDT;
  TString bdtWeightsResolved="weights/bdt_BDT_multiClass_resolved_jan10_I_test2.weights.xml",
  TString bdtWeightsBoosted="weights/bdt_BDT_multiClass_boosted_jan10_I_test2.weights.xml"
) {
  if(dataCardDir!="") system(Form("mkdir -p MitVHBBAnalysis/datacards/%s",dataCardDir.Data()));
  vector<float> MVAbins; TString MVAVarName;
  if(MVAVarType==1) {
    MVAVarName="Higgs p_{T} classifier [GeV]";
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) MVAbins={100,120,140,160,180,200,250,300,350};
    else if(selection==kWHSR)                                                              MVAbins={100,120,140,160,180,200,250,300,350};
    else throw std::runtime_error("bad selection argument");
  } else if(MVAVarType==2) {
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) {
      MVAbins={-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1}; 
      MVAVarName="Discriminator: Lesser CMVA";
    } else if(selection==kWHLightFlavorFJCR || selection==kWHHeavyFlavorFJCR || selection==kWH2TopFJCR) {
      MVAbins={0,20,40,60,80,100,120,140,160,180,200,220,240};
      MVAVarName="Discriminator: Fatjet soft drop mass";
    } else if(selection==kWHSR || selection==kWHFJSR) {
      //MVAbins={0.2,0.4,0.6,0.7,0.75,0.8,0.85,0.9,0.925,0.95};
      MVAbins={0.20,0.40,0.60,0.70,0.80,0.85,0.9,0.95,1.00};
      MVAVarName="Multiclass BDT";
    }
    else throw std::runtime_error("bad selection argument");
  } else if(MVAVarType==3) {
    if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) {
      MVAbins={-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1}; 
      MVAVarName="Discriminator: Lesser CMVA";
    } else if(selection==kWHLightFlavorFJCR || selection==kWHHeavyFlavorFJCR || selection==kWH2TopFJCR) {
      MVAbins={0,20,40,60,80,100,120,140,160,180,200,220,240};
      MVAVarName="Discriminator: Fatjet soft drop mass";
    } else if(selection==kWHSR || selection==kWHFJSR) {
      MVAbins={0.20,0.40,0.60,0.70,0.80,0.85,0.9,0.95,1.00};
      MVAVarName="Single class BDT";
    }
    else throw std::runtime_error("bad selection argument");
  } else throw std::runtime_error("bad MVAVarType");
  
  // Load Files
  gSystem->Load("libCondFormatsJetMETObjects.so");
  system("mkdir -p MitVHBBAnalysis/datacards");
  system("mkdir -p MitVHBBAnalysis/plots");
  TFile *plotTreeFile = TFile::Open(plotTreeFileName, "READ");
  assert(plotTreeFile && plotTreeFile->IsOpen());
  string theJecUncertainties = "PandaAnalysis/data/jec/23Sep2016V4/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt";
  JetCorrectionUncertainty *jcuTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters(theJecUncertainties, "Total")));
  CSVHelper *cmvaReweighter = new CSVHelper("PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root"   , "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root"   , 5);
  
  // Initialize the tree
  TTree *plotTree = (TTree*)plotTreeFile->Get("plotTree"); assert(plotTree);
  unsigned char typeLepSel, theCategory;
  unsigned selectionBits, selectionBits_jesUp, selectionBits_jesDown, nMinusOneBits;
  int nJet, nSoft2, nSoft5, nSoft10, nIsojet;
  float sumEtSoft1;
  float pfmet, pfmetphi, pfmetsig;
  float lepton1Pt, lepton1Eta, lepton1Phi, lepton1RelIso;
  float lepton2Pt, lepton2Eta, lepton2Phi;
  float hbbJet1Pt, hbbJet1Eta, hbbJet1Phi;
  float hbbJet2Pt, hbbJet2Eta, hbbJet2Phi;
  float topWBosonPt, topWBosonEta, topWBosonPhi, topWBosonCosThetaCS;
  float hbbDijetPt, hbbDijetMass, hbbCosThetaJJ, hbbCosThetaCSJ1;
  float bDiscrMin, bDiscrMax;
  float deltaPhiLep1Met, deltaPhiVH;
  float topMassLep1Met;
  float weight;
  float weight_pdfUp, weight_pdfDown;
  float weight_QCDr1f2, weight_QCDr1f5, weight_QCDr2f1, weight_QCDr2f2, weight_QCDr5f1, weight_QCDr5f5, weight_lepSFUp;
  float weight_cmvaJESUp     [5][3] , weight_cmvaJESDown     [5][3] , 
        weight_cmvaLFUp      [5][3] , weight_cmvaLFDown      [5][3] , 
        weight_cmvaHFUp      [5][3] , weight_cmvaHFDown      [5][3] , 
        weight_cmvaHFStats1Up[5][3] , weight_cmvaHFStats1Down[5][3] , 
        weight_cmvaHFStats2Up[5][3] , weight_cmvaHFStats2Down[5][3] , 
        weight_cmvaLFStats1Up[5][3] , weight_cmvaLFStats1Down[5][3] , 
        weight_cmvaLFStats2Up[5][3] , weight_cmvaLFStats2Down[5][3] , 
        weight_cmvaCErr1Up   [5][3] , weight_cmvaCErr1Down   [5][3] , 
        weight_cmvaCErr2Up   [5][3] , weight_cmvaCErr2Down   [5][3] ; 
  float mT;
  float fj1Tau32;
  float fj1Tau21;
  float fj1Tau32SD;
  float fj1Tau21SD;
  float fj1MSD;
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
  float fj1Eta;
  float fj1M;
  float fj1MaxCSV;
  float fj1MinCSV;
  float fj1DoubleCSV;
  float fj1HTTMass;
  float fj1HTTFRec;
  float fj1SDEFrac100;
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
  mvaPsi["012004011003"] = -3393; 
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


  plotTree->SetBranchAddress("selectionBits"   , &selectionBits   );
  plotTree->SetBranchAddress("selectionBits_jesUp"     , &selectionBits_jesUp   );
  plotTree->SetBranchAddress("selectionBits_jesDown"   , &selectionBits_jesDown );
  plotTree->SetBranchAddress("nMinusOneBits"   , &nMinusOneBits   );
  plotTree->SetBranchAddress("theCategory"     , &theCategory     );
  plotTree->SetBranchAddress("pfmet"           , &pfmet           );
  plotTree->SetBranchAddress("pfmetsig"        , &pfmetsig        );
  plotTree->SetBranchAddress("pfmetphi"        , &pfmetphi        );
  plotTree->SetBranchAddress("weight"          , &weight          );
  plotTree->SetBranchAddress("weight_pdfUp"       , &weight_pdfUp        );
  plotTree->SetBranchAddress("weight_pdfDown"     , &weight_pdfDown      );
  plotTree->SetBranchAddress("weight_QCDr1f2"     , &weight_QCDr1f2      );
  plotTree->SetBranchAddress("weight_QCDr1f5"     , &weight_QCDr1f5      );
  plotTree->SetBranchAddress("weight_QCDr2f1"     , &weight_QCDr2f1      );
  plotTree->SetBranchAddress("weight_QCDr2f2"     , &weight_QCDr2f2      );
  plotTree->SetBranchAddress("weight_QCDr5f1"     , &weight_QCDr5f1      );
  plotTree->SetBranchAddress("weight_QCDr5f5"     , &weight_QCDr5f5      );
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
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    plotTree->SetBranchAddress("nJet"               , &nJet                );
    plotTree->SetBranchAddress("typeLepSel"         , &typeLepSel          );
    plotTree->SetBranchAddress("hbbDijetPt"         , &hbbDijetPt          );
    plotTree->SetBranchAddress("hbbDijetMass"       , &hbbDijetMass        );
    plotTree->SetBranchAddress("hbbCosThetaJJ"      , &hbbCosThetaJJ       );
    plotTree->SetBranchAddress("hbbCosThetaCSJ1"    , &hbbCosThetaCSJ1     );
    plotTree->SetBranchAddress("topMassLep1Met"     , &topMassLep1Met      );
    plotTree->SetBranchAddress("bDiscrMin"          , &bDiscrMin           );
    plotTree->SetBranchAddress("bDiscrMax"          , &bDiscrMax           );
    plotTree->SetBranchAddress("hbbJet1Pt"          , &hbbJet1Pt           );
    plotTree->SetBranchAddress("hbbJet1Eta"         , &hbbJet1Eta          );
    plotTree->SetBranchAddress("hbbJet1Phi"         , &hbbJet1Phi          );
    plotTree->SetBranchAddress("hbbJet2Pt"          , &hbbJet2Pt           );
    plotTree->SetBranchAddress("hbbJet2Eta"         , &hbbJet2Eta          );
    plotTree->SetBranchAddress("hbbJet2Phi"         , &hbbJet2Phi          );
    plotTree->SetBranchAddress("lepton1Pt"          , &lepton1Pt           );
    plotTree->SetBranchAddress("lepton1Eta"         , &lepton1Eta          );
    plotTree->SetBranchAddress("lepton1Phi"         , &lepton1Phi          );
    plotTree->SetBranchAddress("topWBosonCosThetaCS", &topWBosonCosThetaCS );
    plotTree->SetBranchAddress("topWBosonPt"        , &topWBosonPt         );
    plotTree->SetBranchAddress("topWBosonEta"       , &topWBosonEta        );
    plotTree->SetBranchAddress("topWBosonPhi"       , &topWBosonPhi        );
    plotTree->SetBranchAddress("mT"                 , &mT                  );
    plotTree->SetBranchAddress("deltaPhiLep1Met"    , &deltaPhiLep1Met     );
    plotTree->SetBranchAddress("deltaPhiVH"         , &deltaPhiVH          );
    plotTree->SetBranchAddress("weight_lepSFUp"     , &weight_lepSFUp      );
    plotTree->SetBranchAddress("nSoft2"             , &nSoft2              );
    plotTree->SetBranchAddress("nSoft5"             , &nSoft5              );
    plotTree->SetBranchAddress("nSoft10"            , &nSoft10             );
    plotTree->SetBranchAddress("sumEtSoft1"         , &sumEtSoft1          );
  } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
    plotTree->SetBranchAddress("topWBosonPt"           , &topWBosonPt           );
    plotTree->SetBranchAddress("topWBosonPhi"          , &topWBosonPhi          );
    plotTree->SetBranchAddress("deltaPhiLep1Met"       , &deltaPhiLep1Met       );
    plotTree->SetBranchAddress("deltaPhiVH"            , &deltaPhiVH            );
    plotTree->SetBranchAddress("typeLepSel"            , &typeLepSel            );
    plotTree->SetBranchAddress("lepton1Pt"             , &lepton1Pt             );
    plotTree->SetBranchAddress("lepton1Eta"            , &lepton1Eta            );
    plotTree->SetBranchAddress("lepton1Phi"            , &lepton1Phi            );
    plotTree->SetBranchAddress("lepton1RelIso"         , &lepton1RelIso         );
    plotTree->SetBranchAddress("nIsojet"               , &nIsojet               );
    plotTree->SetBranchAddress("mT"                    , &mT                    );
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
    plotTree->SetBranchAddress("weight_lepSFUp"        , &weight_lepSFUp           );
  }
  float balanceVH=-3393;
  float dPhiHbbJet1MET=-3393;
  float dPhiHbbJet2MET=-3393;
  float MVAVar=-3393;
  float nAddJetZeroOrOne;
  float mvaNSoft2, mvaNSoft5, mvaNSoft10;
  float dPhil1W, dPhil1b1, dPhil1b2, dPhiWb1, dPhiWb2, dPhib1b2, dEtal1W, dEtal1b1, dEtal1b2, dEtaWb1, dEtaWb2, dEtab1b2;
  float dPhil1fj1, dPhiWfj1, dEtal1fj1;
  float mva_nIsojet;

  // construct histograms
  int nPlots=99;
  vector<float> xmin(nPlots), xmax(nPlots);  vector<int> nbins(nPlots);
  vector<TString> histoNames(nPlots), histoTitles(nPlots);
  
  int pMVAVar=-1;
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) {
    int p=0;
    histoNames[p]="pTH"                    ; histoTitles[p]="Higgs daughter dijet p_{T} [GeV]"           ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="mH"                     ; histoTitles[p]="Higgs daughter dijet mass [GeV]"            ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="WpT"                    ; histoTitles[p]="W boson p_{T} [GeV]"                        ; nbins[p]=  18; xmin[p]=    50; xmax[p]=   350; p++; 
    histoNames[p]="deltaPhiVH"             ; histoTitles[p]="#Delta#phi(H,W) [Rad]"                      ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="pTBalanceDijetW"        ; histoTitles[p]="Higgs daughter dijet p_{T} / W boson p_{T}" ; nbins[p]=  30; xmin[p]=     0; xmax[p]=     2; p++; 
    histoNames[p]="lepton1Pt"              ; histoTitles[p]="Lepton p_{T} [GeV]"                         ; nbins[p]=  23; xmin[p]=    20; xmax[p]=   250; p++; 
    histoNames[p]="pfmet"                  ; histoTitles[p]="p_{T}^{miss} [GeV]"                         ; nbins[p]=  25; xmin[p]=     0; xmax[p]=   250; p++; 
    histoNames[p]="mTW"                    ; histoTitles[p]="W transverse mass [GeV]"                    ; nbins[p]=  20; xmin[p]=     0; xmax[p]=   200; p++; 
    histoNames[p]="Hbjet1Pt"               ; histoTitles[p]="Leading H daughter jet p_{T} [GeV]"         ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="Hbjet2Pt"               ; histoTitles[p]="Trailing H daughter jet p_{T} [GeV]"        ; nbins[p]=  38; xmin[p]=    20; xmax[p]=   400; p++; 
    histoNames[p]="dPhiHbbJet1MET"         ; histoTitles[p]="#Delta#phi(H daughter jet 1, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="dPhiHbbJet2MET"         ; histoTitles[p]="#Delta#phi(H daughter jet 2, E_{T}^{miss})" ; nbins[p]=  32; xmin[p]=     0; xmax[p]= 3.142; p++; 
    histoNames[p]="bDiscrMin"              ; histoTitles[p]="Lesser CMVA of Higgs daughter jets"         ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="bDiscrMax"              ; histoTitles[p]="Greater CMVA of Higgs daughter jets"        ; nbins[p]=  40; xmin[p]=   -1.; xmax[p]=    1.; p++; 
    histoNames[p]="topMassLep1Met"         ; histoTitles[p]="Reco top mass [GeV]"                        ; nbins[p]=  22; xmin[p]=   60.; xmax[p]=  500.; p++; 
    histoNames[p]="nAddJetZeroOrOne"       ; histoTitles[p]="Has more than 2 jets"                       ; nbins[p]=   2; xmin[p]=    0.; xmax[p]=    2.; p++; 
    histoNames[p]="sumEtSoft1"             ; histoTitles[p]="#sum E_{T}(soft 1)"                         ; nbins[p]=  40; xmin[p]=    0.; xmax[p]=  200.; p++; 
    histoNames[p]="nSoft2"                 ; histoTitles[p]="N^{soft}_{2}"                               ; nbins[p]=  25; xmin[p]=    0.; xmax[p]=   25.; p++; 
    histoNames[p]="nSoft5"                 ; histoTitles[p]="N^{soft}_{5}"                               ; nbins[p]=  10; xmin[p]=    0.; xmax[p]=   10.; p++; 
    histoNames[p]="nSoft10"                ; histoTitles[p]="N^{soft}_{10}"                              ; nbins[p]=   5; xmin[p]=    0.; xmax[p]=    5.; p++; 
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
    histoNames[p]="MVAVar"                 ; histoTitles[p]=MVAVarName; pMVAVar=p; p++;
  } else if(selection>=kWHLightFlavorFJCR && selection<=kWHFJSR) {
    int p=0;
    histoNames[p]="nIsojet"                ; histoTitles[p]="N central AK4 jets"                         ; nbins[p]=   7; xmin[p]=    0.; xmax[p]=    7.; p++;
    histoNames[p]="fj1Pt"                  ; histoTitles[p]="FJ p_{T}"                                   ; nbins[p]=  20; xmin[p]=  200.; xmax[p]= 2200.; p++;
    histoNames[p]="fj1Eta"                 ; histoTitles[p]="FJ #eta"                                    ; nbins[p]=  20; xmin[p]=  -2.5; xmax[p]=   2.5; p++;
    histoNames[p]="greaterSubjetCSV"       ; histoTitles[p]="FJ greater subjet B-tag (CSV)"              ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="lesserSubjetCSV"        ; histoTitles[p]="FJ lesser subjet B-tag (CSV)"               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1DoubleCSV"           ; histoTitles[p]="FJ double B-tag"                            ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1HTTMass"             ; histoTitles[p]="FJ HTT mass"                                ; nbins[p]=  13; xmin[p]=  260.; xmax[p]=    5.; p++;
    histoNames[p]="fj1HTTFRec"             ; histoTitles[p]="FJ HTT f_{rec}"                             ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    5.; p++;
    histoNames[p]="fj1MSD"                 ; histoTitles[p]="FJ soft drop mass"                          ; nbins[p]=  30; xmin[p]=    0.; xmax[p]=  300.; p++;
    histoNames[p]="fj1Tau32"               ; histoTitles[p]="FJ #tau_{32}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau32SD"             ; histoTitles[p]="FJ #tau_{32}^{SD}"                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau21"               ; histoTitles[p]="FJ #tau_{21}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="fj1Tau21SD"             ; histoTitles[p]="FJ #tau_{21}^{SD}"                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    1.; p++;
    histoNames[p]="lepton1Pt"              ; histoTitles[p]="Lepton p_{T}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2000.; p++;
    histoNames[p]="deltaPhiVH"             ; histoTitles[p]="#Delta#phi(W,H)"                            ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3.142; p++;
    histoNames[p]="deltaPhiLep1Met"        ; histoTitles[p]="#Delta#phi(lepton,p_{T}^{miss})"            ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    5.; p++;
    histoNames[p]="topWBosonPt"            ; histoTitles[p]="W boson p_{T}"                              ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2000.; p++;
    histoNames[p]="mT"                     ; histoTitles[p]="W boson m_{T}"                              ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  400.; p++;
    histoNames[p]="pfmetsig"               ; histoTitles[p]="E_{T}^{miss} sig."                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=   20.; p++;
    histoNames[p]="pfmet"                  ; histoTitles[p]="p_{T}^{miss}"                               ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 2000.; p++;
    histoNames[p]="dPhil1W"                ; histoTitles[p]="#Delta#phi(lep,W)"                          ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3.142; p++;
    histoNames[p]="dPhil1fj1"              ; histoTitles[p]="#Delta#phi(lep,FJ)"                         ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3.142; p++;
    histoNames[p]="dPhiWfj1"               ; histoTitles[p]="#Delta#phi(W,FJ)"                           ; nbins[p]=  20; xmin[p]=    0.; xmax[p]= 3.142; p++;
    histoNames[p]="dEtal1fj1"              ; histoTitles[p]="|#Delta#eta(lep,FJ)|"                       ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=    5.; p++;
    histoNames[p]="psi021004010502"        ; histoTitles[p]="#psi(2,1.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.5 ; p++;
    histoNames[p]="psi012004010502"        ; histoTitles[p]="#psi(1,2.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.3 ; p++;
    histoNames[p]="psi021003011002"        ; histoTitles[p]="#psi(2,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.4 ; p++;
    histoNames[p]="psi022004011002"        ; histoTitles[p]="#psi(2,2.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.3 ; p++;
    histoNames[p]="psi020503010502"        ; histoTitles[p]="#psi(2,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.5 ; p++;
    histoNames[p]="psi022003012002"        ; histoTitles[p]="#psi(2,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.4 ; p++;
    histoNames[p]="psi012004020503"        ; histoTitles[p]="#psi(1,2.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  2.0 ; p++;
    histoNames[p]="psi021004020503"        ; histoTitles[p]="#psi(2,1.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  3.0 ; p++;
    histoNames[p]="psi032003012002"        ; histoTitles[p]="#psi(3,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  20.0; p++;
    histoNames[p]="psi012003011002"        ; histoTitles[p]="#psi(1,2.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.3 ; p++;
    histoNames[p]="psi011003010502"        ; histoTitles[p]="#psi(1,1.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.4 ; p++;
    histoNames[p]="psi031003011002"        ; histoTitles[p]="#psi(3,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  8.0 ; p++;
    histoNames[p]="psi031003012002"        ; histoTitles[p]="#psi(3,1.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.4 ; p++;
    histoNames[p]="psi030503011002"        ; histoTitles[p]="#psi(3,0.5,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.4 ; p++;
    histoNames[p]="psi022004012002"        ; histoTitles[p]="#psi(2,2.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.0 ; p++;
    histoNames[p]="psi021004011002"        ; histoTitles[p]="#psi(2,1.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi012004011002"        ; histoTitles[p]="#psi(1,2.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.0 ; p++;
    histoNames[p]="psi022004021003"        ; histoTitles[p]="#psi(2,2.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  2.0 ; p++;
    histoNames[p]="psi012003010502"        ; histoTitles[p]="#psi(1,2.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  30.0; p++;
    histoNames[p]="psi022003011002"        ; histoTitles[p]="#psi(2,2.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  50.0; p++;
    histoNames[p]="psi022004031003"        ; histoTitles[p]="#psi(2,2.0,4,3,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi012004030503"        ; histoTitles[p]="#psi(1,2.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi021004030503"        ; histoTitles[p]="#psi(2,1.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi020504010502"        ; histoTitles[p]="#psi(2,0.5,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi022004030503"        ; histoTitles[p]="#psi(2,2.0,4,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  50.0; p++;
    histoNames[p]="psi011004010502"        ; histoTitles[p]="#psi(1,1.0,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi030503010502"        ; histoTitles[p]="#psi(3,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  6.0 ; p++;
    histoNames[p]="psi021003010502"        ; histoTitles[p]="#psi(2,1.0,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  50.0; p++;
    histoNames[p]="psi012004011003"        ; histoTitles[p]="#psi(1,2.0,4,1,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  4.0 ; p++;
    histoNames[p]="psi012004021004"        ; histoTitles[p]="#psi(1,2.0,4,2,1.0,4)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.8 ; p++;
    histoNames[p]="psi021004012004"        ; histoTitles[p]="#psi(2,1.0,4,1,2.0,4)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  4.0 ; p++;
    histoNames[p]="psi011003020503"        ; histoTitles[p]="#psi(1,1.0,3,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.9 ; p++;
    histoNames[p]="psi020503011003"        ; histoTitles[p]="#psi(2,0.5,3,1,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  3.0 ; p++;
    histoNames[p]="psi012003030503"        ; histoTitles[p]="#psi(1,2.0,3,3,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  3.0 ; p++;
    histoNames[p]="psi030503012003"        ; histoTitles[p]="#psi(3,0.5,3,1,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  5.0 ; p++;
    histoNames[p]="psi010503010502"        ; histoTitles[p]="#psi(1,0.5,3,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.2 ; p++;
    histoNames[p]="psi020503011002"        ; histoTitles[p]="#psi(2,0.5,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.2 ; p++;
    histoNames[p]="psi011003011002"        ; histoTitles[p]="#psi(1,1.0,3,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.2 ; p++;
    histoNames[p]="psi020504011002"        ; histoTitles[p]="#psi(2,0.5,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.0 ; p++;
    histoNames[p]="psi012004021003"        ; histoTitles[p]="#psi(1,2.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi012004012002"        ; histoTitles[p]="#psi(1,2.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.0 ; p++;
    histoNames[p]="psi020504020503"        ; histoTitles[p]="#psi(2,0.5,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.2 ; p++;
    histoNames[p]="psi011004011002"        ; histoTitles[p]="#psi(1,1.0,4,1,1.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.0 ; p++;
    histoNames[p]="psi021004012002"        ; histoTitles[p]="#psi(2,1.0,4,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.0 ; p++;
    histoNames[p]="psi021003012002"        ; histoTitles[p]="#psi(2,1.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.2 ; p++;
    histoNames[p]="psi011004020503"        ; histoTitles[p]="#psi(1,1.0,4,2,0.5,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi010504010502"        ; histoTitles[p]="#psi(1,0.5,4,1,0.5,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi021004021003"        ; histoTitles[p]="#psi(2,1.0,4,2,1.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi012003012002"        ; histoTitles[p]="#psi(1,2.0,3,1,2.0,2)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="psi022004022003"        ; histoTitles[p]="#psi(2,2.0,4,2,2.0,3)"                      ; nbins[p]=  20; xmin[p]=    0.; xmax[p]=  0.1 ; p++;
    histoNames[p]="MVAVar"                 ; histoTitles[p]=MVAVarName; pMVAVar=p; p++;
  }
  assert(pMVAVar!=-1);

  TH1D *histos[99][nPlotCategories];
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
  TH1D *histo_pdfUp    [nPlotCategories];
  TH1D *histo_pdfDown  [nPlotCategories];
  TH1D *histo_QCDr1f2  [nPlotCategories];
  TH1D *histo_QCDr1f5  [nPlotCategories];
  TH1D *histo_QCDr2f1  [nPlotCategories];
  TH1D *histo_QCDr2f2  [nPlotCategories];
  TH1D *histo_QCDr5f1  [nPlotCategories];
  TH1D *histo_QCDr5f5  [nPlotCategories];
  TH1D *histo_eleSFUp  [nPlotCategories];
  TH1D *histo_muSFUp   [nPlotCategories];
  TH1D *histo_cmvaJESUp     [5][3][nPlotCategories], *histo_cmvaJESDown     [5][3][nPlotCategories]; 
  TH1D *histo_cmvaLFUp      [5][3][nPlotCategories], *histo_cmvaLFDown      [5][3][nPlotCategories]; 
  TH1D *histo_cmvaHFUp      [5][3][nPlotCategories], *histo_cmvaHFDown      [5][3][nPlotCategories]; 
  TH1D *histo_cmvaHFStats1Up[5][3][nPlotCategories], *histo_cmvaHFStats1Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaHFStats2Up[5][3][nPlotCategories], *histo_cmvaHFStats2Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaLFStats1Up[5][3][nPlotCategories], *histo_cmvaLFStats1Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaLFStats2Up[5][3][nPlotCategories], *histo_cmvaLFStats2Down[5][3][nPlotCategories]; 
  TH1D *histo_cmvaCErr1Up   [5][3][nPlotCategories], *histo_cmvaCErr1Down   [5][3][nPlotCategories]; 
  TH1D *histo_cmvaCErr2Up   [5][3][nPlotCategories], *histo_cmvaCErr2Down   [5][3][nPlotCategories]; 
  TH1D *histo_jesUp   [nPlotCategories];
  TH1D *histo_jesDown [nPlotCategories];
  for(theCategory=0; theCategory<nPlotCategories; theCategory++) {
    histo_pdfUp    [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_pdfUp"    , theCategory)); histo_pdfUp    [theCategory]->SetDirectory(0);
    histo_pdfDown  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_pdfDown"  , theCategory)); histo_pdfDown  [theCategory]->SetDirectory(0);
    histo_QCDr1f2  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr1f2"  , theCategory)); histo_QCDr1f2  [theCategory]->SetDirectory(0);
    histo_QCDr1f5  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr1f5"  , theCategory)); histo_QCDr1f5  [theCategory]->SetDirectory(0);
    histo_QCDr2f1  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr2f1"  , theCategory)); histo_QCDr2f1  [theCategory]->SetDirectory(0);
    histo_QCDr2f2  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr2f2"  , theCategory)); histo_QCDr2f2  [theCategory]->SetDirectory(0);
    histo_QCDr5f1  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr5f1"  , theCategory)); histo_QCDr5f1  [theCategory]->SetDirectory(0);
    histo_QCDr5f5  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_QCDr5f5"  , theCategory)); histo_QCDr5f5  [theCategory]->SetDirectory(0);
    histo_eleSFUp  [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_eleSFUp"  , theCategory)); histo_eleSFUp  [theCategory]->SetDirectory(0);
    histo_muSFUp   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_muSFUp"   , theCategory)); histo_muSFUp   [theCategory]->SetDirectory(0);
    for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
      histo_cmvaLFUp        [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFUp_pt%d_eta%d"        , theCategory, iPt, iEta)); histo_cmvaLFUp         [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFUp        [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFUp_pt%d_eta%d"        , theCategory, iPt, iEta)); histo_cmvaHFUp         [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats1Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats1Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaHFStats1Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats2Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats2Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaHFStats2Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats1Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats1Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaLFStats1Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats2Up  [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats2Up_pt%d_eta%d"  , theCategory, iPt, iEta)); histo_cmvaLFStats2Up   [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr1Up     [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr1Up_pt%d_eta%d"     , theCategory, iPt, iEta)); histo_cmvaCErr1Up      [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr2Up     [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr2Up_pt%d_eta%d"     , theCategory, iPt, iEta)); histo_cmvaCErr2Up      [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFDown      [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFDown_pt%d_eta%d"      , theCategory, iPt, iEta)); histo_cmvaLFDown       [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFDown      [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFDown_pt%d_eta%d"      , theCategory, iPt, iEta)); histo_cmvaHFDown       [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats1Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats1Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaHFStats1Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaHFStats2Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaHFStats2Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaHFStats2Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats1Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats1Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaLFStats1Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaLFStats2Down[iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaLFStats2Down_pt%d_eta%d", theCategory, iPt, iEta)); histo_cmvaLFStats2Down [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr1Down   [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr1Down_pt%d_eta%d"   , theCategory, iPt, iEta)); histo_cmvaCErr1Down    [iPt][iEta][theCategory]->SetDirectory(0);
      histo_cmvaCErr2Down   [iPt][iEta][theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_cmvaCErr2Down_pt%d_eta%d"   , theCategory, iPt, iEta)); histo_cmvaCErr2Down    [iPt][iEta][theCategory]->SetDirectory(0);
    }
    histo_jesUp   [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_jesUp"   , theCategory)); histo_jesUp   [theCategory]->SetDirectory(0);
    histo_jesDown [theCategory] = (TH1D*)histos[pMVAVar][theCategory]->Clone(Form("histo%d_jesDown" , theCategory)); histo_jesDown [theCategory]->SetDirectory(0);
  }
  // done constructing histograms

  // Initialize mva reader
  std::map <string, float*> mvaInputRefs;
  string str_greaterSubjetCSV = "fj1MaxCSV!=fj1MaxCSV?0.:TMath::Min(1.,TMath::Max(0.,fj1MaxCSV))";
  string str_lesserSubjetCSV  = "fj1MinCSV!=fj1MinCSV?0.:TMath::Min(1.,TMath::Max(0.,fj1MinCSV))";
  mvaInputRefs["hbbDijetMass"                                        ]=&hbbDijetMass        ;    
  mvaInputRefs["hbbDijetPt"                                          ]=&hbbDijetPt          ;    
  mvaInputRefs["topWBosonPt"                                         ]=&topWBosonPt         ;    
  mvaInputRefs["bDiscrMin"                                           ]=&bDiscrMin           ;    
  mvaInputRefs["topMassLep1Met"                                      ]=&topMassLep1Met      ;    
  mvaInputRefs["deltaPhiVH"                                          ]=&deltaPhiVH          ;    
  mvaInputRefs["TMath::Min(nJet-2,1)"                                ]=&nAddJetZeroOrOne    ;    
  mvaInputRefs["deltaPhiLep1Met"                                     ]=&deltaPhiLep1Met     ;    
  mvaInputRefs["mT"                                                  ]=&mT                  ;    
  mvaInputRefs["pfmet"                                               ]=&pfmet               ;    
  mvaInputRefs["nSoft5"                                              ]=&mvaNSoft5           ;    
  mvaInputRefs["lepton1Pt"                                           ]=&lepton1Pt           ;    
  mvaInputRefs["bDiscrMax"                                           ]=&bDiscrMax           ;    
  mvaInputRefs["pfmetsig"                                            ]=&pfmetsig            ;    
  mvaInputRefs["sumEtSoft1"                                          ]=&sumEtSoft1          ;    
  mvaInputRefs["nSoft2"                                              ]=&mvaNSoft2           ;    
  mvaInputRefs["nSoft10"                                             ]=&mvaNSoft10          ;    
  mvaInputRefs["topWBosonCosThetaCS"                                 ]=&topWBosonCosThetaCS ;    
  mvaInputRefs["hbbCosThetaJJ"                                       ]=&hbbCosThetaJJ       ;    
  mvaInputRefs["hbbCosThetaCSJ1"                                     ]=&hbbCosThetaCSJ1     ;    
  mvaInputRefs["hbbJet1Pt"                                           ]=&hbbJet1Pt           ;    
  mvaInputRefs["hbbJet2Pt"                                           ]=&hbbJet2Pt           ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-topWBosonPhi))" ]=&dPhil1W             ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-hbbJet1Phi))"   ]=&dPhil1b1            ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-hbbJet2Phi))"   ]=&dPhil1b2            ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet1Phi))" ]=&dPhiWb1             ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-hbbJet2Phi))" ]=&dPhiWb2             ;    
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(hbbJet1Phi-hbbJet2Phi))"   ]=&dPhib1b2            ;    
  mvaInputRefs["fabs(lepton1Eta-topWBosonEta)"                       ]=&dEtal1W             ;    
  mvaInputRefs["fabs(lepton1Eta-hbbJet1Eta)"                         ]=&dEtal1b1            ;    
  mvaInputRefs["fabs(lepton1Eta-hbbJet2Eta)"                         ]=&dEtal1b2            ;    
  mvaInputRefs["fabs(topWBosonEta-hbbJet1Eta)"                       ]=&dEtaWb1             ;    
  mvaInputRefs["fabs(topWBosonEta-hbbJet2Eta)"                       ]=&dEtaWb2             ;    
  mvaInputRefs["fabs(hbbJet1Eta-hbbJet2Eta)"                         ]=&dEtab1b2            ;    
  mvaInputRefs["nIsojet"                                             ]=&mva_nIsojet         ;
  mvaInputRefs["fj1Pt"                                               ]=&fj1Pt               ;
  mvaInputRefs["fj1Eta"                                              ]=&fj1Eta              ;
  mvaInputRefs[str_greaterSubjetCSV                                  ]=&fj1MaxCSV           ;
  mvaInputRefs[str_lesserSubjetCSV                                   ]=&fj1MinCSV           ;
  mvaInputRefs["fj1DoubleCSV"                                        ]=&fj1DoubleCSV        ;
  mvaInputRefs["fj1HTTMass"                                          ]=&fj1HTTMass          ;
  mvaInputRefs["fj1HTTFRec"                                          ]=&fj1HTTFRec          ;
  mvaInputRefs["fj1MSD"                                              ]=&fj1MSD              ;
  mvaInputRefs["fj1Tau32"                                            ]=&fj1Tau32            ;
  mvaInputRefs["fj1Tau32SD"                                          ]=&fj1Tau32SD          ;
  mvaInputRefs["fj1Tau21"                                            ]=&fj1Tau21            ;
  mvaInputRefs["fj1Tau21SD"                                          ]=&fj1Tau21SD          ;
  mvaInputRefs["lepton1Pt"                                           ]=&lepton1Pt           ;
  mvaInputRefs["deltaPhiVH"                                          ]=&deltaPhiVH          ;
  mvaInputRefs["deltaPhiLep1Met"                                     ]=&deltaPhiLep1Met     ;
  mvaInputRefs["topWBosonPt"                                         ]=&topWBosonPt         ;
  mvaInputRefs["mT"                                                  ]=&mT                  ;
  mvaInputRefs["pfmetsig"                                            ]=&pfmetsig            ;
  mvaInputRefs["pfmet"                                               ]=&pfmet               ;
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-topWBosonPhi))" ]=&dPhil1W             ;
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(lepton1Phi-fj1Phi))"       ]=&dPhil1fj1           ;
  mvaInputRefs["fabs(TVector2::Phi_mpi_pi(topWBosonPhi-fj1Phi))"     ]=&dPhiWfj1            ;
  mvaInputRefs["fabs(lepton1Eta-fj1Eta)"                             ]=&dEtal1fj1           ;
  mvaInputRefs["psi021004010502"                                     ]=&mvaPsi["021004010502"]; 
  mvaInputRefs["psi012004010502"                                     ]=&mvaPsi["012004010502"]; 
  mvaInputRefs["psi021003011002"                                     ]=&mvaPsi["021003011002"]; 
  mvaInputRefs["psi022004011002"                                     ]=&mvaPsi["022004011002"]; 
  mvaInputRefs["psi020503010502"                                     ]=&mvaPsi["020503010502"]; 
  mvaInputRefs["psi022003012002"                                     ]=&mvaPsi["022003012002"]; 
  mvaInputRefs["psi012004020503"                                     ]=&mvaPsi["012004020503"]; 
  mvaInputRefs["psi021004020503"                                     ]=&mvaPsi["021004020503"]; 
  mvaInputRefs["psi032003012002"                                     ]=&mvaPsi["032003012002"]; 
  mvaInputRefs["psi012003011002"                                     ]=&mvaPsi["012003011002"]; 
  mvaInputRefs["psi011003010502"                                     ]=&mvaPsi["011003010502"]; 
  mvaInputRefs["psi031003011002"                                     ]=&mvaPsi["031003011002"]; 
  mvaInputRefs["psi031003012002"                                     ]=&mvaPsi["031003012002"]; 
  mvaInputRefs["psi030503011002"                                     ]=&mvaPsi["030503011002"]; 
  mvaInputRefs["psi022004012002"                                     ]=&mvaPsi["022004012002"]; 
  mvaInputRefs["psi021004011002"                                     ]=&mvaPsi["021004011002"]; 
  mvaInputRefs["psi012004011002"                                     ]=&mvaPsi["012004011002"]; 
  mvaInputRefs["psi022004021003"                                     ]=&mvaPsi["022004021003"]; 
  mvaInputRefs["psi012003010502"                                     ]=&mvaPsi["012003010502"]; 
  mvaInputRefs["psi022003011002"                                     ]=&mvaPsi["022003011002"]; 
  mvaInputRefs["psi022004031003"                                     ]=&mvaPsi["022004031003"]; 
  mvaInputRefs["psi012004030503"                                     ]=&mvaPsi["012004030503"]; 
  mvaInputRefs["psi021004030503"                                     ]=&mvaPsi["021004030503"]; 
  mvaInputRefs["psi020504010502"                                     ]=&mvaPsi["020504010502"]; 
  mvaInputRefs["psi022004030503"                                     ]=&mvaPsi["022004030503"]; 
  mvaInputRefs["psi011004010502"                                     ]=&mvaPsi["011004010502"]; 
  mvaInputRefs["psi030503010502"                                     ]=&mvaPsi["030503010502"]; 
  mvaInputRefs["psi021003010502"                                     ]=&mvaPsi["021003010502"]; 
  mvaInputRefs["psi012004011003"                                     ]=&mvaPsi["012004011003"]; 
  mvaInputRefs["psi012004021004"                                     ]=&mvaPsi["012004021004"]; 
  mvaInputRefs["psi021004012004"                                     ]=&mvaPsi["021004012004"]; 
  mvaInputRefs["psi011003020503"                                     ]=&mvaPsi["011003020503"]; 
  mvaInputRefs["psi020503011003"                                     ]=&mvaPsi["020503011003"]; 
  mvaInputRefs["psi012003030503"                                     ]=&mvaPsi["012003030503"]; 
  mvaInputRefs["psi030503012003"                                     ]=&mvaPsi["030503012003"]; 
  mvaInputRefs["psi010503010502"                                     ]=&mvaPsi["010503010502"]; 
  mvaInputRefs["psi020503011002"                                     ]=&mvaPsi["020503011002"]; 
  mvaInputRefs["psi011003011002"                                     ]=&mvaPsi["011003011002"]; 
  mvaInputRefs["psi020504011002"                                     ]=&mvaPsi["020504011002"]; 
  mvaInputRefs["psi012004021003"                                     ]=&mvaPsi["012004021003"]; 
  mvaInputRefs["psi012004012002"                                     ]=&mvaPsi["012004012002"]; 
  mvaInputRefs["psi020504020503"                                     ]=&mvaPsi["020504020503"]; 
  mvaInputRefs["psi011004011002"                                     ]=&mvaPsi["011004011002"]; 
  mvaInputRefs["psi021004012002"                                     ]=&mvaPsi["021004012002"]; 
  mvaInputRefs["psi021003012002"                                     ]=&mvaPsi["021003012002"]; 
  mvaInputRefs["psi011004020503"                                     ]=&mvaPsi["011004020503"]; 
  mvaInputRefs["psi010504010502"                                     ]=&mvaPsi["010504010502"]; 
  mvaInputRefs["psi021004021003"                                     ]=&mvaPsi["021004021003"]; 
  mvaInputRefs["psi012003012002"                                     ]=&mvaPsi["012003012002"]; 
  mvaInputRefs["psi022004022003"                                     ]=&mvaPsi["022004022003"]; 
  TMVA::Reader *reader; TString theBDTWeights="";
  // to do - not hardcoded
  if(MVAVarType==2) {
    reader=new TMVA::Reader(); // =new TMVA::Reader();
    if(selection==kWHSR)        theBDTWeights = bdtWeightsResolved; //"weights/bdt_BDT_multiClass_resolved_jan10_I_test2.weights.xml";
    else if(selection==kWHFJSR) theBDTWeights = bdtWeightsBoosted; // "weights/bdt_BDT_multiClass_boosted_jan10_I_test2.weights.xml";
  } else if(MVAVarType==3) {
    reader=new TMVA::Reader(); // =new TMVA::Reader();
    if(selection==kWHSR)        theBDTWeights = bdtWeightsResolved; //"weights/bdt_BDT_multiClass_resolved_whAll_jan11_singleClass.weights.xml";
    //else if(selection==kWHFJSR) theBDTWeights = "weights/bdt_BDT_multiClass_boosted_jan10_I_test2.weights.xml";
  } 
  if(theBDTWeights!="") {
    vector<pair<string, string> > varLabelExprs=parseTmvaXmlFile(string(theBDTWeights.Data()));
    for(unsigned iVar=0; iVar<(unsigned)varLabelExprs.size(); iVar++)
      if(mvaInputRefs.find( varLabelExprs[iVar].second ) != mvaInputRefs.end())
        reader->AddVariable( varLabelExprs[iVar].second, mvaInputRefs[varLabelExprs[iVar].second]);
    reader->BookMVA("BDT", theBDTWeights.Data());
  }

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
    // physics calculations
    mT = TMath::Sqrt(2.*pfmet*lepton1Pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(lepton1Phi-pfmetphi))));
    if(selection>=kWHLightFlavorCR && selection<=kWHSR) {
      balanceVH = hbbDijetPt/topWBosonPt;  
      dPhiHbbJet1MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet1Phi - pfmetphi  ));
      dPhiHbbJet2MET = TMath::Abs(TVector2::Phi_mpi_pi( hbbJet2Phi - pfmetphi  ));
      nAddJetZeroOrOne = TMath::Min((int)(nJet-2), (int)1);  
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
    } else if (selection>=kWHLightFlavorFJCR && selection>=kWHFJSR) {
      dPhil1W   = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - topWBosonPhi));
      dPhil1fj1 = fabs(TVector2::Phi_mpi_pi(lepton1Phi   - fj1Phi  ));
      dPhiWfj1  = fabs(TVector2::Phi_mpi_pi(topWBosonPhi - fj1Phi  ));
      dEtal1fj1 = fabs(lepton1Eta   - fj1Eta );
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
    }
    float weight_jesUp=weight, weight_jesDown=weight;
    //jcuTotal->setJetPt(hbbJet1Pt); jcuTotal->setJetEta(hbbJet1Eta);
    //weight_jesUp   *= (1+jcuTotal->getUncertainty(true));
    //jcuTotal->setJetPt(hbbJet1Pt); jcuTotal->setJetEta(hbbJet1Eta);
    //weight_jesDown *= (1-jcuTotal->getUncertainty(true));
    //jcuTotal->setJetPt(hbbJet2Pt); jcuTotal->setJetEta(hbbJet2Eta);
    //weight_jesUp   *= (1+jcuTotal->getUncertainty(false));
    //jcuTotal->setJetPt(hbbJet2Pt); jcuTotal->setJetEta(hbbJet2Eta);
    //weight_jesDown *= (1-jcuTotal->getUncertainty(false));
    switch(MVAVarType) {
      case 1:
      default:
        if(selection<=kWHLightFlavorCR && selection<=kWHSR)
          MVAVar=hbbDijetPt;
        else if(selection<=kWHLightFlavorFJCR && selection<=kWHFJSR)
          MVAVar=fj1Pt;
        break;
      case 2:
        if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) 
          MVAVar=bDiscrMin;
        else if(selection==kWHSR)
          MVAVar=(reader->EvaluateMulticlass("BDT")[0]);
        else if(selection==kWHLightFlavorFJCR || selection==kWHHeavyFlavorFJCR || selection==kWH2TopFJCR) 
          MVAVar=fj1MSD;
        else if(selection==kWHFJSR)
          MVAVar=(reader->EvaluateMulticlass("BDT")[0]);
        break;
      case 3:
        if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR) 
          MVAVar=bDiscrMin;
        else if(selection==kWHSR)
          MVAVar=reader->EvaluateMVA("BDT");
        else if(selection==kWHLightFlavorFJCR || selection==kWHHeavyFlavorFJCR || selection==kWH2TopFJCR) 
          MVAVar=fj1MSD;
        else if(selection==kWHFJSR)
          MVAVar=reader->EvaluateMVA("BDT");
        break;
    }
    MVAVar=TMath::Min(MVAVar, (float)(MVAbins[MVAbins.size()-1]-0.001));

    float theVar; 
    bool passFullSel         = (selectionBits & selection)!=0; 
    bool passFullSel_jesUp   = (selectionBits_jesUp & selection)!=0; 
    bool passFullSel_jesDown = (selectionBits_jesDown & selection)!=0; 
    bool makePlot;
    // hack
    //if(typeLepSel==1 && theCategory==kPlotTop) weight*=-1.;
    if(selection>=kWHLightFlavorCR && selection<=kWHSR) for(int p=0; p<nPlots; p++) {
      makePlot=false;
      // Variables -- change the makePlot for n-1 later
      if      (histoNames[p]=="pTH"                 ) { theVar = TMath::Min(hbbDijetPt    , float(xmax[p]-0.001)); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="mH"                  ) { theVar = TMath::Min(hbbDijetMass  , float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="WpT"                 ) { theVar = TMath::Min(topWBosonPt   , float(xmax[p]-0.001)); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && theVar<=100.); }
      else if (histoNames[p]=="deltaPhiVH"          ) { theVar = TMath::Min(deltaPhiVH    , float(xmax[p]-0.001)); makePlot = passFullSel || (selection==kWHSR && (nMinusOneBits & selection)!=0 && theVar<2.5); }
      else if (histoNames[p]=="pTBalanceDijetW"     ) { theVar = TMath::Min(balanceVH     , float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="lepton1Pt"           ) { theVar = TMath::Min(lepton1Pt     , float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmet"               ) { theVar = TMath::Min(pfmet         , float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="pfmetsig"            ) { theVar = TMath::Min(pfmet         , float(xmax[p]-0.001)); makePlot = passFullSel || ((selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR) && (nMinusOneBits & selection)!=0 && theVar<=2.); }
      else if (histoNames[p]=="mTW"                 ) { theVar = TMath::Min(mT            , float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet1Pt"            ) { theVar = TMath::Min(hbbJet1Pt     , float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="Hbjet2Pt"            ) { theVar = TMath::Min(hbbJet2Pt     , float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet1MET"      ) { theVar = TMath::Min(dPhiHbbJet1MET, float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="dPhiHbbJet2MET"      ) { theVar = TMath::Min(dPhiHbbJet2MET, float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="bDiscrMin"           ) { theVar = TMath::Min(bDiscrMin     , float(xmax[p]-0.001)); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && (selection==kWHSR && theVar<bDiscrLoose )); }
      else if (histoNames[p]=="bDiscrMax"           ) { theVar = TMath::Min(bDiscrMax     , float(xmax[p]-0.001)); makePlot = passFullSel || ((nMinusOneBits & selection)!=0 && ((selection==kWHLightFlavorCR && (theVar<bDiscrLoose || theVar>bDiscrMedium) ) || (selection!=kWHLightFlavorCR && theVar<bDiscrTight)) ); }
      else if (histoNames[p]=="topMassLep1Met"      ) { theVar = TMath::Min(topMassLep1Met, float(xmax[p]-0.001)); makePlot = passFullSel; }
      else if (histoNames[p]=="sumEtSoft1"          ) { theVar = TMath::Min(sumEtSoft1    , float(xmax[p]-0.001)); makePlot = passFullSel; }
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
      else if (histoNames[p]=="MVAVar"              ) { theVar = MVAVar                                          ; makePlot = passFullSel; }
      else continue;
      if(!makePlot) continue;
      histos[p][theCategory]->Fill(theVar, weight);
    }
    // Fill systematic shapes
    if(passFullSel && theCategory!=kPlotData) {
      histo_pdfUp    [theCategory]->Fill(MVAVar, weight_pdfUp    ); 
      histo_pdfDown  [theCategory]->Fill(MVAVar, weight_pdfDown  ); 
      histo_QCDr1f2  [theCategory]->Fill(MVAVar, weight_QCDr1f2  ); 
      histo_QCDr1f5  [theCategory]->Fill(MVAVar, weight_QCDr1f5  ); 
      histo_QCDr2f1  [theCategory]->Fill(MVAVar, weight_QCDr2f1  ); 
      histo_QCDr2f2  [theCategory]->Fill(MVAVar, weight_QCDr2f2  ); 
      histo_QCDr5f1  [theCategory]->Fill(MVAVar, weight_QCDr5f1  ); 
      histo_QCDr5f5  [theCategory]->Fill(MVAVar, weight_QCDr5f5  ); 
      histo_muSFUp   [theCategory]->Fill(MVAVar, typeLepSel==1?weight_lepSFUp  : weight); 
      histo_eleSFUp  [theCategory]->Fill(MVAVar, typeLepSel==2?weight_lepSFUp : weight); 
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
      histo_jesUp   [theCategory]->Fill(MVAVar, weight); 
    if(passFullSel_jesDown && theCategory!=kPlotData)
      histo_jesDown [theCategory]->Fill(MVAVar, weight);
  }
  TFile *output_plots = new TFile(outputFileName,"RECREATE","",ROOT::CompressionSettings(ROOT::kZLIB,9));
  for(int p=0; p<nPlots; p++) {
    if(histoNames[p]=="") continue;
    TDirectory *plotDir = output_plots->mkdir(histoNames[p]);
    plotDir->cd();
    for(theCategory=kPlotData; theCategory!=nPlotCategories; theCategory++)
      histos[p][theCategory]->Write();

  }
  output_plots->Close();

  // Prepare the datacard C-('.' Q)
  char leptonChar='l';
  if(theLepSel==1) leptonChar='m';
  if(theLepSel==2) leptonChar='e';
  char shapeType[64];
  if     (MVAVarType==1) sprintf(shapeType,"ptShape");
  else if(MVAVarType==2) sprintf(shapeType,"multiClassBDTShape");
  else if(MVAVarType==3) sprintf(shapeType,"singleClassBDTShape");
  if(selection==kWHLightFlavorCR || selection==kWHHeavyFlavorCR || selection==kWH2TopCR || selection==kWHSR) { for(int nb=1; nb<=histos[pMVAVar][kPlotData]->GetNbinsX(); nb++) {
    float bgSum=0;
    for(int ic=0; ic<nPlotCategories; ic++)
      if(ic!=kPlotVH) bgSum+=histos[pMVAVar][ic]->GetBinContent(nb);
    if(bgSum<1e-3 && histos[pMVAVar][kPlotVH]->GetBinContent(nb)<1e-3) continue;
    char outputLimitsShape[512];
    sprintf(outputLimitsShape,"MitVHBBAnalysis/datacards/%s/histo_limits_W%cnHbb_%s_%s_bin%d.txt",dataCardDir.Data(),leptonChar,selectionNames[selection].Data(), shapeType, nb-1);
    ofstream card; card.open(outputLimitsShape);
    card << Form("imax 1 number of channels\n");
    card << Form("jmax * number of background\n");
    card << Form("kmax * number of nuisance parameters\n");
    card << Form("Observation %d\n", (int)round(isBlinded? bgSum : histos[pMVAVar][kPlotData]->GetBinContent(nb))); 
    card << Form("bin     "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%s ",Form("W%cn%sb%d",leptonChar,selectionNames[selection].Data(),nb-1)); card<<std::endl;
    card << Form("process "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%9s ",plotBaseNames[ic].Data()); card<<std::endl;
    card << Form("process "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%9d ", ic==kPlotVH?0:ic); card<<std::endl;
    card << Form("rate    "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%9.3f ", histos[pMVAVar][ic]->GetBinContent(nb)); card<<std::endl;
    // Systematics - Shape Uncertainties
    // PDF Acceptance
    card<<Form("pdf_qqbar_ACCEPT lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) 
      card << ((ic==kPlotVZbb||ic==kPlotVVLF||ic==kPlotWbb||ic==kPlotWb||ic==kPlotWLF||ic==kPlotZbb||ic==kPlotZb||ic==kPlotZLF||ic==kPlotVH)? Form("%7.5f/%7.5f ",
        histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,histo_pdfDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb)):1,
        histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Min(2.0,histo_pdfUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb)):1
      ): "- "); card<<std::endl;
    card<<Form("pdf_qg_ACCEPT lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) 
      card << ((ic==kPlotTop)? Form("%7.5f/%7.5f ",
        histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,histo_pdfDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb)):1,
        histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Min(2.0,histo_pdfUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb)):1
      ): "- "); card<<std::endl;
    card<<Form("pdf_gg_ACCEPT lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) 
      card << ((ic==kPlotTT)? Form("%7.5f/%7.5f ",
        histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,histo_pdfDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb)):1,
        histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Min(2.0,histo_pdfUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb)):1
      ): "- "); card<<std::endl;
    // QCD Scale
    for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
      if(ic==kPlotTop) continue; // Bug in MiniAOD
      card<<Form("QCDscale_%s lnN ",plotBaseNames[ic].Data());
      for(int jc=kPlotData+1; jc<nPlotCategories; jc++) {
        double scaleUp=1, scaleDown=1;
        if(jc==ic && histos[pMVAVar][ic]->GetBinContent(nb)>0) {
          scaleUp  =histos[pMVAVar][ic]->GetBinContent(nb);
          scaleDown=scaleUp;
          scaleUp=TMath::Max(scaleUp,histo_QCDr1f2[ic]->GetBinContent(nb)); scaleDown=TMath::Min(scaleDown,histo_QCDr1f2[ic]->GetBinContent(nb)); 
          scaleUp=TMath::Max(scaleUp,histo_QCDr1f5[ic]->GetBinContent(nb)); scaleDown=TMath::Min(scaleDown,histo_QCDr1f5[ic]->GetBinContent(nb));
          scaleUp=TMath::Max(scaleUp,histo_QCDr2f1[ic]->GetBinContent(nb)); scaleDown=TMath::Min(scaleDown,histo_QCDr2f1[ic]->GetBinContent(nb));
          scaleUp=TMath::Max(scaleUp,histo_QCDr2f2[ic]->GetBinContent(nb)); scaleDown=TMath::Min(scaleDown,histo_QCDr2f2[ic]->GetBinContent(nb));
          scaleUp=TMath::Max(scaleUp,histo_QCDr5f1[ic]->GetBinContent(nb)); scaleDown=TMath::Min(scaleDown,histo_QCDr5f1[ic]->GetBinContent(nb));
          scaleUp=TMath::Max(scaleUp,histo_QCDr5f5[ic]->GetBinContent(nb)); scaleDown=TMath::Min(scaleDown,histo_QCDr5f5[ic]->GetBinContent(nb));
          scaleUp  /=histos[pMVAVar][ic]->GetBinContent(nb);
          scaleDown/=histos[pMVAVar][ic]->GetBinContent(nb);
          scaleUp  =TMath::Min(scaleUp,2.0);
          scaleDown=TMath::Max(scaleDown,0.5);
        }
        card << (jc!=ic? "- ":Form("%7.5f/%7.5f ", scaleDown, scaleUp));
      }
      card<<std::endl;
    }
    // BTag Shape Uncertainty
    for(unsigned iPt=0; iPt<5; iPt++) for(unsigned iEta=0; iEta<3; iEta++) {
      card<<Form("CMS_bTagWeightJES_Pt%d_Eta%d      lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaJESDown     [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaJESUp     [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightLF_Pt%d_Eta%d       lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaLFDown      [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaLFUp      [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightHF_Pt%d_Eta%d       lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaHFDown      [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaHFUp      [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightHFStats1_Pt%d_Eta%d lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaHFStats1Down[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaHFStats1Up[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightHFStats2_Pt%d_Eta%d lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaHFStats2Down[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaHFStats2Up[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightLFStats1_Pt%d_Eta%d lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaLFStats1Down[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaLFStats1Up[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightLFStats2_Pt%d_Eta%d lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaLFStats2Down[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaLFStats2Up[iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightCErr1_Pt%d_Eta%d    lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaCErr1Down   [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaCErr1Up   [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
      card<<Form("CMS_bTagWeightCErr2_Pt%d_Eta%d    lnN ", iPt, iEta); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ", histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaCErr2Down   [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1, histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_cmvaCErr2Up   [iPt][iEta][ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1 ); card<<std::endl;
    }
    //// Total Jet Energy Scale Shape Uncertainty
    // TODO: check that the jesUp/Down is nonzero
    card<<Form("CMS_scale_j lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f/%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_jesDown[ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1,
      histos[pMVAVar][ic]->GetBinContent(nb)>0? TMath::Max(0.5,TMath::Min(2.0,histo_jesUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb))):1
    ); card<<std::endl;
    // Muon SF Shape
    card<<Form("CMS_eff2016_m lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_muSFUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
    ); card<<std::endl;
    // Electron SF Shape
    card<<Form("CMS_eff2016_e lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<Form("%7.5f ",
      histos[pMVAVar][ic]->GetBinContent(nb)>0? histo_eleSFUp  [ic]->GetBinContent(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1
    ); card<<std::endl;

    // Systematics -- Normalization Uncertainties
    // Muon Energy Scale Norm
    card<<Form("CMS_scale2016_m lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.01 "; card <<std::endl;
    // Electron Energy Scale Norm
    card<<Form("CMS_scale2016_e lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.01 "; card<<std::endl;
    // Trigger
    card<<Form("CMS_trigger2016 lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.01 "; card << std::endl;
    // JES (temporary)
    //card<<Form("CMS_scale_j lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.020 "; card << std::endl;
    // Lumi
    card<<Form("lumi_13TeV2016 lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< "1.025 "; card << std::endl;
    // VH Cross Section
    card<<Form("pdf_qqbar lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< ((ic==kPlotVZbb||ic==kPlotVVLF||ic==kPlotVH)?"1.01 ":"- "); card << std::endl;
    card<<Form("pdf_gg lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< ((ic==kPlotTop)?"1.01 ":"- "); card << std::endl;
    // Single Top Cross Section
    card<<Form("CMS_VH_TopNorm lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< (ic==kPlotTop? "1.30 ":"- "); card<<std::endl;
    // Diboson Cross Section
    card<<Form("CMS_VH_VVNorm lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<< ((ic==kPlotVZbb||ic==kPlotVVLF)? "1.30 ":"- "); card<<std::endl;
    // VH BDT
    if((selection==kWHSR || selection==kWHFJSR) && (MVAVarType==2 || MVAVarType==3)) {
      card<<Form("CMS_VH_BDT lnN "); for(int ic=kPlotData+1; ic<nPlotCategories; ic++) card<<"1.02 "; card<<std::endl;
    }
    // Statistical Uncertainties
    for(int ic=kPlotData+1; ic<nPlotCategories; ic++) {
      card<<Form("W%cn%s_%sStatBounding_13TeV2016_Bin%d lnN ",leptonChar,selectionNames[selection].Data(), plotBaseNames[ic].Data(), nb-1);
      for(int jc=kPlotData+1; jc<nPlotCategories; jc++) 
        card << (jc!=ic? "- ":Form("%7.5f ",histos[pMVAVar][ic]->GetBinContent(nb)>0? 1.+histos[pMVAVar][ic]->GetBinError(nb)/histos[pMVAVar][ic]->GetBinContent(nb):1));
      card<<std::endl;
    }
    if(selection>=kWHLightFlavorCR && selection<=kWHFJSR) {
      card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",plotBaseNames[kPlotTT ].Data(),plotBaseNames[kPlotTT ].Data()) << std::endl;
      card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",plotBaseNames[kPlotWbb].Data(),plotBaseNames[kPlotWbb].Data()) << std::endl;
      card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",plotBaseNames[kPlotWb ].Data(),plotBaseNames[kPlotWb ].Data()) << std::endl;
      card << Form("SF_%s_Wln rateParam  * %s 1 [0.2,5]",plotBaseNames[kPlotWLF].Data(),plotBaseNames[kPlotWLF].Data()) << std::endl;
    }
 
  }} 
}

