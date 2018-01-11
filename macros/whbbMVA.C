#include <TROOT.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TCut.h>
#include <TTree.h>
#include <TString.h>
#include "vhbbPlot.h"

void whbbMVA(
  TString plotTreeFileName, 
  TString extraString="", 
  bool isBoosted=false, 
  bool useGaussDeco=false,
  float eventFrac=1
) {
  gROOT->ProcessLine("TMVA::gConfig().GetVariablePlotting().fMaxNumOfAllowedVariablesForScatterPlots = 100");
  TFile *output_file;
  TMVA::Factory *factory;
  
  // Determine the input trees
  TFile *plotTreeFile = TFile::Open(plotTreeFileName, "READ"); assert(plotTreeFile);
  TTree *plotTree = (TTree*)plotTreeFile->Get("plotTree"); assert(plotTree);
  
  // Initialize the factory
  TString trainName="BDT_multiClass_"+TString(isBoosted?"boosted":"resolved")+(extraString == "" ? "" : "_"+extraString);
  output_file=TFile::Open("/data/t3home000/dhsu/mva/whbb/trainingResult_whbb_"+trainName+".root", "RECREATE");
  //factory = new TMVA::Factory("bdt", output_file, "!V:!Silent:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Multiclass");
  TString factoryOptions="!V:!Silent:DrawProgressBar:AnalysisType=Multiclass";
  if(useGaussDeco) factoryOptions += ":Transformations=G,D";
  else             factoryOptions += ":Transformations=I";
  factory = new TMVA::Factory("bdt", output_file, factoryOptions);
  
  factory->AddTree(plotTree,"Signal",1.0, "theCategory==12", "test train");
  factory->AddTree(plotTree,"W+jets",1.0, "theCategory==6 || theCategory==7 || theCategory==8", "test train");
  factory->AddTree(plotTree,"N-top" ,1.0, "theCategory==4 || theCategory==5 ", "test train");
  factory->AddTree(plotTree,"Other" ,1.0, "theCategory==2 || theCategory==3 || theCategory==9 || theCategory==10 || theCategory==11", "test train"); // no QCD
  factory->SetWeightExpression("weight", "Signal");
  factory->SetWeightExpression("weight", "W+jets");
  factory->SetWeightExpression("weight", "N-top" );
  factory->SetWeightExpression("weight", "Other" );
  
  TCut preselectionCut;
  
  if(isBoosted) {
    preselectionCut = Form("(selectionBits & %d)!=0 && weight<500. && fj1ECFN_2_4_20>0", vhbbPlot::kWHFJSR);
    
    TString dPhil1W   = "dPhil1W   := fabs(TVector2::Phi_mpi_pi(lepton1Phi   - topWBosonPhi))";
    TString dPhil1fj1 = "dPhil1fj1 := fabs(TVector2::Phi_mpi_pi(lepton1Phi   - fj1Phi  ))";
    TString dPhiWfj1  = "dPhiWfj1  := fabs(TVector2::Phi_mpi_pi(topWBosonPhi - fj1Phi  ))";
    TString dEtal1fj1 = "dEtal1fj1 := fabs(lepton1Eta   - fj1Eta )";
    TString greaterSubjetCSV = "greaterSubjetCSV := fj1MaxCSV!=fj1MaxCSV? 0. : TMath::Min(1., TMath::Max(0., fj1MaxCSV))";
    TString lesserSubjetCSV  = "lesserSubjetCSV  := fj1MinCSV!=fj1MinCSV? 0. : TMath::Min(1., TMath::Max(0., fj1MinCSV))";

    //ECF definitions
    std::map<string, TString> psi, psiName;
    psi["021004010502"] = "psi021004010502 := fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)";
    psi["012004010502"] = "psi012004010502 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)";
    psi["021003011002"] = "psi021003011002 := fj1ECFN_2_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)";
    psi["022004011002"] = "psi022004011002 := fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),4.00)";
    psi["020503010502"] = "psi020503010502 := fj1ECFN_2_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)";
    psi["022003012002"] = "psi022003012002 := fj1ECFN_2_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),2.00)";
    psi["012004020503"] = "psi012004020503 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_05),2.00)";
    psi["021004020503"] = "psi021004020503 := fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_2_3_05),2.00)";
    psi["032003012002"] = "psi032003012002 := fj1ECFN_3_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),3.00)";
    psi["012003011002"] = "psi012003011002 := fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)";
    psi["011003010502"] = "psi011003010502 := fj1ECFN_1_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)";
    psi["031003011002"] = "psi031003011002 := fj1ECFN_3_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),3.00)";
    psi["031003012002"] = "psi031003012002 := fj1ECFN_3_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.50)";
    psi["030503011002"] = "psi030503011002 := fj1ECFN_3_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.50)";
    psi["022004012002"] = "psi022004012002 := fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),2.00)";
    psi["021004011002"] = "psi021004011002 := fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)";
    psi["012004011002"] = "psi012004011002 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),2.00)";
    psi["022004021003"] = "psi022004021003 := fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_10),2.00)";
    psi["012003010502"] = "psi012003010502 := fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)";
    psi["022003011002"] = "psi022003011002 := fj1ECFN_2_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_10),4.00)";
    psi["022004031003"] = "psi022004031003 := fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_3_3_10),1.33)";
    psi["012004030503"] = "psi012004030503 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_3_3_05),1.33)";
    psi["021004030503"] = "psi021004030503 := fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_3_3_05),1.33)";
    psi["020504010502"] = "psi020504010502 := fj1ECFN_2_4_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)";
    psi["022004030503"] = "psi022004030503 := fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_3_3_05),2.67)";
    psi["011004010502"] = "psi011004010502 := fj1ECFN_1_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),2.00)";
    psi["030503010502"] = "psi030503010502 := fj1ECFN_3_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),3.00)";
    psi["021003010502"] = "psi021003010502 := fj1ECFN_2_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_05),4.00)";
    psi["012004011003"] = "psi012004011003 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_3_10),2.00)";
    psi["012004021004"] = "psi012004021004 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_2_4_10),1.00)";
    psi["021004012004"] = "psi021004012004 := fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_4_20),1.00)";
    psi["011003020503"] = "psi011003020503 := fj1ECFN_1_3_10/pow(TMath::Max(0.,fj1ECFN_2_3_05),1.00)";
    psi["020503011003"] = "psi020503011003 := fj1ECFN_2_3_05/pow(TMath::Max(0.,fj1ECFN_1_3_10),1.00)";
    psi["012003030503"] = "psi012003030503 := fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_3_3_05),1.33)";
    psi["030503012003"] = "psi030503012003 := fj1ECFN_3_3_05/pow(TMath::Max(0.,fj1ECFN_1_3_20),0.75)";
    psi["010503010502"] = "psi010503010502 := fj1ECFN_1_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),1.00)";
    psi["020503011002"] = "psi020503011002 := fj1ECFN_2_3_05/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)";
    psi["011003011002"] = "psi011003011002 := fj1ECFN_1_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)";
    psi["020504011002"] = "psi020504011002 := fj1ECFN_2_4_05/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)";
    psi["012004021003"] = "psi012004021003 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_10),1.00)";
    psi["012004012002"] = "psi012004012002 := fj1ECFN_1_4_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)";
    psi["020504020503"] = "psi020504020503 := fj1ECFN_2_4_05/pow(TMath::Max(0.,fj1ECFN_2_3_05),1.00)";
    psi["011004011002"] = "psi011004011002 := fj1ECFN_1_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_10),1.00)";
    psi["021004012002"] = "psi021004012002 := fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)";
    psi["021003012002"] = "psi021003012002 := fj1ECFN_2_3_10/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)";
    psi["011004020503"] = "psi011004020503 := fj1ECFN_1_4_10/pow(TMath::Max(0.,fj1ECFN_2_3_05),1.00)";
    psi["010504010502"] = "psi010504010502 := fj1ECFN_1_4_05/pow(TMath::Max(0.,fj1ECFN_1_2_05),1.00)";
    psi["021004021003"] = "psi021004021003 := fj1ECFN_2_4_10/pow(TMath::Max(0.,fj1ECFN_2_3_10),1.00)";
    psi["012003012002"] = "psi012003012002 := fj1ECFN_1_3_20/pow(TMath::Max(0.,fj1ECFN_1_2_20),1.00)";
    psi["022004022003"] = "psi022004022003 := fj1ECFN_2_4_20/pow(TMath::Max(0.,fj1ECFN_2_3_20),1.00)";
    
    factory->AddVariable( "nIsojet"                           , "N central AK4 jets"                , ""           , 'I');
    factory->AddVariable( "fj1Pt"                             , "FJ p_{T}"                          , "GeV"        , 'F');
    factory->AddVariable( "fj1Eta"                            , "FJ #eta"                           , ""           , 'F');
    factory->AddVariable( greaterSubjetCSV                    , "FJ greater subjet B-tag"           , ""           , 'F');
    factory->AddVariable( lesserSubjetCSV                     , "FJ lesser subjet B-tag"            , ""           , 'F');
    factory->AddVariable( "fj1DoubleCSV"                      , "FJ double B-tag"                   , ""           , 'F');
    factory->AddVariable( "fj1HTTMass"                        , "FJ HTT mass"                       , "GeV"        , 'F');
    factory->AddVariable( "fj1HTTFRec"                        , "FJ HTT f_{rec}"                    , ""           , 'F');
    factory->AddVariable( "fj1MSD"                            , "FJ soft drop mass"                 , "GeV"        , 'F');
    factory->AddVariable( "fj1Tau32"                          , "FJ #tau_{32}"                      , ""           , 'F');
    factory->AddVariable( "fj1Tau32SD"                        , "FJ #tau_{32}^{SD}"                 , ""           , 'F');
    factory->AddVariable( "fj1Tau21"                          , "FJ #tau_{21}"                      , ""           , 'F');
    factory->AddVariable( "fj1Tau21SD"                        , "FJ #tau_{21}^{SD}"                 , ""           , 'F');
    factory->AddVariable( "lepton1Pt"                         , "Lepton p_{T}"                      , "GeV"        , 'F');
    factory->AddVariable( "deltaPhiVH"                        , "#Delta#phi(W,H)"                   , "Rad"        , 'F');
    factory->AddVariable( "deltaPhiLep1Met"                   , "#Delta#phi(lepton,p_{T}^{miss})"   , "Rad"        , 'F');
    factory->AddVariable( "topWBosonPt"                       , "W boson p_{T}"                     , "GeV"        , 'F');
    factory->AddVariable( "mT"                                , "W boson m_{T}"                     , "GeV"        , 'F');
    factory->AddVariable( "pfmetsig"                          , "E_{T}^{miss} sig."                 , "#sqrt{GeV}" , 'F');
    factory->AddVariable( "pfmet"                             , "p_{T}^{miss}"                      , "GeV"        , 'F');
    factory->AddVariable( dPhil1W                             , "#Delta#phi(lep,W)"                 , "Rad"        , 'F');
    factory->AddVariable( dPhil1fj1                           , "#Delta#phi(lep,FJ)"                , "Rad"        , 'F');
    factory->AddVariable( dPhiWfj1                            , "#Delta#phi(W,FJ)"                  , "Rad"        , 'F');
    factory->AddVariable( dEtal1fj1                           , "|#Delta#eta(lep,FJ)|"              , ""           , 'F');
    factory->AddVariable( psi["021004010502"]                 , "#psi(2,1.0,4,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004010502"]                 , "#psi(1,2.0,4,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["021003011002"]                 , "#psi(2,1.0,3,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022004011002"]                 , "#psi(2,2.0,4,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["020503010502"]                 , "#psi(2,0.5,3,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022003012002"]                 , "#psi(2,2.0,3,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004020503"]                 , "#psi(1,2.0,4,2,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["021004020503"]                 , "#psi(2,1.0,4,2,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["032003012002"]                 , "#psi(3,2.0,3,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["012003011002"]                 , "#psi(1,2.0,3,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["011003010502"]                 , "#psi(1,1.0,3,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["031003011002"]                 , "#psi(3,1.0,3,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["031003012002"]                 , "#psi(3,1.0,3,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["030503011002"]                 , "#psi(3,0.5,3,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022004012002"]                 , "#psi(2,2.0,4,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["021004011002"]                 , "#psi(2,1.0,4,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004011002"]                 , "#psi(1,2.0,4,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022004021003"]                 , "#psi(2,2.0,4,2,1.0,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["012003010502"]                 , "#psi(1,2.0,3,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022003011002"]                 , "#psi(2,2.0,3,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022004031003"]                 , "#psi(2,2.0,4,3,1.0,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004030503"]                 , "#psi(1,2.0,4,3,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["021004030503"]                 , "#psi(2,1.0,4,3,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["020504010502"]                 , "#psi(2,0.5,4,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022004030503"]                 , "#psi(2,2.0,4,3,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["011004010502"]                 , "#psi(1,1.0,4,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["030503010502"]                 , "#psi(3,0.5,3,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["021003010502"]                 , "#psi(2,1.0,3,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004011003"]                 , "#psi(1,2.0,4,1,1.0,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004021004"]                 , "#psi(1,2.0,4,2,1.0,4)"             , ""           , 'F'); 
    factory->AddVariable( psi["021004012004"]                 , "#psi(2,1.0,4,1,2.0,4)"             , ""           , 'F'); 
    factory->AddVariable( psi["011003020503"]                 , "#psi(1,1.0,3,2,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["020503011003"]                 , "#psi(2,0.5,3,1,1.0,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["012003030503"]                 , "#psi(1,2.0,3,3,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["030503012003"]                 , "#psi(3,0.5,3,1,2.0,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["010503010502"]                 , "#psi(1,0.5,3,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["020503011002"]                 , "#psi(2,0.5,3,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["011003011002"]                 , "#psi(1,1.0,3,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["020504011002"]                 , "#psi(2,0.5,4,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004021003"]                 , "#psi(1,2.0,4,2,1.0,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["012004012002"]                 , "#psi(1,2.0,4,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["020504020503"]                 , "#psi(2,0.5,4,2,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["011004011002"]                 , "#psi(1,1.0,4,1,1.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["021004012002"]                 , "#psi(2,1.0,4,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["021003012002"]                 , "#psi(2,1.0,3,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["011004020503"]                 , "#psi(1,1.0,4,2,0.5,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["010504010502"]                 , "#psi(1,0.5,4,1,0.5,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["021004021003"]                 , "#psi(2,1.0,4,2,1.0,3)"             , ""           , 'F'); 
    factory->AddVariable( psi["012003012002"]                 , "#psi(1,2.0,3,1,2.0,2)"             , ""           , 'F'); 
    factory->AddVariable( psi["022004022003"]                 , "#psi(2,2.0,4,2,2.0,3)"             , ""           , 'F'); 
    //factory->AddVariable( "fj1SDEFrac100"                     , "FJ SD E frac. hardest 100"         , ""           , 'F');
  } else {
    preselectionCut = Form("(selectionBits & %d)!=0 && weight<500.", vhbbPlot::kWHSR);
    
    TString dPhil1W   = "dPhil1W  := fabs(TVector2::Phi_mpi_pi(lepton1Phi   - topWBosonPhi))";
    TString dPhil1b1  = "dPhil1b1 := fabs(TVector2::Phi_mpi_pi(lepton1Phi   - hbbJet1Phi  ))";
    TString dPhil1b2  = "dPhil1b2 := fabs(TVector2::Phi_mpi_pi(lepton1Phi   - hbbJet2Phi  ))";
    TString dPhiWb1   = "dPhiWb1  := fabs(TVector2::Phi_mpi_pi(topWBosonPhi - hbbJet1Phi  ))";
    TString dPhiWb2   = "dPhiWb2  := fabs(TVector2::Phi_mpi_pi(topWBosonPhi - hbbJet2Phi  ))";
    TString dPhib1b2  = "dPhib1b2 := fabs(TVector2::Phi_mpi_pi(hbbJet1Phi   - hbbJet2Phi  ))";
    TString dEtal1W   = "dEtal1W  := fabs(lepton1Eta   - topWBosonEta)";
    TString dEtal1b1  = "dEtal1b1 := fabs(lepton1Eta   - hbbJet1Eta  )";
    TString dEtal1b2  = "dEtal1b2 := fabs(lepton1Eta   - hbbJet2Eta  )";
    TString dEtaWb1   = "dEtaWb1  := fabs(topWBosonEta - hbbJet1Eta  )";
    TString dEtaWb2   = "dEtaWb2  := fabs(topWBosonEta - hbbJet2Eta  )";
    TString dEtab1b2  = "dEtab1b2 := fabs(hbbJet1Eta   - hbbJet2Eta  )";

    factory->AddVariable( "hbbDijetMass"                      , "H(bb) mass"                        , "GeV"        , 'F');
    factory->AddVariable( "hbbDijetPt"                        , "H(bb) p_{T}"                       , "GeV"        , 'F');
    factory->AddVariable( "topWBosonPt"                       , "W boson p_{T}"                     , "GeV"        , 'F');
    factory->AddVariable( "bDiscrMin"                         , "Lesser b-tag"                      , ""           , 'F');
    factory->AddVariable( "topMassLep1Met"                    , "Top mass"                          , "GeV"        , 'F');
    factory->AddVariable( "deltaPhiVH"                        , "#Delta#phi(W,H)"                   , "Rad"        , 'F');
    factory->AddVariable( "nAJ := TMath::Min(nJet-2,1)"       , "Additional jets"                   , ""           , 'I');
    factory->AddVariable( "deltaPhiLep1Met"                   , "#Delta#phi(lepton,p_{T}^{miss})"   , "Rad"        , 'F');
    factory->AddVariable( "mT"                                , "W boson m_{T}"                     , "GeV"        , 'F');
    factory->AddVariable( "pfmet"                             , "p_{T}^{miss}"                      , "GeV"        , 'F');
    factory->AddVariable( "nSoft5"                            , "N^{soft}_{5}"                      , ""           , 'I');
    // unorthodox variables
    factory->AddVariable( "lepton1Pt"                         , "Lepton p_{T}"                      , "GeV"        , 'F');
    factory->AddVariable( "bDiscrMax"                         , "Greater b-tag"                     , ""           , 'F');
    factory->AddVariable( "pfmetsig"                          , "E_{T}^{miss} sig."                 , "#sqrt{GeV}" , 'F');
    factory->AddVariable( "sumEtSoft1"                        , "#sum E_{T}(soft 1)"                , "GeV"        , 'F');
    factory->AddVariable( "nSoft2"                            , "N^{soft}_{2}"                      , ""           , 'I');
    factory->AddVariable( "nSoft10"                           , "N^{soft}_{10}"                     , ""           , 'I');
    factory->AddVariable( "topWBosonCosThetaCS"               , "W boson cos #theta^{CS}"           , ""           , 'F');
    factory->AddVariable( "hbbCosThetaJJ"                     , "cos #theta(bb)"                    , ""           , 'F');
    factory->AddVariable( "hbbCosThetaCSJ1"                   , "cos #theta^{CS} harder b"          , ""           , 'F');
    factory->AddVariable( "hbbJet1Pt"                         , "p_{T}(b1)"                         , "GeV"        , 'F');
    factory->AddVariable( "hbbJet2Pt"                         , "p_{T}(b2)"                         , "GeV"        , 'F');
    factory->AddVariable( dPhil1W                             , "#Delta#phi(lep,W)"                 , "Rad"        , 'F');
    factory->AddVariable( dPhil1b1                            , "#Delta#phi(lep,b1)"                , "Rad"        , 'F');
    factory->AddVariable( dPhil1b2                            , "#Delta#phi(lep,b2)"                , "Rad"        , 'F');
    factory->AddVariable( dPhiWb1                             , "#Delta#phi(W,b1)"                  , "Rad"        , 'F');
    factory->AddVariable( dPhiWb2                             , "#Delta#phi(W,b2)"                  , "Rad"        , 'F');
    factory->AddVariable( dPhib1b2                            , "#Delta#phi(b1,b2)"                 , "Rad"        , 'F');
    factory->AddVariable( dEtal1W                             , "|#Delta#eta(lep,W)|"               , ""           , 'F');
    factory->AddVariable( dEtal1b1                            , "|#Delta#eta(lep,b1)|"              , ""           , 'F');
    factory->AddVariable( dEtal1b2                            , "|#Delta#eta(lep,b2)|"              , ""           , 'F');
    factory->AddVariable( dEtaWb1                             , "|#Delta#eta(W,b1)|"                , ""           , 'F');
    factory->AddVariable( dEtaWb2                             , "|#Delta#eta(W,b2)|"                , ""           , 'F');
    factory->AddVariable( dEtab1b2                            , "|#Delta#eta(b1,b2)|"               , ""           , 'F');
  }
  TString prepareOptions="";
  if(eventFrac<1) {
    unsigned nTrain_Signal = round(eventFrac/2.*plotTree->GetEntries(preselectionCut && "theCategory==12"                                   ));
    unsigned nTrain_Wjets  = round(eventFrac/2.*plotTree->GetEntries(preselectionCut && "theCategory==6 || theCategory==7 || theCategory==8"));
    unsigned nTrain_Ntop   = round(eventFrac/2.*plotTree->GetEntries(preselectionCut && "theCategory==4 || theCategory==5"                  ));
    unsigned nTrain_Other  = round(eventFrac/2.*plotTree->GetEntries(preselectionCut && "theCategory==2 || theCategory==3 || theCategory==9 || theCategory==10 || theCategory==11"));
    prepareOptions=Form("NormMode=None:nTrain_Signal=%d:nTest_Signal=%d:nTrain_W+jets=%d:nTest_W+jets=%d:nTrain_N-top=%d:nTest_N-top=%d:nTrain_Other=%d:nTest_Other=%d",
      nTrain_Signal,nTrain_Signal,
      nTrain_Wjets,nTrain_Wjets,
      nTrain_Ntop,nTrain_Ntop,
      nTrain_Other,nTrain_Other);
  } else {
    prepareOptions=Form("NormMode=None");
  }
  factory->PrepareTrainingAndTestTree(preselectionCut, prepareOptions);
  TString hyperparameters="!H:!V:NTrees=1000:MinNodeSize=5%:MaxDepth=4:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=1000:PruneMethod=NoPruning";
  if(useGaussDeco) hyperparameters += ":VarTransform=G,D";
  factory->BookMethod(
    TMVA::Types::kBDT,
    trainName,
    hyperparameters
  );
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
}
