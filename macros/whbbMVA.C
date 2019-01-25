#include <TROOT.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TCut.h>
#include <TTree.h>
#include <TString.h>
#include "vhbbPlot.h"

using namespace vhbbPlot;
void whbbMVA(
  TString inputFileName, 
  TString extraString="", 
  bool isBoosted=false, 
  bool useGaussDeco=false,
  bool useMulticlass=false,
  bool vzbbMode=false
) {
  gROOT->ProcessLine("TMVA::gConfig().GetVariablePlotting().fMaxNumOfAllowedVariablesForScatterPlots = 50");
  TFile *output_file;
  TMVA::Factory *factory;
  
  // Determine the input trees
  TFile *inputFile = TFile::Open(inputFileName,"READ");
  TTree *mvaTree = (TTree*)inputFile->Get("mvaTree");
  
  // Initialize the factory
  TString trainName="BDT_"+TString(useMulticlass?"multiClass_":"singleClass_")+TString(isBoosted?"boosted":"resolved")+(extraString == "" ? "" : "_"+extraString);
  output_file=TFile::Open(trainName+".root", "RECREATE");
  //factory = new TMVA::Factory("bdt", output_file, "!V:!Silent:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Multiclass");
  TString factoryOptions="!V:!Silent:DrawProgressBar";
  //TString factoryOptions="!V:!Silent:!DrawProgressBar";
  if(useMulticlass) factoryOptions+=":AnalysisType=Multiclass";
  else              factoryOptions+=":AnalysisType=Classification";
  if(useGaussDeco)  factoryOptions += ":Transformations=G,D";

  int signalIdx = vzbbMode? (int)vhbbPlot::kPlotVZbb:(int)vhbbPlot::kPlotWH;
  //else              factoryOptions += ":Transformations=I";
  factory = new TMVA::Factory("bdt", output_file, factoryOptions);
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("MitVHBBAnalysis");
  if(useMulticlass) {
    return;
  } else {
    TCut cutTrainSignal = Form("%s && (category==%d)",trainTreeEventSplitStr.Data(),signalIdx);
    TCut cutTrainBkg    = Form("%s && (category!=%d)",trainTreeEventSplitStr.Data(),signalIdx);
    TCut cutTestSignal  = Form("%s && (category==%d)",testTreeEventSplitStr.Data() ,signalIdx);
    TCut cutTestBkg     = Form("%s && (category!=%d)",testTreeEventSplitStr.Data() ,signalIdx);
    if(isBoosted) {
      cutTrainBkg     = cutTrainBkg   && "psi022004022003>=0";
      cutTrainSignal  = cutTrainSignal&& "psi022004022003>=0";
      cutTestBkg      = cutTestBkg    && "psi022004022003>=0";
      cutTestSignal   = cutTestSignal && "psi022004022003>=0";
    }
    dataloader->AddTree(mvaTree, "Background", 1.0, cutTrainBkg   , "train");
    dataloader->AddTree(mvaTree, "Signal"    , 1.0, cutTrainSignal, "train");
    dataloader->AddTree(mvaTree, "Background", 1.0, cutTestBkg   , "test");
    dataloader->AddTree(mvaTree, "Signal"    , 1.0, cutTestSignal, "test");
    dataloader->SetWeightExpression("abs(weight)", "Signal");
    dataloader->SetWeightExpression("abs(weight)", "Background");
  } 
  
  TCut preselectionCut;
  if(isBoosted) {
    dataloader->AddVariable("lepton1Pt"        ,  "lepton pT"                 , "", 'F' ); 
    dataloader->AddVariable("lepton1Charge"    ,  "lepton charge"             , "", 'F' ); 
    dataloader->AddVariable("nIsojet"          ,  "N iso. AK4 jets"           , "", 'F' ); 
    dataloader->AddVariable("MSD"              ,  "mSD"                       , "", 'F' ); 
    dataloader->AddVariable("Tau21SD"          ,  "#tau2/#tau1 SD"            , "", 'F' ); 
    dataloader->AddVariable("Tau32SD"          ,  "#tau3/#tau2 SD"            , "", 'F' ); 
    dataloader->AddVariable("fjPt"             ,  "Fatjet pT"                 , "", 'F' ); 
    dataloader->AddVariable("psi022004031003"  ,  "#psi(2,2.0,4,3,1.0,3)"     , "", 'F' ); 
    dataloader->AddVariable("psi022004030503"  ,  "#psi(2,2.0,4,3,0.5,3)"     , "", 'F' ); 
    dataloader->AddVariable("ptBalanceWHFJ"    ,  "WH pT balance"             , "", 'F' ); 
    dataloader->AddVariable("dEtaLep1FJ"       ,  "#Delta#eta(lepton,FJ)"     , "", 'F' ); 
    dataloader->AddVariable("dPhiWHFJ"         ,  "#Delta#phi(W,H)"           , "", 'F' ); 
    dataloader->AddVariable("HTTFRec"          ,  "HepTopTagger f_rec"        , "", 'F' ); 
    dataloader->AddVariable("doubleBTag"       ,  "Double B Tag"              , "", 'F' ); 
    //dataloader->AddVariable("WBosonPt"         ,  "W boson pT"                , "", 'F' ); 
    //dataloader->AddVariable("psi022004022003"  ,  "#psi(2,2.0,4,2,2.0,3)"     , "", 'F' ); 
  } else {
    dataloader->AddVariable("lepton1Pt"        ,  "lepton pT"                 , "", 'F' ); 
    dataloader->AddVariable("lepton1Charge"    , "lepton1Charge"              , "", 'F' );
    dataloader->AddVariable("sumEtSoft1"       , "sumEtSoft1"                 , "", 'F' );
    dataloader->AddVariable("nSoft5"           , "nSoft5"                     , "", 'F' );
    dataloader->AddVariable("bjet1Pt"          , "bjet1Pt"                    , "", 'F' );
    dataloader->AddVariable("bjet2Pt"          , "bjet2Pt"                    , "", 'F' );
    dataloader->AddVariable("bjet2btag"        , "bjet2btag"                  , "", 'F' );
    dataloader->AddVariable("hbbpt"            , "hbbpt"                      , "", 'F' );
    dataloader->AddVariable("hbbm"             , "hbbm"                       , "", 'F' );
    dataloader->AddVariable("dPhiWH"           , "dPhiWH"                     , "", 'F' );
    dataloader->AddVariable("ptBalanceWH"      , "ptBalanceWH"                , "", 'F' );
    dataloader->AddVariable("topMass"          , "topMass"                    , "", 'F' );
    dataloader->AddVariable("dEtaLep1H"        , "dEtaLep1H"                  , "", 'F' );
    dataloader->AddVariable("nAddJet"          , "nAddJet"                    , "", 'F' );
    //dataloader->AddVariable("WBosonPt"         , "WBosonPt"                   , "", 'F' );
    //dataloader->AddVariable("dRBjets"          , "dRBjets"                    , "", 'F' );
  }
  TString prepareOptions="NormMode=None";
    prepareOptions+=":SplitMode=Block"; // use e.g. all events selected by trainTreeEventSplitStr for training
    prepareOptions+=":MixMode=Random";
  dataloader->PrepareTrainingAndTestTree("", prepareOptions);
  
  // for resolved
  TString hyperparameters=
  isBoosted?
  "!H:!V:BoostType=AdaBoost:MinNodeSize=5%:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=MisClassificationError:NTrees=500:MaxDepth=3:AdaBoostBeta=0.12:nCuts=10000":
  "!H:!V:BoostType=AdaBoost:MinNodeSize=5%:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=MisClassificationError:NTrees=500:MaxDepth=3:AdaBoostBeta=0.12:nCuts=10000";

  //TString hyperparameters="!H:!V:NTrees=500:MinNodeSize=5%:MaxDepth=3:BoostType=Grad:Shrinkage=0.1:nCuts=30:PruneMethod=CostComplexity";
  //TString hyperparameters="!H:!V:NTrees=500:NegWeightTreatment=Pray:MinNodeSize=5%:MaxDepth=2:BoostType=Grad:Shrinkage=0.1:nCuts=30";
  // for boosted
  //TString hyperparameters="!H:!V:NTrees=1000:NegWeightTreatment=Pray:SeparationType=MisClassificationError:MinNodeSize=5%:MaxDepth=2:BoostType=Grad:Shrinkage=0.05:nCuts=1000";
  //TString hyperparameters="!H:!V:NTrees=1000:NegWeightTreatment=Pray:SeparationType=MisClassificationError:MinNodeSize=5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.12:nCuts=1000";
  //if(useGaussDeco) hyperparameters += ":VarTransform=G,D";
  factory->BookMethod(dataloader, TMVA::Types::kBDT, trainName, hyperparameters);
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  output_file->Close();
}
