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
void zhbbMVA(
  TString inputFileName, 
  TString extraString="", 
  bool isBoosted=false, 
  bool useGaussDeco=false,
  bool useMulticlass=true,
  TString variableSet = "orthodox" // orthodox | fullResolved | everything
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
  //else              factoryOptions += ":Transformations=I";
  factory = new TMVA::Factory("bdt", output_file, factoryOptions);
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("MitVHBBAnalysis");
  if(useMulticlass) {
    return;
  } else {
    TCut cutTrainSignal = Form("%s && (category==%d || category==%d)",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    TCut cutTrainBkg    = Form("%s && (category!=%d && category!=%d)",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    TCut cutTestSignal = Form("%s && (category==%d || category==%d)",testTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    TCut cutTestBkg    = Form("%s && (category!=%d && category!=%d)",testTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    dataloader->AddTree(mvaTree, "Background", 1.0, cutTrainBkg   , "train");
    dataloader->AddTree(mvaTree, "Signal"    , 1.0, cutTrainSignal, "train");
    dataloader->AddTree(mvaTree, "Background", 1.0, cutTestBkg   , "test");
    dataloader->AddTree(mvaTree, "Signal"    , 1.0, cutTestSignal, "test");
    dataloader->SetWeightExpression("weight", "Signal");
    dataloader->SetWeightExpression("weight", "Background");
  } 
  
  TCut preselectionCut;
  if(isBoosted) {
    dataloader->AddVariable("nIsojet"       , "nIsojet"       , "", 'F');
    dataloader->AddVariable("fjPt"          , "fjPt"          , "", 'F');
    dataloader->AddVariable("MSD"           , "MSD"           , "", 'F');
    dataloader->AddVariable("Tau21SD"       , "Tau21SD"       , "", 'F');
    dataloader->AddVariable("ptBalanceZHFJ" , "ptBalanceZHFJ" , "", 'F');
    dataloader->AddVariable("dEtaZHFJ"      , "dEtaZHFJ"      , "", 'F');
    dataloader->AddVariable("dPhiZHFJ"      , "dPhiZHFJ"      , "", 'F');
    dataloader->AddVariable("mTZHFJ"        , "mTZHFJ"        , "", 'F');
    dataloader->AddVariable("ptBalanceL1L2" , "ptBalanceL1L2" , "", 'F');
    dataloader->AddVariable("dRL1L2"        , "dRL1L2"        , "", 'F');
    dataloader->AddVariable("ZBosonPt"      , "ZBosonPt"      , "", 'F'); 
    dataloader->AddVariable("lepton1Pt"     , "lepton1Pt"     , "", 'F');
    dataloader->AddVariable("lepton2Pt"     , "lepton2Pt"     , "", 'F');
    //dataloader->AddVariable("lepton1Eta"    , "lepton1Eta"    , "", 'F');
    //dataloader->AddVariable("lepton2Eta"    , "lepton2Eta"    , "", 'F');
    dataloader->AddVariable("deltaM"        , "deltaM"        , "", 'F');
    dataloader->AddVariable("CosThetaCS"  ,"CosThetaCS"  , "", 'F'); 
  } else {
    dataloader->AddVariable("sumEtSoft1"  ,"sumEtSoft1"  , "", 'F'); 
    //dataloader->AddVariable("nSoft2"      ,"nSoft2"      , "", 'F'); 
    //dataloader->AddVariable("nSoft5"      ,"nSoft5"      , "", 'F'); 
    //dataloader->AddVariable("nSoft10"     ,"nSoft10"     , "", 'F'); 
    dataloader->AddVariable("bjet1Pt"     ,"bjet1Pt"     , "", 'F'); 
    dataloader->AddVariable("bjet2Pt"     ,"bjet2Pt"     , "", 'F'); 
    dataloader->AddVariable("bjet1btag"   ,"bjet1btag"   , "", 'F'); 
    dataloader->AddVariable("bjet2btag"   ,"bjet2btag"   , "", 'F'); 
    dataloader->AddVariable("ZBosonPt"    ,"ZBosonPt"    , "", 'F'); 
    dataloader->AddVariable("ZBosonM"     ,"ZBosonM"     , "", 'F'); 
    dataloader->AddVariable("CosThetaCS"  ,"CosThetaCS"  , "", 'F'); 
    dataloader->AddVariable("CosThetaStar","CosThetaStar", "", 'F'); 
    dataloader->AddVariable("hbbpt"       ,"hbbpt"       , "", 'F'); 
    dataloader->AddVariable("hbbm"        ,"hbbm"        , "", 'F'); 
    dataloader->AddVariable("dPhiZH"      ,"dPhiZH"      , "", 'F'); 
    dataloader->AddVariable("ptBalanceZH" ,"ptBalanceZH" , "", 'F'); 
    dataloader->AddVariable("dRBjets"     ,"dRBjets"     , "", 'F'); 
    dataloader->AddVariable("dEtaBjets"   ,"dEtaBjets"   , "", 'F'); 
    //dataloader->AddVariable("dRZH"        ,"dRZH"        , "", 'F'); 
    //dataloader->AddVariable("dEtaZH"      ,"dEtaZH"      , "", 'F'); 
    dataloader->AddVariable("nAddJet"     ,"nAddJet"     , "", 'F'); 
  }
  TString prepareOptions="NormMode=None";
    prepareOptions+=":SplitMode=Block"; // use e.g. all events selected by trainTreeEventSplitStr for training
    prepareOptions+=":MixMode=Random";
  dataloader->PrepareTrainingAndTestTree("", prepareOptions);
  
  // for resolved
  TString hyperparameters=
  isBoosted?
  "!H:!V:BoostType=AdaBoost:MinNodeSize=5%:NegWeightTreatment=Pray:SeparationType=MisClassificationError:NTrees=400:MaxDepth=2:AdaBoostBeta=0.10:nCuts=10000":
  "!H:!V:BoostType=AdaBoost:MinNodeSize=5%:NegWeightTreatment=Pray:SeparationType=MisClassificationError:NTrees=200:MaxDepth=3:AdaBoostBeta=0.12:nCuts=10000";

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
