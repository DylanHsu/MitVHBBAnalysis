#include <TROOT.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TCut.h>
#include <TTree.h>
#include <TString.h>
#include "vhbbPlot.h"

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
  output_file=TFile::Open("/data/t3home000/dhsu/mva/zhbb/trainingResult_zllhbb_"+trainName+".root", "RECREATE");
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
    TCut cutSignal = Form("category==%d || category==%d",(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    TCut cutBkg    = Form("category!=%d && category!=%d",(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    dataloader->AddTree(mvaTree, "Background", 1.0, cutBkg   , "test train");
    dataloader->AddTree(mvaTree, "Signal"    , 1.0, cutSignal, "test train");
    dataloader->SetWeightExpression("weight", "Signal");
    dataloader->SetWeightExpression("weight", "Background");
  } 
  
  TCut preselectionCut;
  if(isBoosted) {
  } else {
    dataloader->AddVariable("sumEtSoft1"  ,"sumEtSoft1"  , "", 'F'); 
    dataloader->AddVariable("nSoft2"      ,"nSoft2"      , "", 'F'); 
    dataloader->AddVariable("nSoft5"      ,"nSoft5"      , "", 'F'); 
    dataloader->AddVariable("nSoft10"     ,"nSoft10"     , "", 'F'); 
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
    dataloader->AddVariable("dRZH"        ,"dRZH"        , "", 'F'); 
    dataloader->AddVariable("dEtaZH"      ,"dEtaZH"      , "", 'F'); 
    dataloader->AddVariable("nAddJet"     ,"nAddJet"     , "", 'F'); 
    
  }
  TString prepareOptions="NormMode=None";
    //prepareOptions+=":SplitMode=Block"; // use e.g. all events selected by trainTreeEventSplitStr for training
    prepareOptions+=":SplitMode=Random"; // use e.g. all events selected by trainTreeEventSplitStr for training
    prepareOptions+=":MixMode=Random";
  dataloader->PrepareTrainingAndTestTree("", prepareOptions);
  
  // for resolved
  TString hyperparameters="!H:!V:BoostType=AdaBoost:MinNodeSize=2.5%:NegWeightTreatment=Pray:SeparationType=MisClassificationError:NTrees=200:MaxDepth=3:AdaBoostBeta=0.12::nCuts=1000";

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
