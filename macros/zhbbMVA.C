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
  bool isZHSel=true,
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
    TCut cutTrainSignal = Form("%s && (category==%d || category==%d) && MSD >= 80 && MSD < 150",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    TCut cutTrainBkg    = Form("%s && (category!=%d && category!=%d) && MSD >= 80 && MSD < 150",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    TCut cutTestSignal  = Form("%s && (category==%d || category==%d) && MSD >= 80 && MSD < 150",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    TCut cutTestBkg     = Form("%s && (category!=%d && category!=%d) && MSD >= 80 && MSD < 150",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    if     ( isBoosted &&  isZHSel) {
      // nothing to do
    }
    else if( isBoosted && !isZHSel) {
      cutTrainSignal = Form("%s && (category==%d) && MSD >= 50 && MSD < 120",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotVZbb);
      cutTrainBkg    = Form("%s && (category!=%d) && MSD >= 50 && MSD < 120",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotVZbb);
      cutTestSignal  = Form("%s && (category==%d) && MSD >= 50 && MSD < 120",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotVZbb);
      cutTestBkg     = Form("%s && (category!=%d) && MSD >= 50 && MSD < 120",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotVZbb);
    }
    else if(!isBoosted &&  isZHSel) {
      cutTrainSignal = Form("%s && (category==%d || category==%d) && hbbm >= 90 && hbbm < 150",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
      cutTrainBkg    = Form("%s && (category!=%d && category!=%d) && hbbm >= 90 && hbbm < 150",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
      cutTestSignal  = Form("%s && (category==%d || category==%d) && hbbm >= 90 && hbbm < 150",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
      cutTestBkg     = Form("%s && (category!=%d && category!=%d) && hbbm >= 90 && hbbm < 150",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotZH,(int)vhbbPlot::kPlotGGZH);
    
    }
    else if(!isBoosted && !isZHSel) {
      cutTrainSignal = Form("%s && (category==%d) && hbbm >= 60 && hbbm < 120",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotVZbb);
      cutTrainBkg    = Form("%s && (category!=%d) && hbbm >= 60 && hbbm < 120",trainTreeEventSplitStr.Data(),(int)vhbbPlot::kPlotVZbb);
      cutTestSignal  = Form("%s && (category==%d) && hbbm >= 60 && hbbm < 120",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotVZbb);
      cutTestBkg     = Form("%s && (category!=%d) && hbbm >= 60 && hbbm < 120",testTreeEventSplitStr.Data(), (int)vhbbPlot::kPlotVZbb);
    
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
    dataloader->AddVariable("nIsojet"       , "nIsojet"       , "", 'F');
    dataloader->AddVariable("fjPt"          , "fjPt"          , "", 'F');
    dataloader->AddVariable("MSD"           , "MSD"           , "", 'F');
    dataloader->AddVariable("Tau21SD"       , "Tau21SD"       , "", 'F');
    dataloader->AddVariable("ptBalanceZHFJ" , "ptBalanceZHFJ" , "", 'F');
    dataloader->AddVariable("dEtaZHFJ"      , "dEtaZHFJ"      , "", 'F');
    dataloader->AddVariable("dPhiZHFJ"      , "dPhiZHFJ"      , "", 'F');
    dataloader->AddVariable("dRL1L2"        , "dRL1L2"        , "", 'F');
    dataloader->AddVariable("lepton1Pt"     , "lepton1Pt"     , "", 'F');
    dataloader->AddVariable("lepton2Pt"     , "lepton2Pt"     , "", 'F');
    dataloader->AddVariable("CosThetaCS"    , "CosThetaCS"    , "", 'F'); 
    dataloader->AddVariable("doubleBTag"    ,  "Double B Tag" , "", 'F' ); 
    //dataloader->AddVariable("ptBalanceL1L2" , "ptBalanceL1L2" , "", 'F');
    //dataloader->AddVariable("nIsoBjet"      , "nIsoBjet"      , "", 'F');
    //dataloader->AddVariable("ZBosonPt"      , "ZBosonPt"      , "", 'F'); 
  } else {
    dataloader->AddVariable("sumEtSoft1"    , "sumEtSoft1"    , "", 'F'); 
    dataloader->AddVariable("nSoft5"        , "nSoft5"        , "", 'F'); 
    dataloader->AddVariable("bjet1Pt"       , "bjet1Pt"       , "", 'F'); 
    dataloader->AddVariable("bjet2Pt"       , "bjet2Pt"       , "", 'F'); 
    dataloader->AddVariable("bjet1btag"     , "bjet1btag"     , "", 'F'); 
    dataloader->AddVariable("bjet2btag"     , "bjet2btag"     , "", 'F'); 
    dataloader->AddVariable("lepton1Pt"     , "lepton1Pt"     , "", 'F');
    dataloader->AddVariable("lepton2Pt"     , "lepton2Pt"     , "", 'F');
    dataloader->AddVariable("ZBosonM"       , "ZBosonM"       , "", 'F'); 
    dataloader->AddVariable("CosThetaCS"    , "CosThetaCS"    , "", 'F'); 
    dataloader->AddVariable("CosThetaStar"  , "CosThetaStar"  , "", 'F'); 
    dataloader->AddVariable("hbbpt"         , "hbbpt"         , "", 'F'); 
    dataloader->AddVariable("hbbm"          , "hbbm"          , "", 'F'); 
    dataloader->AddVariable("dPhiZH"        , "dPhiZH"        , "", 'F'); 
    dataloader->AddVariable("ptBalanceZH"   , "ptBalanceZH"   , "", 'F'); 
    dataloader->AddVariable("dRZH"          , "dRZH"          , "", 'F'); 
    dataloader->AddVariable("nAddJet"       , "nAddJet"       , "", 'F'); 
    //dataloader->AddVariable("dRBjets"       , "dRBjets"       , "", 'F'); 
    //dataloader->AddVariable("ZBosonPt"      , "ZBosonPt"      , "", 'F'); 
  }
  TString prepareOptions="NormMode=None";
    prepareOptions+=":SplitMode=Block"; // use e.g. all events selected by trainTreeEventSplitStr for training
    prepareOptions+=":MixMode=Random";
  dataloader->PrepareTrainingAndTestTree("", prepareOptions);

  //Pray/IgnoreNegWeightsInTraining

  // for resolved
  TString hyperparameters=
  isBoosted?
  "!H:!V:BoostType=AdaBoost:MinNodeSize=5%:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=MisClassificationError:NTrees=500:MaxDepth=3:AdaBoostBeta=0.12:nCuts=10000":
  "!H:!V:BoostType=AdaBoost:MinNodeSize=5%:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=MisClassificationError:NTrees=500:MaxDepth=3:AdaBoostBeta=0.12:nCuts=10000";

  //TString hyperparameters="!H:!V:NTrees=500:MinNodeSize=5%:MaxDepth=3:BoostType=Grad:Shrinkage=0.1:nCuts=30:PruneMethod=CostComplexity";
  //TString hyperparameters="!H:!V:NTrees=500:NegWeightTreatment=IgnoreNegWeightsInTraining:MinNodeSize=5%:MaxDepth=2:BoostType=Grad:Shrinkage=0.1:nCuts=30";
  // for boosted
  //TString hyperparameters="!H:!V:NTrees=1000:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=MisClassificationError:MinNodeSize=5%:MaxDepth=2:BoostType=Grad:Shrinkage=0.05:nCuts=1000";
  //TString hyperparameters="!H:!V:NTrees=1000:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=MisClassificationError:MinNodeSize=5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.12:nCuts=1000";
  //if(useGaussDeco) hyperparameters += ":VarTransform=G,D";
  factory->BookMethod(dataloader, TMVA::Types::kBDT, trainName, hyperparameters);
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  output_file->Close();
}
